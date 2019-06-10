#include "cramore.h"
#include "bcf_filtered_reader.h"
#include "sam_filtered_reader.h"
#include "sc_drop_seq.h"
#include "louvain.h"
#include <ctime>

///////////////////////////////////////////////////////////////////
// Freemuxlet : Genotype-free deconvolution of scRNA-seq doublets
//////////////////////////////////////////////////////////////////
int32_t cmdCramFreemux2(int32_t argc, char** argv) {
  //std::string gtfFile;
  std::string outPrefix;
  std::string plpPrefix;
  std::string initClusterFile;
  int32_t capBQ = 40;
  int32_t minBQ = 13;
  //std::vector<double> gridAlpha;
  double doublet_prior = 0.5;
  double geno_error = 0.0;
  std::string groupList;
  int32_t minTotalReads = 0;
  int32_t minUniqReads = 0;
  int32_t minCoveredSNPs = 0;
  int32_t nSamples = 0;
  double bfThres = 5.41;
  double fracInitClust = 1.00; // use 50% of cells for initial clustering
  bool auxFiles = false;
  int32_t initIteration = 10;
  bool keepInitMissing = false;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input pileup", NULL)
    LONG_STRING_PARAM("plp",&plpPrefix, "Prefix of input files generated by dsc-pileup")
    LONG_STRING_PARAM("init-cluster",&initClusterFile, "Input file containing the initial cluster information")    

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outPrefix,"Output file prefix")
//  LONG_MULTI_DOUBLE_PARAM("alpha",&gridAlpha, "Grid of alpha to search for (default is 0, 0.5)")
    LONG_INT_PARAM("nsample",&nSamples,"Number of samples multiplexed together")
    LONG_PARAM("aux-files", &auxFiles, "Turn on writing auxilary output files")
    LONG_INT_PARAM("verbose", &globalVerbosityThreshold, "Turn on verbose mode with specific verbosity threshold. 0: fully verbose, 100 : no verbose messages")

    LONG_PARAM_GROUP("Options for statistical inference", NULL)
    LONG_DOUBLE_PARAM("doublet-prior",&doublet_prior, "Prior of doublet")
    LONG_DOUBLE_PARAM("geno-error",&geno_error, "Genotype error parameter per cluster")    
    LONG_DOUBLE_PARAM("bf-thres",&bfThres,"Bayes Factor Threshold used in the initial clustering")
    LONG_DOUBLE_PARAM("frac-init-clust",&fracInitClust,"Fraction of droplets to be clustered in the very first round of initial clustering procedure")
    LONG_INT_PARAM("iter-init",&initIteration, "Iteration for initial cluster assignment (set to zero to skip the iterations)")
    LONG_PARAM("keep-init-missing",&keepInitMissing, "Keep missing cluster assignment as missing in the initial iteration")    
    
    LONG_PARAM_GROUP("Read filtering Options", NULL)
    LONG_INT_PARAM("cap-BQ", &capBQ, "Maximum base quality (higher BQ will be capped)")
    LONG_INT_PARAM("min-BQ", &minBQ, "Minimum base quality to consider (lower BQ will be skipped)")

    LONG_PARAM_GROUP("Cell/droplet filtering options", NULL)
    LONG_STRING_PARAM("group-list",&groupList, "List of tag readgroup/cell barcode to consider in this run. All other barcodes will be ignored. This is useful for parallelized run")    
    LONG_INT_PARAM("min-total", &minTotalReads, "Minimum number of total reads for a droplet/cell to be considered")
    LONG_INT_PARAM("min-uniq", &minUniqReads, "Minimum number of unique reads (determined by UMI/SNP pair) for a droplet/cell to be considered")
    LONG_INT_PARAM("min-snp", &minCoveredSNPs, "Minimum number of SNPs with coverage for a droplet/cell to be considered")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( plpPrefix.empty() || outPrefix.empty() || ( nSamples == 0 ) )
    error("Missing required option(s) : --plp, --out, --nsample");

  std::set<std::string> bcdSet;
  sc_dropseq_lib_t scl;

  scl.load_from_plp(plpPrefix.c_str());

  std::map<std::string, int32_t> initCluster;

  // if initial clusters are provided, use them here
  if ( !initClusterFile.empty() ) {
    tsv_reader tsv_clustf(initClusterFile.c_str());
    while ( tsv_clustf.read_line() > 0 ) {
      if ( tsv_clustf.nfields != 2 )
	error("ERROR: Initial clustering file %d has to have 2 columnes", initClusterFile.c_str());
      int32_t iclust = tsv_clustf.int_field_at(1);
      if ( iclust >= 0 ) {
	if ( iclust >= nSamples )
	  error("ERROR: --nsample %d parameter was set. The cluster ID must be between 0 to %d, or use negative values to not assign initial cluster (not implemented yet)", nSamples, nSamples-1);
	initCluster[tsv_clustf.str_field_at(0)] = iclust;
      }
    }
  }

  // sort cells based on the number of SNP-overlapping unique reads using singlet scores
  htsFile* wmix = NULL;
  std::vector<int32_t> nSNPs(scl.nbcs,0);
  std::vector<int32_t> nReads(scl.nbcs,0);  

  wmix = hts_open((outPrefix+".lmix").c_str(),"w");
  hprintf(wmix, "INT_ID\tBARCODE\tNSNPs\tNREADs\tDBL.LLK\tSNG.LLK\tLOG.BF\tBFpSNP\n");
  
  std::vector< std::map<int32_t,snp_droplet_pileup*> > cell_snp_plps(scl.nbcs);
  std::vector< std::map<int32_t,snp_droplet_pileup*> > snp_cell_plps(scl.nsnps);
    
  for(int32_t i=0; i < scl.nbcs; ++i) {
    int32_t si = i; // drops_srted[i];
    if (i % 1000 == 0 )
      notice("Processing singlet scores for %d droplets..", i+1);

    // likelihood calculation across the overlapping SNPs
    std::map<int32_t,sc_snp_droplet_t* >::iterator it = scl.cell_umis[si].begin();

    double llk0 = 0, llk2 = 0; // LLK of IBD0, IBD1, IBD2     
    
    while( it != scl.cell_umis[si].end() ) {
      //double gls[9] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
      double af = scl.snps[it->first].af;

      // calculate genotype likelihoods
      // calculate_snp_droplet_doublet_GL(it->second, gls, 0.5);
      if ( cell_snp_plps[i][it->first] == NULL ) 
	cell_snp_plps[i][it->first] = snp_cell_plps[it->first][i] = new snp_droplet_pileup;
      calculate_snp_droplet_pileup(it->second, cell_snp_plps[i][it->first], 0.5);
      double* gls = cell_snp_plps[i][it->first]->gls;

      double lk0 = 0, lk2 = 0;
      double gps[3];
      gps[0] = (1.0-af) * (1.0-af);
      gps[1] = 2.0 * af * (1.0-af);
      gps[2] = af * af;
	
      for(int32_t gi=0; gi < 3; ++gi) {
	lk2 += ( gls[gi*3 + gi] * gps[gi] );
	for(int32_t gj=0; gj < 3; ++gj) {
	  lk0 += ( gls[gi*3 + gj] * gps[gi] * gps[gj] );
	}
      }
      nReads[i] += (int32_t)it->second->size();
      ++nSNPs[i];
      
      ++it;

      llk0 += log(lk0);
      llk2 += log(lk2);
    }

    scl.cell_scores[si] = llk2 - llk0;  // score of being singlet

    hprintf(wmix,"%d\t%s\t%d\t%d\t%.2lf\t%.2lf\t%.2lf\t%.4lf\n", si, scl.bcs[si].c_str(), nSNPs[i], nReads[i], llk0, llk2, llk0-llk2, (llk0-llk2)/nSNPs[i]);
  }
  hts_close(wmix);

  // sort droplets by singlet scores
  std::vector<int32_t> drops_srted(scl.nbcs);
  for(int32_t i=0; i < scl.nbcs; ++i) {
    drops_srted[i] = i;
  }  
  sc_drop_comp_t sdc(&scl);
  std::sort( drops_srted.begin(), drops_srted.end(), sdc );


  std::vector<int32_t> clusts(scl.nbcs,-1);
  std::vector<int32_t> ccounts(nSamples,0);  
  // initial clustering
  // calculate pairwise genetic distances on demand while clustering
  // use the assigned clusters if already provided
  if ( !initClusterFile.empty() ) {
    int32_t nmiss = 0;
    for(int32_t i=0; i < scl.nbcs; ++i) {
      std::map<std::string,int32_t>::iterator it =
	initCluster.find(scl.bcs[i]);
      if ( it == initCluster.end() ) {
	++nmiss;
	//error("ERROR: Cannot find the initial cluster for barcode %s", scl.bcs[i].c_str());
      }
      else {
	clusts[i] = it->second;
	++ccounts[it->second];
      }
    }
    if ( nmiss > 0 ) {
      warning("WARNING: %d of %d droplets do not have initial cluster assignment", nmiss, scl.nbcs);
    }    
  }
  else { // greedy clustering
    // maintains GLs for each cluster, initially, everything is set up to be empty for each variant.
    std::vector< std::map<int32_t,snp_droplet_pileup> > clustPileup(nSamples);

    // greedy initial clustering
    double sumMaxScore = 0;
    for(int32_t i=0; i < scl.nbcs; ++i) {
      int32_t si = drops_srted[i];
      if ( i > scl.nbcs * fracInitClust ) continue; // skip the droplet if exceed initial fraction of cells to be clustered.

      // compute the distance with each clusters
      std::vector<dropD> dropDs;
      for(int32_t j=0; j < nSamples; ++j) {
	// compute genetic distance with each droplet
	dropDs.push_back(scl.calculate_droplet_clust_distance(cell_snp_plps[si], clustPileup[j]));
      }

      int32_t maxClust = 0;
      double maxScore = dropDs[0].llk2 - dropDs[0].llk0;
      for(int32_t j=1; j < nSamples; ++j) {
	if ( dropDs[j].llk2 - dropDs[j].llk0 > maxScore ) {
	  maxClust = j;
	  maxScore = dropDs[j].llk2 - dropDs[j].llk0;
	}
      }
      clusts[si] = maxClust;
      ++ccounts[maxClust];
      sumMaxScore += maxScore;

      for(std::map<int32_t,snp_droplet_pileup*>::const_iterator it = cell_snp_plps[si].begin();
	  it != cell_snp_plps[si].end(); ++it) {
	clustPileup[clusts[si]][it->first].merge(*it->second);
      }

      if ( i % 100 == 0 ) {
	std::string s;
	catprintf(s, "Processing %d droplets. Avg maxScore = %.5lg. Cluster counts:", i+1, sumMaxScore/(i+1));
	for(int32_t j=0; j < nSamples; ++j) 
	  catprintf(s, " %d",ccounts[j]);
	notice(s.c_str());
      }
    }
  }

  notice("Finished assigning initial identity of the cluster..");

  if ( auxFiles ) {
    htsFile* wc0 = hts_open((outPrefix+".clust0.samples.gz").c_str(),"wz");
    hprintf(wc0, "INT_ID\tBARCODE\tCLUST0\n");
    //std::vector< std::vector<int32_t> > iclusts(nSamples);
    for(int32_t i=0; i < scl.nbcs; ++i) {
      hprintf(wc0, "%d\t%s\t%d\n", i, scl.bcs[i].c_str(), clusts[i]);
      //iclusts[clusts[i]].push_back(i);
    }
    hts_close(wc0);
  }

  std::vector< std::map<int32_t,snp_droplet_pileup> > clustPileup(nSamples);

  std::vector<bool> snps_observed(scl.nsnps,false);
  for(int32_t i=0; i < scl.nbcs; ++i) {
    std::map<int32_t,snp_droplet_pileup*>::const_iterator it = cell_snp_plps[i].begin();
    while(it != cell_snp_plps[i].end()) {
      if ( clusts[i] >= 0 ) 
	clustPileup[clusts[i]][it->first].merge(*it->second);
      snps_observed[it->first] = true;
      ++it;
    }
  }
  
  time_t now = std::time(NULL);
  tm *ltm = localtime(&now);
  
  // write initial clusters
  if ( auxFiles ) {
    htsFile* vc0 = hts_open((outPrefix+".clust0.vcf.gz").c_str(),"wz");
    hprintf(vc0,"##fileformat=VCFv4.2\n");
    hprintf(vc0,"##fileDate=%04d%02d%02d\n",1970+ltm->tm_year,1+ltm->tm_mon,ltm->tm_mday);
    hprintf(vc0,"##source=cramore-freemuxlet\n");
    for(int32_t i=0; i < (int32_t)scl.rid2chr.size(); ++i)
      hprintf(vc0, "##contig=<ID=%s>\n", scl.rid2chr[i].c_str());
    hprintf(vc0,"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n");
    hprintf(vc0,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    hprintf(vc0,"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Phred-scale Genotype Quality\">\n");    
    hprintf(vc0,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
    
    hprintf(vc0,"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic Read Depth\">\n");
    hprintf(vc0,"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scale genotype likelihood\">\n");
    hprintf(vc0,"##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Posterior probability using pooled allele frequencies\">\n");              
    hprintf(vc0,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for(int32_t i=0; i < nSamples; ++i) hprintf(vc0, "\tCLUST%d", i);
    hprintf(vc0, "\n");
    for(int32_t v=0; v < scl.nsnps; ++v) {
      if ( !snps_observed[v] ) continue;
      sc_snp_t& s = scl.snps[v];
      hprintf(vc0,"%s\t%d\t.\t%c\t%c\t.\tPASS\tAF=%.5lf\tGT:GQ:DP:AD:PL:GP",scl.rid2chr[s.rid].c_str(),s.pos,s.ref,s.alt,s.af);
      double gps[3] = { (1-s.af)*(1-s.af), 2.*s.af*(1-s.af), s.af*s.af };
      double pps[3], sumPP;
      int32_t pls[3];
      int bestG, gq;
      for(int32_t i=0; i < nSamples; ++i) {
	snp_droplet_pileup& sdp = clustPileup[i][v];
	double maxGL = sdp.gls[0];
	if ( maxGL < sdp.gls[4] ) maxGL = sdp.gls[4];
	if ( maxGL < sdp.gls[8] ) maxGL = sdp.gls[8];
	pls[0] = (int)(-10.0*log10(sdp.gls[0]/maxGL));
	pls[1] = (int)(-10.0*log10(sdp.gls[4]/maxGL));
	pls[2] = (int)(-10.0*log10(sdp.gls[8]/maxGL));

	pps[0] = gps[0] * sdp.gls[0] / maxGL + 1e-100;
	pps[1] = gps[1] * sdp.gls[4] / maxGL + 1e-100;
	pps[2] = gps[2] * sdp.gls[8] / maxGL + 1e-100;
	sumPP = pps[0] + pps[1] + pps[2];

	pps[0] /= sumPP;
	pps[1] /= sumPP;
	pps[2] /= sumPP;

	bestG = ( pps[0] > pps[1] ) ? ( pps[0] > pps[2] ? 0 : 2 ) : ( pps[1] > pps[2] ? 1 : 2 );
	gq = (int)(-0.1*log10(1-pps[bestG]+1e-100));
	if ( gq > 255 ) gq = 255;
	//hprintf(vc0,"\t%d/%d:%d:%d:%d,%d:%d,%d,%d",bestG == 2 ? 1 : 0, bestG > 0 ? 1 : 0, gq, sdp.nreads,sdp.nref,sdp.nalt,pls[0],pls[1],pls[2]);
	hprintf(vc0,"\t%d/%d:%d:%d:%d,%d:%d,%d,%d:%.3lg,%.3lg,%.3lg",bestG == 2 ? 1 : 0, bestG > 0 ? 1 : 0, gq, sdp.nreads,sdp.nref,sdp.nalt,pls[0], pls[1], pls[2], pps[0], pps[1], pps[2]);	
      }
      hprintf(vc0,"\n");
    }
    hts_close(vc0);
  }


  std::vector<int32_t> jBests(scl.nbcs,-1);
  std::vector<int32_t> kBests(scl.nbcs,-1);
  std::vector<int32_t> jNexts(scl.nbcs,-1);
  std::vector<int32_t> kNexts(scl.nbcs,-1);  
  std::vector<double> bestLLKs(scl.nbcs,-1e300);
  std::vector<double> nextLLKs(scl.nbcs,-1e300);
  std::vector<double> sngBestLLKs(scl.nbcs,-1e300);
  std::vector<double> sngNextLLKs(scl.nbcs,-1e300);
  std::vector<int32_t> sBests(scl.nbcs,-1);
  std::vector<int32_t> sNexts(scl.nbcs,-1);    
  std::vector<int32_t> dBest1s(scl.nbcs,-1);
  std::vector<int32_t> dBest2s(scl.nbcs,-1);
  std::vector<int32_t> dNext1s(scl.nbcs,-1);
  std::vector<int32_t> dNext2s(scl.nbcs,-1);      
  std::vector<double> dblBestLLKs(scl.nbcs,-1e300);
  std::vector<double> dblNextLLKs(scl.nbcs,-1e300);
  // posterior probs
  std::vector<double> bestPPs(scl.nbcs,-1e300);
  std::vector<double> sngPPs(scl.nbcs,-1e300);
  std::vector<double> sngOnlyPPs(scl.nbcs,-1e300);
  std::vector<double> sumLLKs(scl.nbcs,-1e300);  
  std::vector<int32_t> types(scl.nbcs,-1);
      
  // calculate probabilities of singlets/doublets
  int32_t max_iter = 10;  
  for(int32_t iter=0; iter < max_iter; ++iter) {
    notice("Inferring doublets and refining clusters.., iter = %d", iter+1);
    
    double gp1s[3], gp2s[3], gp0s[3], sum1, sum2;
    int32_t npairs = nSamples*(nSamples+1)/2;
    double log_single_prior = log((1.0-doublet_prior)/nSamples);
    double log_double_prior = log(doublet_prior/nSamples/(nSamples-1)*2.0);

    for(int32_t i=0; i < scl.nbcs; ++i) {
      std::vector<double> llks(npairs, 0);
      std::map<int32_t,snp_droplet_pileup*>::iterator it;
      for(it = cell_snp_plps[i].begin(); it != cell_snp_plps[i].end(); ++it) {
	double af = scl.snps[it->first].af;
	gp0s[0] = (1.0-af)*(1.0-af);
	gp0s[1] = 2*af*(1.0-af);
	gp0s[2] = af*af;
	std::vector<double> lks(npairs, 0);
	double lk;
	double* glis = it->second->gls;	
	for(int32_t j=0; j < nSamples; ++j) {
	  snp_droplet_pileup& sdp1 = clustPileup[j][it->first];
	  gp1s[0] = (1.0-af)*(1.0-af)*sdp1.gls[0];
	  gp1s[1] = 2*af*(1.0-af)*sdp1.gls[4];
	  gp1s[2] = af*af*sdp1.gls[8];
	  sum1 = gp1s[0]+gp1s[1]+gp1s[2];
	  gp1s[0] /= sum1;
	  gp1s[1] /= sum1;
	  gp1s[2] /= sum1;
	  if ( geno_error > 0 ) {
	  //if ( ( geno_error > 0 ) && ( iter + 1 == max_iter ) ) {	    
	    gp1s[0] = (1-geno_error)*gp1s[0] + geno_error*gp0s[0];
	    gp1s[1] = (1-geno_error)*gp1s[1] + geno_error*gp0s[1];
	    gp1s[2] = (1-geno_error)*gp1s[2] + geno_error*gp0s[2];	    
	  }
	  for(int32_t k=0; k < j; ++k) {
	    snp_droplet_pileup& sdp2 = clustPileup[k][it->first];	  
	    // Pr(D|g1,g2)Pr(g1|C1)Pr(g2|C2)Pr(C1)Pr(C2)
	    gp2s[0] = (1.0-af)*(1.0-af)*sdp2.gls[0];
	    gp2s[1] = 2*af*(1.0-af)*sdp2.gls[4];
	    gp2s[2] = af*af*sdp2.gls[8];
	    sum2 = gp2s[0]+gp2s[1]+gp2s[2];
	    gp2s[0] /= sum2;
	    gp2s[1] /= sum2;
	    gp2s[2] /= sum2;
	    if ( geno_error > 0 ) {
	    //if ( ( geno_error > 0 ) && ( iter + 1 == max_iter ) ) {	      
	      gp2s[0] = (1-geno_error)*gp2s[0] + geno_error*gp0s[0];
	      gp2s[1] = (1-geno_error)*gp2s[1] + geno_error*gp0s[1];
	      gp2s[2] = (1-geno_error)*gp2s[2] + geno_error*gp0s[2];	    
	    }	    
	    
	    lk = 0;
	    for(int32_t g1=0; g1 < 3; ++g1) {
	      for(int32_t g2=0; g2 < 3; ++g2) {
		lk += ( glis[g1*3+g2] * gp1s[g1] * gp2s[g2] );
	      }
	    }
	    lks[j*(j+1)/2+k] = lk;
	  }
	  lk = 0;
	  for(int32_t g1=0; g1 < 3; ++g1) {
	    lk += ( glis[g1*3+g1] * gp1s[g1] );
	  }
	  lks[j*(j+1)/2+j] = lk;
	}
	for(int32_t i=0; i < npairs; ++i)
	  llks[i] += log(lks[i]);
      }

      //int32_t jBest = -1, kBest = -1, jNext = -1, kNext = -1;
      int32_t sBest = -1, sNext = -1, dBest1 = -1, dBest2 = -1, dNext1 = -1, dNext2 = -1;
      //double bestLLK = -1e300;
      //double nextLLK = -1e300;      
      double sngBestLLK = -1e300;
      double sngNextLLK = -1e300;
      double dblBestLLK = -1e300;
      double dblNextLLK = -1e300;            
      double sumLLK = -1e300;
      double sngLLK = -1e300;
      double tmpLLK;
      for(int32_t j=0; j < nSamples; ++j) {
	for(int32_t k=0; k < j; ++k) {
	  tmpLLK = llks[j*(j+1)/2+k]; // + log_double_prior;
	  if ( tmpLLK > dblBestLLK ) {
	    dNext1 = dBest1; dNext2 = dBest2;
	    dblNextLLK = dblBestLLK;
	    dBest1 = j; dBest2 = k;
	    dblBestLLK = tmpLLK;
	  }
	  else if ( tmpLLK > dblNextLLK ) {
	    dNext1 = j; dNext2 = k;
	    dblNextLLK = tmpLLK;
	  }
	  sumLLK = logAdd(sumLLK,tmpLLK + log_double_prior);
	}
	
	tmpLLK = llks[j*(j+1)/2+j]; //+ log_single_prior;
	if ( tmpLLK > sngBestLLK ) {
	  sNext = sBest;
	  sngNextLLK = sngBestLLK;
	  sBest = j;
	  sngBestLLK = tmpLLK;
	}
	else if ( tmpLLK > sngNextLLK ) {
	  sNext = j;
	  sngNextLLK = tmpLLK;
	}	
	sumLLK = logAdd(sumLLK,tmpLLK + log_single_prior);
	sngLLK = logAdd(sngLLK,tmpLLK + log_single_prior);	
      }

      sBests[i] = sBest;
      sngBestLLKs[i] = sngBestLLK;
      sNexts[i] = sNext;
      sngNextLLKs[i] = sngNextLLK;
      dBest1s[i] = dBest1;
      dBest2s[i] = dBest2;      
      dblBestLLKs[i] = dblBestLLK;
      dNext1s[i] = dNext1;
      dNext2s[i] = dNext2;      
      dblNextLLKs[i] = dblNextLLK;
      sngPPs[i]  = exp(sngLLK - sumLLK);
      sngOnlyPPs[i] = exp(sngBestLLK + log_single_prior - sngLLK);
      sumLLKs[i] = sumLLK;
    }

    // re-assign sample identities
    clustPileup.clear();
    clustPileup.resize(nSamples);
    int32_t nsingle = 0, namb = 0, nchanged = 0;
    for(int32_t i=0; i < scl.nbcs; ++i) {
      if ( dblBestLLKs[i] > sngBestLLKs[i] + 2 ) { // best call is doublet
	// consider as changed only when the assignment category was changed.
	if ( types[i] != 1 ) ++nchanged;
	
	types[i] = 1; // doublet
	bestPPs[i] = ( dblBestLLKs[i] + log_double_prior - sumLLKs[i] );

	jBests[i] = dBest1s[i];
	kBests[i] = dBest2s[i];
	bestLLKs[i] = dblBestLLKs[i];

	if ( dblNextLLKs[i] > sngBestLLKs[i] + 2 ) { // next best is doublet
	  jNexts[i] = dNext1s[i];
	  kNexts[i] = dNext2s[i];
	  nextLLKs[i] = dblNextLLKs[i];
	}
	else {
	  jNexts[i] = kNexts[i] = sBests[i];  // next best is singlet
	  nextLLKs[i] = sngBestLLKs[i];
	}
      }
      else if ( sngBestLLKs[i] > sngNextLLKs[i] + 2 ) { // double call is singlet
	if ( ( types[i] != 0 ) || ( jBests[i] != sBests[i] ) || ( kBests[i] != sBests[i] ) )
	  ++nchanged;
	
	types[i] = 0; // singlet
	++nsingle;

	bestPPs[i] = ( sngBestLLKs[i] + log_single_prior - sumLLKs[i] );
	jBests[i] = kBests[i] = sBests[i];
	bestLLKs[i] = sngBestLLKs[i];

	if ( dblBestLLKs[i] > sngNextLLKs[i] + 2 ) { // next best is doublet
	  jNexts[i] = dBest1s[i];
	  kNexts[i] = dBest2s[i];
	  nextLLKs[i] = dblBestLLKs[i];	  
	}
	else {
	  jNexts[i] = kNexts[i] = sNexts[i];  // next best is also singlet
	  nextLLKs[i] = sngNextLLKs[i];	  
	}
      }
      else {  // ambiguous calls, use singlet as the best call
	if ( types[i] != 2 ) ++nchanged;
	
	types[i] = 2; // ambiguous
	++namb;

	bestPPs[i] = ( sngBestLLKs[i] + log_single_prior - sumLLKs[i] );
	jBests[i] = kBests[i] = sBests[i];
	bestLLKs[i] = sngBestLLKs[i];

	if ( dblBestLLKs[i] > sngNextLLKs[i] + 2 ) {
	  jNexts[i] = dBest1s[i];
	  kNexts[i] = dBest2s[i];
	  nextLLKs[i] = dblNextLLKs[i];	  
	}
	else {
	  jNexts[i] = kNexts[i] = sNexts[i];
	  nextLLKs[i] = sngNextLLKs[i];	  
	}	
      }
      
      // old criteria
      //if ( bestPPs[i] < 0.8 ) ++namb;
      //else if ( jBests[i] == kBests[i] ) ++nsingle;

      std::map<int32_t,snp_droplet_pileup*>::const_iterator it = cell_snp_plps[i].begin();
      while(it != cell_snp_plps[i].end()) {
	if ( ( jBests[i] == kBests[i] ) && ( types[i] == 0 ) ) {
	  clustPileup[jBests[i]][it->first].merge(*it->second);
	}
	++it;
      }      
    }

    notice("Refining per-cluster genotype likelihoods.... %d singlets, %d doublets, %d ambiguous, and %d changed", nsingle, scl.nbcs-nsingle-namb, namb, nchanged);

    if ( nchanged == 0 ) {
      notice("No more changes in cluster assginment and singlet identities. Finishing iterations early");
      break;
    }
  }


  htsFile* vc1 = hts_open((outPrefix+".clust1.vcf.gz").c_str(),"wz");
  hprintf(vc1,"##fileformat=VCFv4.2\n");
  hprintf(vc1,"##fileDate=%04d%02d%02d\n",1970+ltm->tm_year,1+ltm->tm_mon,ltm->tm_mday);
  hprintf(vc1,"##source=cramore-freemuxlet\n");
  for(int32_t i=0; i < (int32_t)scl.rid2chr.size(); ++i)
    hprintf(vc1, "##contig=<ID=%s>\n", scl.rid2chr[i].c_str());
  hprintf(vc1,"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n");
  hprintf(vc1,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
  hprintf(vc1,"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Phred-scale Genotype Quality\">\n");    
  hprintf(vc1,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
  hprintf(vc1,"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic Read Depth\">\n");
  hprintf(vc1,"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scale genotype likelihood\">\n");
  hprintf(vc1,"##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Posterior probability using pooled allele frequencies\">\n");          
  hprintf(vc1,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");  
  for(int32_t i=0; i < nSamples; ++i) hprintf(vc1, "\tCLUST%d", i);
  hprintf(vc1, "\n");
  for(int32_t v=0; v < scl.nsnps; ++v) {
    if ( !snps_observed[v] ) continue;
    sc_snp_t& s = scl.snps[v];
    hprintf(vc1,"%s\t%d\t.\t%c\t%c\t.\tPASS\tAF=%.5lf\tGT:GQ:DP:AD:PL:GP",scl.rid2chr[s.rid].c_str(),s.pos,s.ref,s.alt,s.af);
    double gps[3] = { (1.-s.af)*(1.-s.af), 2. * s.af* (1.-s.af), s.af * s.af };
    double pps[3], sumPP;
    int32_t pls[3];    
    int bestG, gq;

    for(int32_t i=0; i < nSamples; ++i) {
      snp_droplet_pileup& sdp = clustPileup[i][v];
      double maxGL = sdp.gls[0];
      if ( maxGL < sdp.gls[4] ) maxGL = sdp.gls[4];
      if ( maxGL < sdp.gls[8] ) maxGL = sdp.gls[8];
      pls[0] = (int32_t)(-10.0*log10(sdp.gls[0]/maxGL));
      pls[1] = (int32_t)(-10.0*log10(sdp.gls[4]/maxGL));
      pls[2] = (int32_t)(-10.0*log10(sdp.gls[8]/maxGL));
      
      pps[0] = gps[0] * ( sdp.gls[0] / maxGL ) + 1e-100;
      pps[1] = gps[1] * ( sdp.gls[4] / maxGL ) + 1e-100;
      pps[2] = gps[2] * ( sdp.gls[8] / maxGL ) + 1e-100;
      sumPP = pps[0] + pps[1] + pps[2];
      
      pps[0] /= sumPP;
      pps[1] /= sumPP;
      pps[2] /= sumPP;
      
      bestG = ( pps[0] > pps[1] ) ? ( pps[0] > pps[2] ? 0 : 2 ) : ( pps[1] > pps[2] ? 1 : 2 );
      gq = (int32_t)(-10*log10(1.0-pps[bestG]+1e-100));
      if ( gq > 255 ) gq = 255;
      hprintf(vc1,"\t%d/%d:%d:%d:%d,%d:%d,%d,%d:%.3lg,%.3lg,%.3lg", bestG == 2 ? 1 : 0, bestG > 0 ? 1 : 0, gq, sdp.nreads, sdp.nref, sdp.nalt, pls[0], pls[1], pls[2], pps[0], pps[1], pps[2]);
    }
    hprintf(vc1,"\n");
  }
  hts_close(vc1);

  htsFile* wc1 = hts_open((outPrefix+".clust1.samples.gz").c_str(),"wz");
  hprintf(wc1, "INT_ID\tBARCODE\tNUM.SNPS\tNUM.READS\tDROPLET.TYPE\tBEST.GUESS\tBEST.LLK\tNEXT.GUESS\tNEXT.LLK\tDIFF.LLK.BEST.NEXT\tBEST.POSTERIOR\tSNG.POSTERIOR\tSNG.BEST.GUESS\tSNG.BEST.LLK\tSNG.NEXT.GUESS\tSNG.NEXT.LLK\tSNG.ONLY.POSTERIOR\tDBL.BEST.GUESS\tDBL.BEST.LLK\tDIFF.LLK.SNG.DBL\n");
  for(int32_t i=0; i < scl.nbcs; ++i) {
    hprintf(wc1, "%d\t%s\t%d\t%d\t%s\t%d,%d\t%.2lf\t%d,%d\t%.2lf\t%.2lf\t%.5lf\t%.2lg\t%d\t%.2lf\t%d\t%.2lf\t%.5lf\t%d,%d\t%.2lf\t%.2lf\n", i, scl.bcs[i].c_str(), nSNPs[i], nReads[i], (types[i] == 2) ? "AMB" : ((types[i] == 0) ? "SNG" : "DBL"), jBests[i], kBests[i], bestLLKs[i], jNexts[i], kNexts[i], nextLLKs[i], bestLLKs[i]-nextLLKs[i], bestPPs[i], sngPPs[i], sBests[i], sngBestLLKs[i], sNexts[i], sngNextLLKs[i], sngOnlyPPs[i], dBest1s[i], dBest2s[i], dblBestLLKs[i], sngBestLLKs[i]-dblBestLLKs[i]);
  }
  hts_close(wc1);  

  return 0;
}
