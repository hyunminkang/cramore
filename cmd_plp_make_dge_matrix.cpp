#include "cramore.h"
#include "gtf.h"
#include "tsv_reader.h"
#include <ctime>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>

///////////////////////////////////////////////////////////////////
// plp-make-dge-matrix : Make DGE matrix from digital pileup file 
//////////////////////////////////////////////////////////////////
int32_t cmdPlpMakeDGEMatrix(int32_t argc, char** argv) {
  std::string gtfFile;
  std::string plpPrefix;  
  std::string outPrefix;
  std::string groupList;
  int32_t minTotalReads = 0;
  int32_t minUniqReads = 0;
  int32_t minCoveredSNPs = 0;
  bool removeChrPrefix = false;
  bool addChrPrefix = false;
  std::vector<std::string> genetypes;
  bool commonGenetypes = false;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("plp",&plpPrefix, "Prefix of input files generated by dsc-pileup")
    LONG_STRING_PARAM("gtf",&gtfFile, "GTF-formatted file for gene/transcript annotation")
    LONG_PARAM("gtf-remove-chr",&removeChrPrefix, "Remove 'chr' prefix from input GTF file")
    LONG_PARAM("gtf-add-chr"   ,&addChrPrefix,    "Add 'chr' prefix from input GTF file")
    LONG_MULTI_STRING_PARAM("gene-type", &genetypes, "Gene types to include to produce DGE matrix (e.g. protein-coding)")
    LONG_PARAM("common-gene-types", &commonGenetypes, "Load only common gene types, searching for specific gene types - protein_coding, lincRNA, antisense, IG_ and TR_ genes")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outPrefix,"Output file prefix")
    LONG_INT_PARAM("verbose", &globalVerbosityThreshold, "Turn on verbose mode with specific verbosity threshold. 0: fully verbose, 100 : no verbose messages")

    LONG_PARAM_GROUP("Cell/droplet filtering options", NULL)
    LONG_STRING_PARAM("group-list",&groupList, "List of tag readgroup/cell barcode to consider in this run. All other barcodes will be ignored. This is useful for parallelized run")
    LONG_INT_PARAM("min-total", &minTotalReads, "Minimum number of total reads for a droplet/cell to be considered")
    LONG_INT_PARAM("min-uniq",  &minUniqReads,  "Minimum number of unique reads (determined by UMI/SNP pair) for a droplet/cell to be considered")
    LONG_INT_PARAM("min-snp",   &minCoveredSNPs,"Minimum number of SNPs with coverage for a droplet/cell to be considered")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( plpPrefix.empty() || outPrefix.empty() || gtfFile.empty() )
    error("Missing required option(s) : --plp, --out, --gtf");

  // write out directories every possible elements
  int32_t ret = mkdir(outPrefix.c_str(), 0777); // create output directory
  if ( ret )
    error("FATAL ERROR -- Cannot create a directory %s", outPrefix.c_str());

  notice("Analysis Started");    

  // read pileups
  char fname[65535];
  sprintf(fname, "%s.cel.gz", plpPrefix.c_str());
  notice("Reading barcode information from %s.cel.gz..", plpPrefix.c_str());
  tsv_reader tsv_bcdf(fname);

  // read the header line
  int32_t n_expected_toks = 6;
  if ( tsv_bcdf.read_line() > 0 ) {
    if ( ( tsv_bcdf.nfields == 5 ) &&   // for backward compatibility
	 ( ( strcmp("#DROPLET_ID",tsv_bcdf.str_field_at(0)) != 0 ) ||
	   ( strcmp("BARCODE",tsv_bcdf.str_field_at(1)) != 0 ) ||
	   ( strcmp("NUM.READ",tsv_bcdf.str_field_at(2)) != 0 ) ||
	   ( strcmp("NUM.UMI",tsv_bcdf.str_field_at(3)) != 0 ) ||
	   ( strcmp("NUM.SNP",tsv_bcdf.str_field_at(4)) != 0 ) ) ) {
      error("The header line of %s.cel.gz is malformed or outdated. Expecting #DROPLET_ID BARCODE NUM.READ NUM.UMI NUM.SNP", plpPrefix);
    }
    else if ( ( tsv_bcdf.nfields == 6 ) &&
	 ( ( strcmp("#DROPLET_ID",tsv_bcdf.str_field_at(0)) != 0 ) ||
	   ( strcmp("BARCODE",tsv_bcdf.str_field_at(1)) != 0 ) ||
	   ( strcmp("NUM.READ",tsv_bcdf.str_field_at(2)) != 0 ) ||
	   ( strcmp("NUM.UMI",tsv_bcdf.str_field_at(3)) != 0 ) ||
	   ( strcmp("NUM.UMIwSNP",tsv_bcdf.str_field_at(4)) != 0 ) ||	   
	   ( strcmp("NUM.SNP",tsv_bcdf.str_field_at(5)) != 0 ) ) ) {
      error("The header line of %s.cel.gz is malformed or outdated. Expecting #DROPLET_ID BARCODE NUM.READ NUM.UMI NUM.UMIwSNP NUM.SNP", plpPrefix);
    }
    else if ( ( tsv_bcdf.nfields < 5 ) || ( tsv_bcdf.nfields > 6 ) ) {
      error("The header line of %s.cel.gz is malformed or outdated. Expecting #DROPLET_ID BARCODE NUM.READ NUM.UMI (NUM.UMIwSNP-optional) NUM.SNP", plpPrefix);  
    }
    n_expected_toks = tsv_bcdf.nfields;
  }
  else error("Cannot read the first line of %s.cel.gz", plpPrefix);

  // read the barcode information
  std::map<int32_t,int32_t> id_cel2dge;
  std::vector<std::string> dgebcds;
  std::vector<int32_t> cell_totl_reads;
  std::vector<int32_t> cell_totl_umis;  
  std::vector<int32_t> cell_umi_w_snps;
  std::vector<int32_t> cell_num_snps;      

  int32_t ibcd = 0;
  int32_t skipbcd = 0;
  for( int32_t i=0; tsv_bcdf.read_line() > 0; ++i ) {
    if ( i != tsv_bcdf.int_field_at(0) ) 
      error("[E:%s] Observed DROPLET_ID %d is different from expected DROPLET_ID %d. Did you modify the digital pileup files by yourself?", __PRETTY_FUNCTION__, tsv_bcdf.int_field_at(0), i);
    int32_t numisnp = 0, numi = 0, nsnp = 0;
    int32_t nread = numisnp = tsv_bcdf.int_field_at(2);          
    if ( n_expected_toks == 5 ) {
      numi = numisnp = tsv_bcdf.int_field_at(3);
      nsnp = tsv_bcdf.int_field_at(4);
    }
    else {
      numi = tsv_bcdf.int_field_at(3);
      numisnp = tsv_bcdf.int_field_at(4);
      nsnp = tsv_bcdf.int_field_at(5);
    }
    if ( ( nread < minTotalReads ) || ( numi < minUniqReads ) || ( nsnp < minCoveredSNPs )  ) {
      ++skipbcd;
      continue;
    }
    id_cel2dge[i] = ibcd;
    dgebcds.push_back(tsv_bcdf.str_field_at(1));
    cell_totl_reads.push_back(nread);
    cell_totl_umis.push_back(numi);
    cell_umi_w_snps.push_back(numisnp);
    cell_num_snps.push_back(nsnp);
    ++ibcd;
  }
  notice("Finished loading %d droplets, skipping %d.", ibcd, skipbcd);
  //tsv_bcdf.close();

  if ( commonGenetypes ) {
    genetypes.push_back("protein_coding");
    genetypes.push_back("lincRNA");
    genetypes.push_back("antisense");
    genetypes.push_back("IG_LV_gene");
    genetypes.push_back("IG_V_gene");
    genetypes.push_back("IG_LV_pseudogene");
    genetypes.push_back("IG_D_gene");
    genetypes.push_back("IG_J_gene");
    genetypes.push_back("IG_J_pseudogene");
    genetypes.push_back("IG_C_gene");
    genetypes.push_back("IG_C_pseudogene");
    genetypes.push_back("TR_V_gene");
    genetypes.push_back("TR_V_pseudogene");
    genetypes.push_back("TR_D_gene");
    genetypes.push_back("TR_J_gene");
    genetypes.push_back("TR_J_pseudogene");
    genetypes.push_back("TR_C_gene");    
  }

  // read GTF file
  notice("Opening GTF file %s...", gtfFile.c_str());
  gtf inGTF(gtfFile.c_str(), &genetypes, addChrPrefix, removeChrPrefix);
  notice("Finished reading GTF file %s...", gtfFile.c_str());

  // read the UMI information per each barcode
  sprintf(fname, "%s.umi.gz", plpPrefix.c_str());
  notice("Reading UMI information from %s.umi.gz..", plpPrefix.c_str());
  tsv_reader tsv_umif(fname);
  std::string chrom;
  int32_t beg1, end0;

  // count every possible GTF elements
  std::map<gtfElement*, std::map<int32_t,int32_t> > dgeMap;
  std::map<std::string, int64_t> typeCount;
  
  for (int64_t line = 1; tsv_umif.read_line() > 0; ++line) {
    int32_t old_id = tsv_umif.int_field_at(0);
    if ( id_cel2dge.find(old_id) == id_cel2dge.end() ) {
      if ( skipbcd > 0 ) continue; // if anything was skipped, missing a specific ID is fine.
      else error("Cannot find barcode ID %d", old_id);
    }
    int32_t new_id = id_cel2dge[old_id];

    if ( line % 1000000 == 0 ) {
      notice("Processing %d UMIs over %d barcodes", line, new_id);
    }
    
    // parse the current UMI to add to the current gene count profile
    // for processing the current UMI, genes, transcripts, and exons are only counted once each time
    std::set<gtfElement*> sElems; // store every element in GTF field for the UMI
    for(int32_t j=3; j < tsv_umif.nfields; ++j) {
      posLocus::parseRegion(tsv_umif.str_field_at(j), chrom, beg1, end0);
      inGTF.findOverlappingElements(chrom.c_str(), beg1, end0, sElems);
    }

    // focus only on exons
    std::set<gtfElement*> umiElems;
    for(std::set<gtfElement*>::iterator it = sElems.begin(); it != sElems.end(); ++it) {
      if ( (*it)->type == "exon" ) {
	umiElems.insert(*it);                   // insert the exon
	umiElems.insert((*it)->parent);         // insert the transcript
	umiElems.insert((*it)->parent->parent); // insert the gene
      }
    }

    for(std::set<gtfElement*>::iterator it = umiElems.begin(); it != umiElems.end(); ++it) {
      ++(dgeMap[*it][new_id]);
      ++typeCount[(*it)->type];
    }
  }

  // write outputs for every possible elements
  // build a dictionary of gene, transcript, and exons
  std::map<std::string, std::vector<gtfElement*> >       vecElems;
  std::map<std::string, std::map<gtfElement*, int32_t> > mapElems;
  for(gtf::gtf_chr_it_t it = inGTF.mmap.begin(); it != inGTF.mmap.end(); ++it) { // iterate over chromosomes
    for(gtf::gtf_elem_it_t jt = it->second.begin(); jt != it->second.end(); ++jt) {
      vecElems[jt->second->type].push_back(jt->second);
      mapElems[jt->second->type][jt->second] = (int32_t)vecElems[jt->second->type].size()-1;
    }
  }

  std::map<std::string, htsFile*> mtxFiles;  
  std::map<std::string,int64_t>::iterator it;
  for(it = typeCount.begin(); it != typeCount.end(); ++it) {
    std::string subdir = (outPrefix + "/" + it->first);
    ret = mkdir(subdir.c_str(), 0777);
    if ( ret )
      error("Cannot create directory %s", subdir.c_str());

    // write barcodes file
    htsFile* wbcd = hts_open( (subdir+"/barcodes.tsv").c_str(), "w");
    if ( wbcd == NULL )
      error("Cannot open file %s/barcodes.tsv for writing", subdir.c_str());
    for(int32_t i=0; i < (int32_t)dgebcds.size(); ++i) 
      hprintf(wbcd,"%s\n",dgebcds[i].c_str());
    hts_close(wbcd);

    // write genes.tsv file
    htsFile* wgene = hts_open( (subdir+"/genes.tsv").c_str(), "w");
    if ( wgene == NULL )
      error("Cannot open file %s/genes.tsv for writing", subdir.c_str());
    std::vector<gtfElement*>& vElems = vecElems[it->first];
    if ( it->first == "gene" ) {
      for(int32_t i=0; i < (int32_t)vElems.size(); ++i) {
	gtfGene* g = (gtfGene*)vElems[i];
	hprintf(wgene,"%s\t%s\n", g->geneId.c_str(), g->geneName.c_str());
      }
    }
    else if ( it->first == "transcript" ) {
      for(int32_t i=0; i < (int32_t)vElems.size(); ++i) {
	gtfTranscript* t = (gtfTranscript*)vElems[i];
	gtfGene*       g = (gtfGene*)t->parent;
	hprintf(wgene,"%s\t%s:%s\n", t->transcriptId.c_str(), g->geneName.c_str(), t->transcriptId.c_str());
      }      
    }
    else {
      for(int32_t i=0; i < (int32_t)vElems.size(); ++i) {
	gtfElement* e = vElems[i];
	gtfTranscript* t = (gtfTranscript*)e->parent;
	gtfGene*       g = (gtfGene*)t->parent;
	hprintf(wgene,"%s:%s:%d-%d\t%s:%s:%s:%d-%d\n", t->transcriptId.c_str(), g->seqname.c_str(), e->locus.beg1, e->locus.end0, g->geneName.c_str(), t->transcriptId.c_str(), g->seqname.c_str(), e->locus.beg1, e->locus.end0);
      }
    }
    hts_close(wgene);

    // open matrix.mtx file and write headers
    htsFile* wmtx = hts_open( (subdir+"/matrix.mtx").c_str(), "w");
    if ( wmtx == NULL )
      error("Cannot open file %s/matrix.mtx for writing", subdir.c_str());
    hprintf(wmtx, "%%%%MatrixMarket matrix coordinate integer general\n%%\n");
    hprintf(wmtx, "%u %u %lld\n", vElems.size(), dgebcds.size(), it->second);
    mtxFiles[it->first] = wmtx;
  }

  notice("Writing sparse matices of UMI counts");
  int32_t nelems = 0;
  for(std::map<gtfElement*, std::map<int32_t,int32_t> >::iterator it = dgeMap.begin();
      it != dgeMap.end(); ++it) {
    if ( nelems % 10000 == 0 )
      notice("Processing %d / %u GTF elements...", nelems, dgeMap.size());
    gtfElement* e = it->first;
    gtfElement* root = e;
    while( root->parent != NULL ) root = root->parent;
    std::map<gtfElement*, int32_t>& igenes = mapElems[e->type];
    if ( igenes.find(e) == igenes.end() ) {
      error("Cannot find element %p with type %s at gene %s at chromosome %s", (void*)e, e->type.c_str(), ((gtfGene*)root)->geneId.c_str(), ((gtfGene*)root)->seqname.c_str());
    }
    int32_t igene = mapElems[e->type][e] + 1;
    htsFile* hf = mtxFiles[e->type];
    for(std::map<int32_t,int32_t>::iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {
												    hprintf(hf, "%d %d %d\n", igene, jt->first + 1, jt->second);
    }
    ++nelems;
  }

  notice("Finished processing %u / %u GTF elements...", dgeMap.size(), dgeMap.size());

  for(it = typeCount.begin(); it != typeCount.end(); ++it) {
    hts_close(mtxFiles[it->first]);
  }

  notice("Analysis Finished");  
  
  return 0;
}
