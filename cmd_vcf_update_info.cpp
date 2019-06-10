#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"

int32_t cmdVcfUpdateInfo(int32_t argc, char** argv) {
  std::vector<std::string> inVcfs;
  std::string inFasta;
  std::string out;
  bcf_vfilter_arg vfilt;
  int32_t verbose = 10000;
  bool siteOnly = false;

  // update INFO field using sample level information
  bool updateAC   = false;
  bool updateAF   = false;
  bool updateGC   = false;
  bool updateFMIS = false;
  bool updateMAC  = false;
  bool updateMAF  = false;
  
  // update INFO field requiring the reference FASTA file
  bool updateCpG   = false;
  bool updateKmer  = false;
  bool updateExAWS = false;
  bool updateFzAWS = false;
  std::vector<int32_t> kmerLengths;  

  bool updateALL   = false;  

  // Wish list
  // Variation of depth
  // Hemizygosity?
  // HWE test statistics and ISAF?

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input VCFs", NULL)
    LONG_MULTI_STRING_PARAM("in",&inVcfs, "Input VCF/BCF files")
    LONG_STRING_PARAM("ref",&inFasta, "Input FASTA file")    

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)    
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")    

    LONG_PARAM_GROUP("Update INFO fields that require sample-level genotype information", NULL)
    LONG_PARAM("update-AC",&updateAC, "Update/annotate allele count (AC) information, (together with AN)")
    LONG_PARAM("update-AF",&updateAF, "Update/annotate allele frequency information (based on AC). Must be used with --update-AC")
    LONG_PARAM("update-GC",&updateGC, "Update/annotate genotype frequency (based on AC and AN). Must be used with --update-AC")
    LONG_PARAM("update-FMIS",&updateFMIS, "Update/annotate fraction of missing genotypes")
    LONG_PARAM("update-MAC",&updateMAC, "Update/annotate minor allele count. Must be used with --update-AC")
    LONG_PARAM("update-MAF",&updateMAF, "Update/annotate minor allele frequency (based on MAC and AN). Must be used with --update-AC and --update-AF")

    LONG_PARAM_GROUP("Update INFO fields that require reference FASTA file", NULL)
    LONG_PARAM("update-CpG",  &updateCpG,    "Update/annotate CpG information")
    LONG_PARAM("update-Kmer", &updateKmer,   "Update/annotate k-mer context information (must be used with --kmer-length)")
    LONG_PARAM("update-ExAWS",&updateExAWS,  "Update/annotate exact allelic window shift")
    LONG_PARAM("update-FzAWS",&updateFzAWS,  "Update/annotate fuzzy allelic window shift")

    LONG_PARAM_GROUP("Update all available INFO fields", NULL)
    LONG_PARAM("update-all",  &updateALL,    "Update/annotate all available INFO fields")    

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_PARAM("site-only", &siteOnly,   "Do not write genotype information, and writes only site information (up to INFO field) in output VCF/BCF")    
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")    
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // sanity check of input arguments
  if ( inVcfs.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  if ( !( updateAC || updateAF || updateGC || updateFMIS || updateMAC || updateMAF || updateCpG || updateKmer || updateExAWS || updateFzAWS || updateALL ) )
    error("[E:%s:%d %s] At least one --update-*** argument is required",__FILE__,__LINE__,__FUNCTION__);

  if ( inFasta.empty() && ( updateCpG || updateKmer || updateExAWS || updateFzAWS ) )
    error("[E:%s:%d %s] At least one option requires --ref parameter",__FILE__,__LINE__,__FUNCTION__);

  if ( ( !updateAC ) && ( updateAF || updateGC || updateMAC || updateMAF ) )
    error("[E:%s:%d %s] --update-AC is required to enable the specified options ",__FILE__,__LINE__,__FUNCTION__) ;

  if ( ( !updateAF ) && ( updateMAF ) )
    error("[E:%s:%d %s] --update-AF is required to enable the specified options ",__FILE__,__LINE__,__FUNCTION__) ;

  for(int32_t i=0; i < (int32_t)kmerLengths.size(); ++i) {
    if ( kmerLengths[i] % 2 == 0 )
      error("[E:%s:%d %s] --kmer-length %d is invalud. Must be an odd number",__FILE__,__LINE__,__FUNCTION__, kmerLengths[i]) ;      
  }


  // read input BCF/VCF file(s)
  std::vector<GenomeInterval> intervals;
  BCFOrderedReader* odr = new BCFOrderedReader(inVcfs[0], intervals);
  bcf1_t* iv = bcf_init();

  // handle filter string
  std::string filter_str;
  int32_t filter_logic = 0;
  if ( vfilt.include_expr.empty() ) {
    if ( vfilt.exclude_expr.empty() ) {
      // do nothing
    }
    else {
      filter_str = vfilt.exclude_expr;
      filter_logic |= FLT_EXCLUDE;
    }
  }
  else {
    if ( vfilt.exclude_expr.empty() ) {
      filter_str = vfilt.include_expr;
      filter_logic |= FLT_INCLUDE;      
    }
    else {
      error("[E:%s:%d %s] Cannot use both --include-expr and --exclude-expr options",__FILE__,__LINE__,__FUNCTION__);
    }    
  }

  filter_t* filt = NULL;
  if ( filter_logic != 0 )
    filter_init(odr->hdr, filter_str.c_str());

  // handle --apply-filtrs
  std::vector<int32_t> req_flt_ids;
  if ( !vfilt.required_filters.empty() ) {
    for(int32_t i=0; i < (int32_t)vfilt.required_filters.size(); ++i) {
      req_flt_ids.push_back(bcf_hdr_id2int(odr->hdr, BCF_DT_ID, vfilt.required_filters[i].c_str()));
    }
  }

  int32_t nsamples = bcf_hdr_nsamples(odr->hdr);
  int32_t n_gts = 0; //, n_flds = 0;
  int32_t* gts = (int32_t*) calloc( nsamples * 2, sizeof(int32_t) );
  //int32_t* flds = (int32_t*) calloc( nsamples, sizeof(int32_t) );

  if ( updateALL ) {
    if ( nsamples > 0 ) {  // if genotypes are available
      updateAC = true;
      updateAF = true;
      updateGC = true;
      updateFMIS = true;
      updateMAC = true;
      updateMAF = true;
    }
    
    if ( !inFasta.empty() ) {
      updateCpG   = true;
      updateKmer  = true;
      updateExAWS = true;
      updateFzAWS = true;
      if ( kmerLengths.empty() ) {
	kmerLengths.push_back(3);
	kmerLengths.push_back(5);	
      }
    }
  }

  notice("The following options are enabled:");
  if ( updateAC ) notice("\t--update-AC is enabled");
  if ( updateAF ) notice("\t--update-AF is enabled");
  if ( updateGC ) notice("\t--update-GC is enabled");
  if ( updateFMIS ) notice("\t--update-FMIS is enabled");
  if ( updateMAC ) notice("\t--update-MAC is enabled");
  if ( updateMAF ) notice("\t--update-MAF is enabled");
  if ( updateCpG ) notice("\t--update-CpG is enabled");
  if ( updateKmer ) notice("\t--update-Kmer is enabled");
  if ( updateExAWS ) notice("\t--update-ExAWS is enabled");
  if ( updateFzAWS ) notice("\t--update-FzAWS is enabled");  

  BCFOrderedWriter odw(out.c_str(),0);

  if ( siteOnly ) {
    bcf_hdr_t* hnull = bcf_hdr_subset(odr->hdr, 0, 0, 0);
    bcf_hdr_remove(hnull, BCF_HL_FMT, NULL);
    odw.set_hdr(hnull);
  }
  else {
    odw.set_hdr(odr->hdr);
  }

  char buffer[65536];  

  // check the existence of header and create one if needed
  if ( updateAC && ( bcf_hdr_id2int(odr->hdr, BCF_DT_ID, "AC" ) < 0 ) ) {
    sprintf(buffer,"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count of each non-reference allele\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( updateAC && ( bcf_hdr_id2int(odr->hdr, BCF_DT_ID, "AN" ) < 0 ) ) {
    sprintf(buffer,"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total allele count\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }  
  if ( updateAF && ( bcf_hdr_id2int(odr->hdr, BCF_DT_ID, "AF" ) < 0 ) ) {
    sprintf(buffer,"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency of each non-reference allele\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( updateGC && ( bcf_hdr_id2int(odr->hdr, BCF_DT_ID, "GC" ) < 0 ) ) {
    sprintf(buffer,"##INFO=<ID=GC,Number=G,Type=Integer,Description=\"Diploid genotype count\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( updateGC && ( bcf_hdr_id2int(odr->hdr, BCF_DT_ID, "GN" ) < 0 ) ) {
    sprintf(buffer,"##INFO=<ID=GC,Number=1,Type=Integer,Description=\"Total diploid genotype count\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }  
  if ( updateFMIS && ( bcf_hdr_id2int(odr->hdr, BCF_DT_ID, "FMIS" ) < 0 ) ) {
    sprintf(buffer,"##INFO=<ID=FMIS,Number=1,Type=Float,Description=\"Fraction of missing genotype\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( updateMAC && ( bcf_hdr_id2int(odr->hdr, BCF_DT_ID, "MAC" ) < 0 ) ) {
    sprintf(buffer,"##INFO=<ID=MAC,Number=1,Type=Integer,Description=\"Minor allele count of non-reference allele. Available only for biallelic variants\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( updateMAF && ( bcf_hdr_id2int(odr->hdr, BCF_DT_ID, "MAF" ) < 0 ) ) {
    sprintf(buffer,"##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency of non-reference allele. Available only for biallelic variants\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }   
  if ( updateCpG && ( bcf_hdr_id2int(odr->hdr, BCF_DT_ID, "CPG" ) < 0 ) ) {
    sprintf(buffer,"##INFO=<ID=CPG,Number=0,Type=Flag,Description=\"Indicate that the position is a CpG site\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( updateKmer ) {
    for(int32_t i=0; i < (int32_t)kmerLengths.size(); ++i) {
      std::string key = "KMER";
      catprintf(key,"%d",kmerLengths[i]);
      if ( bcf_hdr_id2int(odr->hdr, BCF_DT_ID, key.c_str() ) < 0 ) {
	sprintf(buffer,"##INFO=<ID=KMER%d,Number=1,Type=String,Description=\"Kmer context in %d centered at the variant position\">\n", kmerLengths[i], kmerLengths[i]);
	bcf_hdr_append(odw.hdr, buffer);
      }
    }
  }
  if ( updateCpG && ( bcf_hdr_id2int(odr->hdr, BCF_DT_ID, "ExAWS" ) < 0 ) ) {
    sprintf(buffer,"##INFO=<ID=ExAWS,Number=1,Type=Integer,Description=\"Exact allelic window shift length\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( updateCpG && ( bcf_hdr_id2int(odr->hdr, BCF_DT_ID, "FzAWS" ) < 0 ) ) {
    sprintf(buffer,"##INFO=<ID=FzAWS,Number=1,Type=Integer,Description=\"Fuzzy allelic window shift length\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }    
  odw.write_hdr();    
  

  // read a specific marker position
  int32_t k = 0, nskip = 0;
  int32_t acs[256];
  int32_t gcs[65536];
 
  for(int32_t l=0; l < (int32_t)inVcfs.size(); ++l) {
    notice("PROCESSING %d-th VCF file %s", l+1, inVcfs[l].c_str());
    
    for(; odr->read(iv); ++k) {  // read marker
      if ( k % verbose == 0 )
	notice("Processing %d markers at %s:%d. Skipped %d markers", k, bcf_hdr_id2name(odr->hdr, iv->rid), iv->pos+1, nskip);
      
      bcf_unpack(iv, BCF_UN_ALL);
      
      // check --apply-filters
      bool has_filter = req_flt_ids.empty() ? true : false;
      if ( ! has_filter ) {
	//notice("%d %d", iv->d.n_flt, (int32_t)req_flt_ids.size());
	for(int32_t i=0; i < iv->d.n_flt; ++i) {
	  for(int32_t j=0; j < (int32_t)req_flt_ids.size(); ++j) {
	    if ( req_flt_ids[j] == iv->d.flt[i] )
	      has_filter = true;
	  }
	}
      }
      
      if ( ! has_filter ) { ++nskip; continue; }
      
      // check filter logic
      if ( filt != NULL ) {
	int32_t ret = filter_test(filt, iv, NULL);
	if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
	else if ( ret ) { has_filter = false; }
      }
      
      if ( ! has_filter ) { ++nskip; continue; }

      // extract genotype and apply genotype level filter
      if ( bcf_get_genotypes(odr->hdr, iv, &gts, &n_gts) < 0 ) {
	error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr->hdr, iv->rid), iv->pos+1);
      }

      int32_t an = 0, nmis = 0, ngeno = 0;
      memset(acs, 0, sizeof(int32_t)*256);
      memset(gcs, 0, sizeof(int32_t)*65535);  

      bcf1_t* nv = bcf_dup(iv);
      bcf_unpack(nv, BCF_UN_ALL);

      // calculate AN, AC, NMIS, AC, GC
      for(int32_t i=0; i < nsamples; ++i) {
	int32_t a1 = bcf_gt_is_missing(gts[2*i]) ? -1 : bcf_gt_allele(gts[2*i]);
	int32_t a2 = bcf_gt_is_missing(gts[2*i+1]) ? -1 : bcf_gt_allele(gts[2*i+1]);

	if ( ( a1 < 0 ) && ( a2 < 0 ) ) {
	  ++nmis;
	  ++ngeno;
	}
	else if ( a1 < 0 ) {
	  ++an;
	  ++acs[a2];
	}
	else if ( a2 < 0 ) {
	  ++an;
	  ++acs[a1];
	}
	else {
	  an += 2;
	  ++acs[a1];
	  ++acs[a2];
	  ++ngeno;
	  if ( a1 < a2 ) { ++gcs[a2*(a2+1)/2+a1]; }
	  else { ++gcs[a1*(a1+1)/2+a2]; }   
	}
      }

      if ( updateAC ) bcf_update_info_int32(odw.hdr, nv, "AN", &an, 1); // update AN;
      if ( updateAC ) bcf_update_info_int32(odw.hdr, nv, "AC", &acs[1], nv->n_allele-1); // update AC
      if ( updateAF ) {
	float afs[256];
	if ( an == 0 ) { memset(afs,0,nv->n_allele*sizeof(int32_t)); }
	else {
	  for(int32_t i=0; i < nv->n_allele; ++i) {
	    afs[i] = (float)acs[i]/(float)an;
	  }
	  bcf_update_info_float(odw.hdr, nv, "AF", &afs[1], nv->n_allele-1); // update AF
	}
	if ( updateMAF && ( nv->n_allele == 2 ) ) {
	  float maf = afs[1] > 0.5 ? afs[0] : afs[1];
	  bcf_update_info_float(odw.hdr, nv, "MAF", &maf, 1); // update MAF
	}
      }
      if ( updateMAC && ( nv->n_allele == 2 ) ) {
	int32_t mac = acs[0] > acs[1] ? acs[1] : acs[0];
	bcf_update_info_int32(odw.hdr, nv, "MAC", &mac, 1); // update MAC
      }
      if ( updateGC ) {
	int32_t nG = nv->n_allele * (nv->n_allele+1) / 2;
	bcf_update_info_int32(odw.hdr, nv, "GC", &gcs, nG); // update GC
	int32_t gn = ngeno - nmis;
	bcf_update_info_int32(odw.hdr, nv, "GN", &gn, 1); // update GN
      }
      if ( updateFMIS ) {
	float fmis = (ngeno == 0) ? 0 : ( (float)nmis/(float)ngeno );
	bcf_update_info_float(odw.hdr, nv, "FMIS", &fmis, 1); // update FMIS	
      }

      if ( siteOnly ) {
	bcf_subset(odw.hdr, nv, 0, 0);
      }  

      odw.write(nv);
      bcf_destroy(nv);
    }
    
    delete odr;  // delete BCFReader if finished reading
    if ( l+1 < (int32_t)inVcfs.size() ) 
      odr = new BCFOrderedReader(inVcfs[l+1], intervals);  // read next file if available
    else
      odr = NULL;
  }
  odw.close();

  notice("Analysis finished");

  return 0;
}

