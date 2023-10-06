#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "sex_ploidy_map.h"

int32_t cmdVcfSqueeze(int32_t argc, char** argv) {
  std::vector<std::string> inVcfs;
  std::string out;
  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;
  int32_t verbose = 1000;
  int32_t xStart = 2699520;
  int32_t xStop = 154931044;
  int32_t minDPmaleX = -1;
  std::string xLabel("X");
  std::string yLabel("Y");
  std::string mtLabel("MT");
  std::string sexMap;  

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input VCFs", NULL)
    LONG_MULTI_STRING_PARAM("in",&inVcfs, "Input VCF/BCF files")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)    
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")    

    LONG_PARAM_GROUP("Genotype Filtering Options", NULL)    
    LONG_INT_PARAM("minDP",&gfilt.minDP,"Minimum depth threshold for filtering genotypes")
    LONG_INT_PARAM("minGQ",&gfilt.minGQ,"Minimum depth threshold for filtering genotypes")
    LONG_INT_PARAM("minDP-male-X",&minDPmaleX,"Minimum depth threshold male X chromosomes for filtering genotypes")        

    LONG_PARAM_GROUP("Sex Chromosomes (for --minDP option)",NULL)
    LONG_STRING_PARAM("sex-map",&sexMap, "Sex map file, containing ID and sex (1 for male and 2 for female) for each individual")    
    LONG_STRING_PARAM("x-label", &xLabel, "Contig name for X chromosome")
    LONG_STRING_PARAM("y-label", &yLabel, "Contig name for Y chromosome")
    LONG_STRING_PARAM("mt-label", &mtLabel, "Contig name for MT chromosome")
    LONG_INT_PARAM("x-start", &xStart, "Start base position of non-PAR region in X chromosome")
    LONG_INT_PARAM("x-stop",  &xStop,  "End base position of non-PAR region in X chromosome")        

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")    
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // sanity check of input arguments
  if ( inVcfs.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // if --minDP-male-X was not set, make it identical to minDP
  if ( minDPmaleX < 0 ) minDPmaleX = gfilt.minDP;

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

  BCFOrderedWriter odw(out.c_str(),0);
  odw.set_hdr(odr->hdr);
  odw.write_hdr();  

  int32_t nsamples = bcf_hdr_nsamples(odr->hdr);
  int32_t n_gts = 0, n_flds =0, max_n_flds = 0;
  int32_t* gts = NULL; //(int32_t*) calloc( nsamples * 2, sizeof(int32_t) );
  int32_t* flds = NULL; //(int32_t*) calloc( nsamples, sizeof(int32_t) );  

  // read a specific marker position
  int32_t k = 0, nskip = 0;

  sex_ploidy_map spmap(xLabel, yLabel, mtLabel, xStart, xStop);
  spmap.load_sex_map_file(sexMap.empty() ? NULL : sexMap.c_str(), odr->hdr);      
 
  for(int32_t l=0; l < (int32_t)inVcfs.size(); ++l) {
    notice("PROCESSING %d-th VCF file %s", l+1, inVcfs[l].c_str());
    
    for(; odr->read(iv); ++k) {  // read marker
      if ( k % verbose == 0 )
	notice("Processing %d markers at %s:%d. Skipped %d markers", k, bcf_hdr_id2name(odr->hdr, iv->rid), iv->pos+1, nskip);
      int8_t* ploidies = spmap.get_ploidies(iv);
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

      if ( gfilt.minDP > 0 ) {
	if ( (n_flds = bcf_get_format_int32(odr->hdr, iv, "AD", &flds, &max_n_flds)) < 0 ) {
	  error("[E:%s:%d %s] Cannot find the field AD from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr->hdr, iv->rid), iv->pos+1);
	  
	}

        int32_t nallele = n_flds / nsamples;
        for(int32_t i=0; i < nsamples; ++i) {
          flds[i] = flds[nallele*i];
          for (int32_t j=1; j < nallele; ++j) {
            flds[i] += flds[nallele*i+j];
          } 
        }
        n_flds = nsamples;

	for(int32_t i=0; i < nsamples; ++i) {
	  if ( ( ploidies[i] == 0 ) ||
	       ( ( ploidies[i] == 1 ) && ( flds[i] < minDPmaleX ) ) ||
	       ( ( ploidies[i] == 2 ) && ( flds[i] < gfilt.minDP ) ) ) {
	    gts[2*i] = bcf_gt_missing;
	    gts[2*i+1] = bcf_gt_missing;	    
	  }
	}
      }
      
      if ( gfilt.minGQ > 0 ) {
	if ( bcf_get_format_int32(odr->hdr, iv, "GQ", &flds, &max_n_flds) < 0 ) {
	  error("[E:%s:%d %s] Cannot find the field GQ from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr->hdr, iv->rid), iv->pos+1);
	  
	}
	for(int32_t i=0; i < nsamples; ++i) {
	  if ( flds[i] < gfilt.minGQ ) {
	    gts[2*i] = bcf_gt_missing;
	    gts[2*i+1] = bcf_gt_missing;
	  }
	}
      }

      bcf1_t* nv = bcf_init();

      nv->n_sample = iv->n_sample;
      nv->rid = iv->rid;
      nv->pos = iv->pos;
      nv->rlen = iv->rlen;
      nv->qual = iv->qual;
      nv->n_allele = iv->n_allele;

      bcf_update_alleles(odw.hdr, nv, (const char**)iv->d.allele, iv->n_allele);
      bcf_update_filter(odw.hdr, nv, iv->d.flt, iv->d.n_flt);
      
      bcf_unpack(nv, BCF_UN_ALL);

      // transfer INFO fields
      for(int32_t i=0; i < iv->n_info; ++i) {
	bcf_info_t& info = iv->d.info[i];
	if ( info.type != BCF_BT_NULL ) {
	  const char* tag = bcf_hdr_int2id(odr->hdr,BCF_DT_ID,info.key);
	  int32_t htype = bcf_hdr_id2type(odr->hdr,BCF_HL_INFO,info.key);
	  int32_t ntmp_arr = 0;
	  void* tmp_arr = NULL;
	  int32_t ret = bcf_get_info_values(odr->hdr, iv, tag, &tmp_arr, &ntmp_arr, htype);
	  if ( ret > 0 ) {
	    if ( bcf_update_info(odw.hdr, nv, tag, tmp_arr, ntmp_arr, htype) < 0 ) {
	      fprintf(stderr,"Cannot write INFO field %s\n",tag);
	      abort();
	    }
	  }
	  else {
	    fprintf(stderr,"Cannot retrieve INFO field %s\n",tag);
	    abort();
	  }
	  free(tmp_arr);
	}
      }

      bcf_update_format_int32(odw.hdr, nv, "GT", gts, nsamples * 2);
      
      odw.write(nv);
      bcf_destroy(nv);
    }
    
    delete odr;
    if ( l+1 < (int32_t)inVcfs.size() ) 
      odr = new BCFOrderedReader(inVcfs[l+1], intervals);
    else
      odr = NULL;
  }
  odw.close();

  notice("Analysis finished");

  return 0;
}

