#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_variant_key.h"
#include "tsv_reader.h"

int32_t cmdVcfIBSMatrix(int32_t argc, char** argv) {
  std::string listVcf; // list of small VCF files
  std::string inVcf;   
  std::string refVcf;  // large VCF file containing many samples
  std::string out;     // output file name containing IBS matrix
  std::string reg;
  int32_t verbose = 20000;

  bcf_vfilter_arg vfilt;
  
  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("list",&listVcf, "List of VCF files, typically single sample")
    LONG_STRING_PARAM("vcf",&inVcf, "Input VCF")    
    LONG_STRING_PARAM("panel",&refVcf,"VCF panel containing genotypes to compare")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)    
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file to store IBS matrix")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")    
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( ( listVcf.empty() + inVcf.empty() != 1 ) || refVcf.empty() || out.empty() ) {
    error("[E:%s:%d %s] --list or --vcf, --panel --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  notice("Analysis started");
  
  // for now, let's load all genotypes together
  std::map<variantKeyS, std::map<int32_t,int32_t> > mGeno;

  std::vector<GenomeInterval> intervals;
  std::vector<std::string> listIDs;
  std::vector<int32_t> nsamplesPerVcf;
  int32_t isample = 0;
  int32_t* gts = NULL;
  int32_t n_gts = 0, geno = 0;

  // read the list of VCFs
  std::vector<std::string> vcfs;
  if ( !listVcf.empty() ) {
    tsv_reader rdr(listVcf.c_str());
    for(int32_t i=0; rdr.read_line() > 0; ++i) {
      vcfs.push_back(rdr.str_field_at(0));
    }
    rdr.close();
  }
  else {
    vcfs.push_back(inVcf);
  }
  
  for(int32_t i=0; i < vcfs.size(); ++i) {
    BCFOrderedReader odr(vcfs[i].c_str(), intervals);
    notice("Processing BCF file #%d - %s",i, vcfs[i].c_str());
    bcf1_t* iv = bcf_init();
    int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
    nsamplesPerVcf.push_back(nsamples);
    for(int32_t j=0; j < nsamples; ++j) {
      listIDs.push_back(bcf_hdr_get_sample_name(odr.hdr,j));
    }

    while(odr.read(iv)) {
      variantKeyS key(odr.hdr, iv);
      std::map<int32_t,int32_t>& var = mGeno[key]; // will create new if needed
      bcf_get_genotypes(odr.hdr, iv, &gts, &n_gts);
      for(int32_t j=0; j < nsamples; ++j) {
	geno = bcf_gt_allele(gts[2*j]) + bcf_gt_allele(gts[2*j+1]);
	var[isample+j] = geno;
      }
      free(gts);
      gts = NULL;
      n_gts = 0;
    }
    odr.close();

    isample += nsamples;
  }

  notice("Finished loading input VCF files.. Loading the reference VCF..");

  // Process the reference BCF
  // handle filter string
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }  
  BCFOrderedReader odr(refVcf.c_str(), intervals);
  
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
    filter_init(odr.hdr, filter_str.c_str());

  // handle --apply-filtrs
  std::vector<int32_t> req_flt_ids;
  if ( !vfilt.required_filters.empty() ) {
    for(int32_t i=0; i < (int32_t)vfilt.required_filters.size(); ++i) {
      req_flt_ids.push_back(bcf_hdr_id2int(odr.hdr, BCF_DT_ID, vfilt.required_filters[i].c_str()));
    }
  }  

  int32_t nVariant = 0;
  int32_t nRefSamples = bcf_hdr_nsamples(odr.hdr);
  std::vector<std::string> refIDs;
  for(int32_t j=0; j < nRefSamples; ++j) {
    refIDs.push_back(bcf_hdr_get_sample_name(odr.hdr,j));
  }  

  notice("Started Reading site information from VCF file, identifying %d samples", nRefSamples);

  if ( nRefSamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

  int32_t* p_gt = NULL;
  int32_t n_gt = 0;
  int32_t nskip = 0, nparsed = 0;

  bcf1_t* iv = bcf_init();

  // tensors of IBS
  int32_t nListSamples = (int32_t)listIDs.size();
  std::vector< std::vector< std::vector<int32_t> > > list2ref2ibs(nListSamples);
  for(int32_t i=0; i < nListSamples; ++i) {
    list2ref2ibs[i].resize(nRefSamples);
    for(int32_t j=0; j < nRefSamples; ++j) {
      list2ref2ibs[i][j].resize(9);
    }
  }

  std::vector<int32_t> genos(nRefSamples);
  for(int32_t k=0; odr.read(iv); ++k) {  // read marker
    // periodic message to user
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d. Skipped %d filtered markers, retaining %d variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nskip, nVariant);

    variantKeyS key(odr.hdr, iv);
    // check if the variant need to be parsed
    std::map<variantKeyS, std::map<int32_t,int32_t> >::iterator itG = mGeno.find(key);
    if ( itG == mGeno.end() ) {
      ++nskip;
      continue; // no need to parse the genotypes
    }
    
    // unpack FILTER column
    bcf_unpack(iv, BCF_UN_FLT);

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

    ++nVariant;

    // extract genotype and apply genotype level filter
    if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }
    for(int32_t i=0; i < nRefSamples; ++i) {
      int32_t g1 = p_gt[2*i];
      int32_t g2 = p_gt[2*i+1];
      if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
	genos[i] = -1;
	//geno = 0;
      }
      else {
	genos[i] = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
      }
    }

    // iterate sparse genotypes
    for(std::map<int32_t,int32_t>::iterator it = itG->second.begin(); it != itG->second.end(); ++it) {
      if ( it->second >= 0 ) {
	for(int32_t i=0; i < nRefSamples; ++i) {
	  if ( genos[i] >= 0 )
	    ++list2ref2ibs[it->first][i][it->second*3 + genos[i]];
	}
      }
    }
    free(p_gt);
    p_gt = NULL;
    n_gt = 0;
  }

  notice("Finished Processing %d markers across %d samples, Skipping %d filtered markers", nVariant, nskip);

  htsFile* wf = hts_open(out.c_str(), (out.substr(out.size()-3).compare(".gz") == 0) ? "wz" : "w");
  hprintf(wf, "LIST_ID\tREF_ID\tREF_REF\tREF_HET\tREF_HOM\tHET_REF\tHET_HOM\tHOM_REF\tHOM_HET\tHOM_HOM\tCONC_ALL\tCONC_NREF\n");
  for(int32_t i=0; i < nListSamples; ++i) {
    for(int32_t j=0; j < nRefSamples; ++j) {
      std::vector<int32_t>& v = list2ref2ibs[i][j];
      int32_t nconc = v[0] + v[4] + v[8];
      int32_t ndisc = v[1] + v[2] + v[3] + v[5] + v[6] + v[7];
      hprintf(wf, "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.5lf\t%.5lf\n", listIDs[i].c_str(), refIDs[j].c_str(), v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], nconc/(nconc+ndisc+1e-6), (nconc-v[0])/(nconc+ndisc-v[0]+1e-6));
    }
  }
  hts_close(wf);

  notice("Analysis Finished");

  return 0;
}

