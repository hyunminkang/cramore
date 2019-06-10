#include "bcf_filter_arg.h"
#include "cramore.h"
#include "bcf_ordered_reader.h"
#include "compact_matrix.h"

// goal -- for a given chunk (e.g. 100kb), identify pairs with no IBS0 pairs
int32_t cmdVcfIBS0Summary(int32_t argc, char** argv) {
  std::string inVcf;
  std::string out;
  int32_t min_hom_gts = 1;
  int32_t verbose = 1000;
  int32_t batch_size = 10000;
  std::string reg;
  int32_t min_variant = 1;

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;
  
  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)    
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("min-hom",&min_hom_gts, "Minimum number of homozygous genotypes to be counted for IBS0")
    LONG_INT_PARAM("batch-size",&batch_size, "Size of batches (in # of samples) to calculate the no-IBS0 pairs")
    LONG_INT_PARAM("min-variant",&min_variant, "Minimum number of variants to present to have output file")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")    
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }
  BCFOrderedReader odr(inVcf, intervals);
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
    filter_init(odr.hdr, filter_str.c_str());

  // handle --apply-filtrs
  std::vector<int32_t> req_flt_ids;
  if ( !vfilt.required_filters.empty() ) {
    for(int32_t i=0; i < (int32_t)vfilt.required_filters.size(); ++i) {
      req_flt_ids.push_back(bcf_hdr_id2int(odr.hdr, BCF_DT_ID, vfilt.required_filters[i].c_str()));
    }
  }  

  int32_t nVariant = 0;
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);

  std::vector<int32_t> nRRs, nAAs;

  notice("Started Reading site information from VCF file, identifying %d samples", nsamples);

  if ( nsamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

  bitmatrix bmatRR(nsamples);
  bitmatrix bmatAA(nsamples);

  int32_t* p_gt = NULL;
  int32_t n_gt = 0;
  int32_t nskip = 0, nmono = 0;
  uint8_t* gtRR = (uint8_t*)calloc(nsamples, sizeof(uint8_t));
  uint8_t* gtAA = (uint8_t*)calloc(nsamples, sizeof(uint8_t));  

  for(int32_t k=0; odr.read(iv); ++k) {  // read marker
    // periodic message to user
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d. Skipped %d filtered markers and %d uninformative markers, retaining %d variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nskip, nmono, nVariant);

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

    // extract genotype and apply genotype level filter
    if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }
    memset(gtRR, 0, nsamples);
    memset(gtAA, 0, nsamples);
    int32_t ac = 0, an = 0;
    int32_t gcs[3] = {0,0,0};
    for(int32_t i=0; i < nsamples; ++i) {
      int32_t g1 = p_gt[2*i];
      int32_t g2 = p_gt[2*i+1];
      int32_t geno;
      if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
	//geno = 0;
      }
      else {
	geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
	if ( geno == 0 )    { gtRR[i] = 1; }
	else if ( geno == 2 ) { gtAA[i] = 1; }
	ac += geno;
	an += 2;
	++gcs[geno];
      }
    }
    if ( ( gcs[0] < min_hom_gts ) || ( gcs[2] < min_hom_gts ) ) { ++nmono; }
    else {
      bmatRR.add_row_bytes(gtRR);
      bmatAA.add_row_bytes(gtAA);
      nRRs.push_back(gcs[0]);
      nAAs.push_back(gcs[2]);  
      ++nVariant;
    }
  }
  notice("Finished Processing %d markers across %d samples, Skipping %d filtered markers and %d uninformative markers", nVariant, nsamples, nskip, nmono);

  free(gtRR);
  free(gtAA);  

  if ( nVariant < min_variant ) {
    notice("Observed only %d informative markers. Skipping IBD segment detection for this chunk...", nVariant);
    return 0;
  }

  bmatRR.transpose();
  bmatAA.transpose();

  notice("Searching for potential IBD segments..");
  int32_t nibds = 0, k = 0;
  std::vector<int32_t> byte2cnt;
  for(int32_t i=0; i < 256; ++i) {
    int32_t sum = 0;
    for(int32_t j=0; j < 8; ++j) {
      if ( ( (0x00ff & i) >> j ) & 0x01 ) ++sum;
    }
    byte2cnt.push_back(sum);
  }

  //return 0;

  htsFile* wf = hts_open(out.c_str(), "w");
  //hprintf(wf,"ID1\tID2\n");

  for(int32_t i=1; i < nsamples; ++i) {
    uint8_t* iRR = bmatRR.get_row_bits(i);
    uint8_t* iAA = bmatAA.get_row_bits(i);
    for(int32_t j=0; j < i; ++j) {
      uint8_t* jRR = bmatRR.get_row_bits(j);
      uint8_t* jAA = bmatAA.get_row_bits(j);
      //int32_t nibs0 = 0;
      for(k=0; k < bmatRR.nbytes_col; ++k) {
	//nibs0 += byte2cnt[( iRR[k] ^ jRR[k] ) & ( iAA[k] ^ jAA[k] ) & 0x0ff];
	//if ( ( i == 5219 ) && ( j == 65 ) ) {
	//  uint8_t idx = ( iRR[k] ^ jRR[k] ) & ( iAA[k] ^ jAA[k] ) & 0x0ff;
	//  printf("%d %d %02x %02x %02x %02x %02x\n", k, nibs0, iRR[k], iAA[k], jRR[k], jAA[k], idx);
	//}
	if ( ( iRR[k] ^ jRR[k] ) & ( iAA[k] ^ jAA[k] ) ) { // IBS0 exists
	  break;
	}
      }
      if ( k == bmatRR.nbytes_col ) { // no IBS0 observed
	hprintf(wf,"%s\t%s\n",odr.hdr->id[BCF_DT_SAMPLE][i].key, odr.hdr->id[BCF_DT_SAMPLE][j].key);
	++nibds;
      }
    }
  }
  notice("Finished searching for potential IBD segments, identifying %d pairs", nibds);
  hts_close(wf);

  return 0;
}

