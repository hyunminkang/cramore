#include "bcf_filter_arg.h"
#include "cramore.h"
#include "bcf_ordered_reader.h"
#include "compact_matrix.h"
#include "emPhaser.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>

// goal -- for a given chunk (e.g. 100kb), identify pairs with no IBS0 pairs
int32_t cmdVcfRefPhase(int32_t argc, char** argv) {
  std::string refVcf;
  std::string targetVcf;
  std::string outVcf;
  
  int32_t min_hom_gts = 1;   // threshold to be considered as 'common' variants
  int32_t min_ac_common = 1; // AC threshold to be considered as 'common' variants
  int32_t min_af_common = 0; // AF threshold to be considered as 'common' variants
  int32_t verbose = 1000;
  std::string reg;
  int32_t seed = 0;

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  /////////////////////////////////
  //// Parameter Handling
  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("ref",&refVcf, "VCF/BCF file of reference haplotype")
    LONG_STRING_PARAM("target",&targetVcf, "VCF/BCF file of targets")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)    
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("min-hom", &min_hom_gts, "Minimum number of homozygous genotypes to be considered as common")
    LONG_INT_PARAM("min-ac", &min_ac_common, "Minimum allele count to be considered as common")
    LONG_INT_PARAM("min-af", &min_ac_common, "Minimum allele frequencies to be considered as common")    
    LONG_INT_PARAM("seed", &seed, "Random seed")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outVcf, "Output VCF/BCF file")    
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")    
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( refVcf.empty() || targetVcf.empty() || outVcf.empty() ) {
    error("[E:%s:%d %s] --ref, --target, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // set random seed
  if ( seed > 0 ) srand(seed);
  else srand(std::time(NULL));

  /////////////////////////////////////
  //// Prepare processing reference VCF
  
  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }
  BCFOrderedReader odr(refVcf, intervals);
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

  //////////////////////////////////////////
  //// Read reference VCFs and fill in the
  //// haplotype and frequency matrix

  int32_t nVariant = 0;
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  int32_t nhaps = nsamples + nsamples;

  notice("Started Reading site information from VCF file, identifying %d samples", nsamples);

  if ( nsamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

  // data structure to store haplotypes and its frequencies
  std::vector<uint8_t*>  refCommon; // reference haplotypes for common variants (8-SNP blocks)
  std::vector<uint32_t*> cntCommon; // haplotype counts of common variants (8-SNP blocks)
  std::vector< std::vector< std::set<uint32_t> > > idxRare; // index of haplotype carriers, grouped by 8-SNP block

  int32_t ncommon = 0, nrare = 0;
  int32_t* p_gt = NULL;
  int32_t n_gt = 0;
  int32_t nskip = 0, nmono = 0;
  uint8_t* bytes = (uint8_t*) malloc(nhaps);

  for(int32_t k=0; odr.read(iv); ++k) {  // read marker
    // periodic message to user
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d. Skipped %d filtered markers and %d uninformative markers, processing %d common and %d rare variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nskip, nmono, ncommon, nrare);

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

    if ( n_gt != nhaps )
      error("[E:%s:%d %s] Number of haplotypes %d is different from expected %d",__FILE__,__LINE__,__FUNCTION__, n_gt, nhaps);
    
    int32_t ac = 0;
    int32_t gcs[2] = {0,0}; // RA, AA
    int32_t cidx = ncommon % 8;
    int32_t ridx = 0;
    for(int32_t i=0; i < n_gt; ++i) {
      int32_t g = p_gt[i];

      //if ( bcf_gt_allele(g) > 0 )
      //fprintf(stderr, "%d %d %d %d %d\n", iv->pos+1, i, bcf_gt_is_missing(g), bcf_gt_is_phased(g), bcf_gt_allele(g));

      if ( bcf_gt_is_missing(g) ) error("Missing genotype is observed at index %d", i);
      if ( ( i % 2 == 1 ) && ( ! bcf_gt_is_phased(g) ) ) error("Unphased haplotype is observed at index %d", i);

      if ( bcf_gt_allele(g) > 0 ) {
	bytes[i] = ( 0x01 << cidx );
	++ac;
	if ( i % 2 == 1 ) ++gcs[bcf_gt_allele(p_gt[i-1])];
      }
      else {
	bytes[i] = 0;
      }
    }
    //if ( iv->pos % 100 == 0 ) break;

    if ( ac == 0 ) { ++nmono; continue; }  // skip monomorphic variants

    // check if common or rare
    if ( ( gcs[1] < min_hom_gts ) || ( ac < min_ac_common ) || ( ac < min_af_common * nhaps ) ) { // rare
      // store indices of carriers
      if ( idxRare.size() <= ncommon / 8 ) {
	idxRare.resize(ncommon/8+1); // need to create another block
	ridx = 0;
      }

      std::vector< std::set<uint32_t> > & ir = idxRare.back();

      ir.resize(ridx+1);  // add a rare variant within the block
      for(int32_t i=0; i < nhaps; ++i) {
	if ( bytes[i] ) ir[ridx].insert(i);
      }
      ++nrare;
    }
    else { // common
      if ( cidx == 0 ) { // create a new block
	refCommon.push_back( (uint8_t*)  calloc( nhaps, sizeof(uint8_t) ) );
	cntCommon.push_back( (uint32_t*) calloc( 256, sizeof(uint32_t) ) );
      }
      
      uint8_t* haps = refCommon.back();
      for(int32_t i=0; i < nhaps; ++i) {
	haps[i] |= bytes[i];
      }

      if ( cidx == 7 ) { // end of the block, calculate haplotype frequency
	uint32_t* cnts = cntCommon.back();
	for(int32_t i=0; i < nhaps; ++i) 
	  ++cnts[haps[i]];
      }
      ++ncommon;
    }
  }
  notice("Finished Processing %d common and %d rare variants markers across %d samples, Skipping %d filtered markers and %d uninformative markers", ncommon, nrare, nsamples, nskip, nmono);

  return 0;
}
