#include "bcf_filter_arg.h"
#include "cramore.h"
#include "bcf_ordered_reader.h"
#include "compact_matrix.h"
#include "emPhaser.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>

// goal -- for a given chunk (e.g. 100kb), identify pairs with no IBS0 pairs
int32_t cmdVcfIBS0Phase(int32_t argc, char** argv) {
  std::string inVcf;
  std::string out;
  int32_t min_hom_gts = 1;
  int32_t verbose = 1000;
  std::string reg;
  int32_t seed = 0;

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;
  
  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)    
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("min-hom", &min_hom_gts, "Minimum number of homozygous genotypes to be counted for IBS0")
    LONG_INT_PARAM("seed", &seed, "Random seed")

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

  if ( seed > 0 ) srand(seed);
  else srand(std::time(NULL));

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

  notice("Started Reading site information from VCF file, identifying %d samples", nsamples);

  if ( nsamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

  //std::vector<uint16_t> genos(nsamples,0);
  //std::vector<uint8_t>  haps(nsamples*2,0);
  int bidx = 0;
  int32_t* p_gt = NULL;
  int32_t n_gt = 0;
  int32_t nskip = 0, nmono = 0;

  // read each variant
  hap8phaser_t phaser(nsamples);
  
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
    int32_t ac = 0, an = 0;
    int32_t gcs[3] = {0,0,0};
    std::vector<uint8_t> genos(nsamples,0);
    for(int32_t i=0; i < nsamples; ++i) {
      int32_t g1 = p_gt[2*i];
      int32_t g2 = p_gt[2*i+1];
      uint8_t geno;
      if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
	geno = 0; // missing GT = 0
      }
      else {
	geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0) + 1; // GT = 1,2,3
	ac += (geno-1);
	an += 2;
	++gcs[geno-1];
      }
      genos[i] = geno;
      //phaser.geno8s[i] |= ((geno & 0x03) << (bidx*2)); // 2 bit encoding of genotypes
    }

    if ( ( gcs[0] < min_hom_gts ) || ( gcs[2] < min_hom_gts ) ) continue;
    else {
      for(int32_t i=0; i < nsamples; ++i)
	phaser.geno8s[i] |= ((genos[i] & 0x03) << (bidx*2)); // 2 bit encoding of genotypes
      
      if ( bidx == 7 ) { // block is complete
	//std::copy(phaser.geno8s.begin(), phaser.geno8s.end(), genos);
	phaser.runEMPhase();
	phaser.printFreq();
	phaser.printHaplotypes();
	std::fill(phaser.geno8s.begin(), phaser.geno8s.end(), 0);
	break;
      }
      ++bidx;
    }
  }
  notice("Finished Processing %d markers across %d samples, Skipping %d filtered markers and %d uninformative markers", nVariant, nsamples, nskip, nmono);

  return 0;
}
