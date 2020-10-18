#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_variant_key.h"
#include "bcf_filtered_reader.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "tsv_reader.h"

int32_t cmdVcfFillPlink2(int32_t argc, char** argv) {
  BCFFilteredReader bfr;
  std::string out;
  std::string fill_dosage_field("DS");
  std::string add_af_field("AF");
  std::string add_mach_rsq_field("RSQ");
  std::string add_fimp_field("FIMP");
  double thres_imp = 0.0;
  int32_t verbose = 10000;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input VCFs", NULL)
    LONG_STRING_PARAM("vcf",&bfr.bcf_file_name, "Input VCF/BCF file")

    LONG_PARAM_GROUP("Options to modify input VCF", NULL)
    LONG_STRING_PARAM("region",&bfr.target_region, "Target region to focus on")    
    LONG_STRING_PARAM("fill-dosage",&fill_dosage_field, "Fill dosage from GT field if not exists. Fill zero if missing (with warning)")
    LONG_STRING_PARAM("add-af-field", &add_af_field, "Add allele frequency field from dosage. Must be used with --fill-dosage")
    LONG_STRING_PARAM("add-mach-rsq-field", &add_mach_rsq_field, "Add a field to represent imputation quality (mach-rsq). Must be used with --fill-dosage")
    LONG_STRING_PARAM("add-fimp-field", &add_fimp_field, "Add a field to represent fraction of imperfect imputation (dosage is not exactly 0,1,2). Must be used with --fill-dosage")
    LONG_DOUBLE_PARAM("thres-fimp", &thres_imp, "Threshold of allowance to determine imperfect imputation")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")    
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // sanity check of input arguments
  if ( out.empty() ) {
    error("[E:%s:%d %s] --vcf and --out is a require parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  std::vector<GenomeInterval> intervals;
  if ( !bfr.target_region.empty() ) {
    intervals.push_back( GenomeInterval(bfr.target_region) );
  }

  bfr.init_params();  

  BCFOrderedWriter odw(out.c_str());
  odw.set_hdr(bfr.cdr.hdr);

  char buffer[65535];

  if ( bcf_hdr_id2int(bfr.cdr.hdr, BCF_DT_ID, fill_dosage_field.c_str() ) < 0 ) {
    error("Cannot find FORMAT field %s", fill_dosage_field.c_str());
  }

  if ( bcf_hdr_id2int(odw.hdr, BCF_DT_ID, add_af_field.c_str() ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=%s,Number=1,Type=Float,Description=\"Allele Frequency\">\n",
	    add_af_field.c_str());
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odw.hdr, BCF_DT_ID, add_mach_rsq_field.c_str() ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=%s,Number=1,Type=Float,Description=\"Imputation Quality RSQ\">\n",
	    add_mach_rsq_field.c_str());
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odw.hdr, BCF_DT_ID, add_fimp_field.c_str() ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=%s,Number=1,Type=Float,Description=\"Fraction of imperfectly imputed dosages at tolerance %f\">\n",
	    add_fimp_field.c_str(), thres_imp);
    bcf_hdr_append(odw.hdr, buffer);
  }
  odw.write_hdr();

  int32_t* gts = NULL;
  int32_t n_gts = 0;
  float* dss = NULL;
  int32_t n_dss = 0;

  int32_t ns = bfr.get_nsamples();
  bool no_dosage = false;

  for( int32_t i=0; bfr.read(); ++i ) {
    bcf1_t* iv = bfr.cursor();
    
    //if ( i % verbose == 0 )
    //  notice("Reading %d variants at %s:%d", i, bcf_get_chrom(bfr.cdr.hdr, bfr.cursor()), bfr.cursor()->pos);

    // get GT fields
    if ( bcf_get_genotypes(bfr.cdr.hdr, iv, &gts, &n_gts) < 0 ) {
      error("Cannot find GT field at %s:%d", bcf_get_chrom(bfr.cdr.hdr, bfr.cursor()), bfr.cursor()->pos);
    }

    // get dosages
    if ( bcf_get_format_float(bfr.cdr.hdr, iv, fill_dosage_field.c_str(), &dss, &n_dss) < 0 ) {
      //error("Cannot find %s field at %s:%d", fill_dosage_field.c_str(), bcf_get_chrom(bfr.cdr.hdr, bfr.cursor()), bfr.cursor()->pos);
      no_dosage = true;
      //dss = (float*)calloc(ns, sizeof(float));
    }
    else no_dosage = false;

    bcf1_t* nv = bcf_dup(bfr.cursor());    

    // now check if dosage is missing and fill it with GTs
    double sumds = 0, ssqds = 0;
    int32_t nimp = 0;
    for(int32_t j=0; j < ns; ++j) {
      //printf("%d %d %d %d %f\n", i, j, gts[2*j], gts[2*j+1], dss[j]);
      if ( no_dosage || isnanf(dss[j]) ) {
	if ( ( bcf_gt_is_missing(gts[2*j+1]) ) || ( bcf_gt_is_missing(gts[2*j]) ) ) { // haploid or missing
	  warning("Missing genotype and dosage at %s:%d for individual %d. Assuming zeros", bcf_get_chrom(bfr.cdr.hdr, bfr.cursor()), bfr.cursor()->pos, j);
	  dss[j] = 0.0;
	  ++nimp;
	}
	else {
	  int32_t geno = bcf_gt_allele(gts[2*j])+bcf_gt_allele(gts[2*j+1]);
	  dss[j] = (float)geno;
	  sumds += (double)geno;
	  ssqds += (double)(geno*geno);
	}
      }
      else {
	sumds += (double)dss[j];
	ssqds += (double)(dss[j] * dss[j]);
	if ( ( ( dss[j] > thres_imp ) && ( dss[j] < 1.0-thres_imp ) ) ||
	     ( ( dss[j] > 1.0+thres_imp ) && ( dss[j] < 2.0-thres_imp ) ) ) {
	  ++nimp;
	}
      }
    }

    bcf_update_format_float(odw.hdr, nv, fill_dosage_field.c_str(), dss, n_dss);

    float af = sumds / ns / 2.0;
    double varExp = (double)sumds * ((double)(ns + ns) - sumds) / (double)ns / (double)ns / 2.0;
    double varObs = (double)ssqds / ns - (double)sumds * (double)sumds / (double)ns / (double)ns;
    float machR2 = (float)(varObs / (varExp + 1e-20));
    float fimp = (float)nimp / (float)ns;

    //notice("%d\t%.5f\t%.5f\t%.3lg\t%.3lg\t%.5f", ns, sumds, ssqds, varExp, varObs, machR2);

    bcf_update_info_float(odw.hdr, nv, add_af_field.c_str(), &af, 1);
    bcf_update_info_float(odw.hdr, nv, add_mach_rsq_field.c_str(), &machR2, 1);
    bcf_update_info_float(odw.hdr, nv, add_fimp_field.c_str(), &fimp, 1);

    odw.write(nv);
    bcf_destroy(nv);
    
    //free(gts);
    //free(dss);

    //notice("Reading %d variants at %s:%d", i, bcf_get_chrom(bfr.cdr.hdr, bfr.cursor()), bfr.cursor()->pos);    

    //if ( i > 10 ) break;
  }
  
  odw.close();
  bfr.cdr.close();

  return 0;
}

