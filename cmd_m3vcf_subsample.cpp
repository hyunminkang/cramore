#include "m3vcf.h"
#include "tsv_reader.h"
#include "cramore.h"

#include <set>
#include <map>
#include <cstring>

// Software tool for subsetting M3VCF file
int32_t cmdM3vcfSubsample(int32_t argc, char** argv) {
  std::string in_m3vcf;
  std::string samplef;
  std::string out_m3vcf;
  bool removeMono = false;

  //// Parameter Handling
  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Required Input Files", NULL)
    LONG_STRING_PARAM("m3vcf", &in_m3vcf, "Input M3VCF file")
    LONG_STRING_PARAM("samples", &samplef, "File containing sample IDs")

    LONG_PARAM_GROUP("Required Output Files", NULL)    
    LONG_STRING_PARAM("out", &out_m3vcf, "Output M3VCF file")    

    LONG_PARAM_GROUP("Additional Options", NULL)    
    LONG_PARAM("remove-mono", &removeMono, "Remove monomorphic variants") 
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( in_m3vcf.empty() || samplef.empty() || out_m3vcf.empty() ) {
    error("[E:%s:%d %s] --m3vcf, --samples, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  tsv_reader tr_samples(samplef.c_str());
  std::set<std::string> subset_ids;
  
  while( tr_samples.read_line() > 0 ) {
    if ( subset_ids.insert(tr_samples.str_field_at(0)).second == false ) // there are duplicates
      error("Observed duplicated sample ID %s", tr_samples.str_field_at(0));
  }

  m3vcf_reader_t m3r;
  m3r.open(in_m3vcf.c_str()); // open M3VCF file

  // check the overlaps between the subset_ids
  std::map<std::string,int32_t> id2idx;
  for(int32_t i=0; i < (int32_t)m3r.m3vcf.inds.size(); ++i) {
    id2idx[m3r.m3vcf.inds[i]] = i;
  }

  std::vector<int32_t> idxs;
  std::vector<std::string> ids;
  for(std::set<std::string>::iterator it = subset_ids.begin(); it != subset_ids.end(); ++it) {
    std::map<std::string,int32_t>::iterator it2 = id2idx.find(*it);
    if ( it2 == id2idx.end() )
      error("Sample ID %s not found from M3VCF", it->c_str());
    idxs.push_back(it2->second);
    ids.push_back(*it);
  }

  m3vcf_t m3sub; // new M3VCF to be used for subsetting
  m3sub.inds = ids;
  m3sub.hdrs = m3r.m3vcf.hdrs;

  // open the output file
  htsFile* wf = hts_open(out_m3vcf.c_str(), "wg");
  m3sub.write_hdrs(wf); // write the header

  for(int32_t i=0; m3r.read_block() > 0; ++i) {
    m3vcf_block_t* p_orig = m3r.m3vcf.blocks.back();
    if ( i % 100 == 0 )
      notice("Processing %d blocks at %s:%d-%d", i+1, p_orig->chrom.c_str(), p_orig->pos_beg, p_orig->pos_end);
    
    m3vcf_block_t* p_block = p_orig->create_subset(idxs, removeMono);
    m3sub.add_block(p_block);
    m3sub.write_block(wf);
  }
  hts_close(wf); // close the output file

  return 0;
}
