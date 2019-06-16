#ifndef __DSC_DGE_H
#define __DSC_DGE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <stdint.h>

#include "Error.h"

class dsc_dge {
public:
  std::vector<std::string> bcds;
  std::vector<std::string> gene_ids;
  std::vector<std::string> gene_symbs;
  std::map<std::string,int32_t> bcd2idx;
  std::map< int32_t,std::map<int32_t,int32_t> > bcd_gene_umi;
  std::map< int32_t,std::map<int32_t,int32_t> > gene_bcd_umi;
  std::vector<int32_t> bcd_umis;
  std::vector<int32_t> gene_umis;
  
  dsc_dge(const char* bcdf, const char* genef, const char* mtxf);
};

#endif
