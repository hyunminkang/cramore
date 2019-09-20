#ifndef __DSC_LIB_H
#define __DSC_LIB_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <stdint.h>

#include "Error.h"
#include "bcf_filtered_reader.h"


#define MIN_NORM_GL 1e-6

class biallelic_var_t {
public:
  int32_t rid;
  int32_t pos;
  std::string ref;
  std::string alt;
  double af;
  double* gps;
}

double logAdd(double la, double lb);

class biallelic_pileup_t {
  int32_t nreads;
  int32_t nref;
  int32_t nalt;
  double  gls[9];
  double  logdenom;

  biallelic_pileup_t : nreads(0), nref(0), nalt(0) {
    std::fill(gls, gls+9, 1.0);    
    logdenom = 0;
  }

  void merge(const biallelic_pileup_t& other) {
    nreads += other.nreads;
    nref += other.nref;
    nalt += other.nalt;
    logdenom += other.logdenom;
    
    for(int i=0; i < 9; ++i) gls[i] *= other.gls[i];
    
    double tmp = 0;
    for(int i=0; i < 9; ++i) tmp += gls[i];
    
    logdenom += log(tmp);
    
    for(int i=0; i < 9; ++i) gls[i] /= tmp;
    
    for(int i=0; i < 9; ++i) {
      if ( gls[i] < MIN_NORM_GL ) {
	gls[i] = MIN_NORM_GL;
      }
    }
    tmp = 0;
    for(int i=0; i < 9; ++i) tmp += gls[i];
    logdenom += log(tmp);
    for(int i=0; i < 9; ++i) gls[i] /= tmp;    
  }

  void remove(const snp_droplet_pileup& other) {
    nreads -= other.nreads;
    nref -= other.nref;
    nalt -= other.nalt;
    logdenom -= other.logdenom;
    
    for(int i=0; i < 9; ++i) gls[i] /= other.gls[i];
    
    double tmp = 0;
    for(int i=0; i < 9; ++i) tmp += gls[i];
    
    logdenom += log(tmp);
    
    for(int i=0; i < 9; ++i) gls[i] /= tmp;
    
    for(int i=0; i < 9; ++i) {
      if ( gls[i] < MIN_NORM_GL ) {
	gls[i] = MIN_NORM_GL;
      }
    }
    tmp = 0;
    for(int i=0; i < 9; ++i) tmp += gls[i];
    logdenom += log(tmp);
    for(int i=0; i < 9; ++i) gls[i] /= tmp;    
  }  
};

class dsc_lib_t {
public:
  int32_t n_tmp_files;
  std::string n
}

class sc_dropseq_lib_t {
 public:
  // vector containing SNP & genotype info, index is snp_id
  std::map<std::string, int32_t> chr2rid;
  std::vector<std::string>       rid2chr;
  
  std::vector<sc_snp_t> snps;

  // mapper between barcode -> bcd_id  
  std::map<std::string,int32_t> bc_map;
  std::vector<std::string> bcs;
  
  // cell_umis[i]->[j] contains the map of UMIs overlapping with snp j in cell i
  std::vector< std::map<int32_t,sc_snp_droplet_t*> > cell_umis;

  // Number of pass-filtered reads and unique reads
  std::vector<int32_t> cell_totl_reads;  
  std::vector<int32_t> cell_pass_reads;
  std::vector<int32_t> cell_uniq_reads;
  std::vector<double>  cell_scores;
  
  std::vector< std::map<int32_t,sc_snp_droplet_t*> > snp_umis;
  int32_t nbcs;
  int32_t nsnps;
  int32_t add_snp(int32_t _rid, int32_t _pos, char _ref, char _alt, double _af, double* _gps);
  int32_t add_cell(const char* barcode);
  bool add_read(int32_t snpid, int32_t cellid, const char* umi, char allele, char qual);

  int32_t load_from_plp(const char* plpPrefix, BCFFilteredReader* pvr = NULL, const char* field = "GP", double genoErrorOffset = 0.1, double genoErrorCoeffR2 = 0.0, const char* r2info = "R2", bool loadUMI = false);

  sc_dropseq_lib_t() : nbcs(0), nsnps(0) {}

  dropD calculate_droplet_clust_distance(std::map<int32_t,snp_droplet_pileup*> dropletPileup, std::map<int32_t,snp_droplet_pileup>& clustPileup);
  
};

double calculate_snp_droplet_GL(sc_snp_droplet_t* ssd, double* gls);
double calculate_snp_droplet_doublet_GL(sc_snp_droplet_t* ssd, double* gls, double alpha);
double calculate_snp_droplet_pileup(sc_snp_droplet_t* ssd, snp_droplet_pileup* sdp, double alpha);

struct sc_drop_comp_t {
  sc_dropseq_lib_t* pscl;
  sc_drop_comp_t(sc_dropseq_lib_t* p) : pscl(p) {}
  bool operator()(const int32_t& lhs, const int32_t& rhs) const {
    double cmp = pscl->cell_scores[lhs] - pscl->cell_scores[rhs];
    if ( cmp != 0 ) return cmp > 0;
    else return lhs > rhs;
  }
};

#endif
