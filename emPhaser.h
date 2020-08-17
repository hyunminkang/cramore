#ifndef __EM_PHASER_H__
#define __EM_PHASER_H__

#include <stdint.h>
#include <vector>
#include <set>

class hap8map_t {
public:
  std::vector<uint16_t> geno2hap2s[65536];
  uint8_t geno8to1s[65536][8];
  hap8map_t();
};

extern hap8map_t hap8map;

// phasing 8 SNPs together
class hap8phaser_t {
public:
  int32_t nsamples;
  double  alpha; // initial randomness parameter : 0 - observed counts only, 1 - random assignments
  hap8phaser_t(int32_t _ns, double _alpha = 0.2);

  std::vector<uint16_t> geno8s;
  uint32_t hap8countsResolved[256]; // counts of observed haplotypes only
  double   hap8FreqU[256];          // uniformly assigned frequency of possible haplotypes
  double   hap8Freq[256];           // frequency of all haplotypes
  //std::vector<uint8_t> hap8s; // predicted haplotypes
  //std::vector<bool> resolved;
  std::vector<int32_t> idxUnresolved;
  //std::vector<int32_t> acs;
  //std::vector<int32_t> ans;

  void initPhase();
  double emIteration();
  int32_t runEMPhase(int32_t maxiter = 20, double thres = 1e-8);
  void printFreq();
  void printHaplotypes();
};

class rareblock_t {
public:
  int32_t nsamples;
  std::vector< std::set<int32_t> > var2ind;
  std::vector< std::set<int32_t> > ind2var;
  std::vector<int32_t> acs;

  rareblock_t(int32_t _ns) : nsamples(_ns) {}

  int32_t addVariant(int32_t ac) {
    acs.push_back(ac);
    return (int32_t)acs.size();
  }

  void add(int32_t iind, int32_t ivar = -1) {
    if ( ivar < 0 ) { ivar = (int32_t)acs.size() - 1; }
    var2ind[ivar].insert(iind);
    ind2var[iind].insert(ivar);
  }
};

#endif // EM_PHASER_H
