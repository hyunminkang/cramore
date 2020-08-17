#ifndef __M3VCF_H
#define __M3VCF_H

#include <cstring>
#include "tsv_reader.h"

// class to store variant site info. can be used for VCF too.
class m3vcf_site_t {
 public:
  std::string chrom;
  int32_t pos;
  std::string varid;
  std::string ref;
  std::string alt;
  std::string qual;
  std::string filt;
  std::string info;

 m3vcf_site_t(const char* _chrom, int32_t _pos, const char* _varid, const char* _ref, const char* _alt, const char* _qual, const char* _filt, const char* _info) : chrom(_chrom), pos(_pos), varid(_varid), ref(_ref), alt(_alt), qual(_qual), filt(_filt), info(_info) {}
};


// class to represent M3VCF block
class m3vcf_block_t {
 public:
  std::string chrom;  // chromosome
  int32_t pos_beg;    // beginning position of the block
  int32_t pos_end;    // ending position of the block
  int32_t nsamples;   // number of samples assumes diploids
  int32_t nreps;      // number of unique haplotypes.
  std::vector<int32_t> ireps;      // [0..(nhaps-1)]
  std::vector<m3vcf_site_t> sites; // site information for each variant
  std::vector< std::vector<int32_t> > variants; // diff-index of non-ref carriers of the haplotypes

  m3vcf_block_t() : nsamples(-1), nreps(-1) {}  // default constructor
  m3vcf_block_t(const std::vector<int32_t>& ihaps) { init(ihaps); } // constructor with ireps specified
  void    init(const std::vector<int32_t>& ihaps);    // ihaps should be a size of 2*nsamples
  int32_t add_variant(const m3vcf_site_t& site, const std::vector<int32_t>& inrefs); // haps should be a size of nhaps

  bool write(htsFile* wf);  // write the block into a file
  
  m3vcf_block_t* create_subset(const std::vector<int32_t>& isamples, bool removeMono = true); // create a new block in the heap based on the sample indices
};

class m3vcf_t {
 public:
  int32_t nblocks;     // number of total blocks read
  bool keep_in_memory;  // if true, blocks are not replaced, but appeneded
  std::vector<std::string> inds;      // sample IDs
  std::vector<std::string> hdrs;      // header lines;
  std::vector<m3vcf_block_t*> blocks; // blocks in memory

  m3vcf_t() : nblocks(0), keep_in_memory(false) {} // default constructor 
  ~m3vcf_t() { // destructor to explictly delete the blocks in the heap
    for(int32_t i=0; i < (int32_t)blocks.size(); ++i) {
      delete blocks[i];
    }
  }

  bool write_hdrs(htsFile* wf); // write headers to file, return true if succeed
  bool write_block(htsFile* wf, int32_t idx = 0); // write a block to file, return true if succeed
  bool add_block(m3vcf_block_t* new_block); // add or replace a block. return true if added false if replaced existing
};


class m3vcf_reader_t {
 public:
  //std::string filename;
  m3vcf_t m3vcf;
  tsv_reader tr;

  int32_t open(const char* filename); // open a m3vcf file
  int32_t read_block();  // read a block
};


#endif // __M3VCF_H
