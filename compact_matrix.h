#ifndef __COMPACT_MATRIX_H
#define __COMPACT_MATRIX_H

#include <stdint.h>
#include <stddef.h>
#include <cstdlib>
#include <cstring>

class byte_pair_operation {
public:
  uint8_t byte_pair_xor[65536];  
  byte_pair_operation() {
    int32_t i, j, k;
    for(i=0; i < 256; ++i) {
      for(j=0; j < 256; ++j) {
	uint8_t x = ((uint8_t)i ^ (uint8_t)j) & 0x0ff;
	uint8_t s = 0;
	for(k = 0; k < 8; ++k) { s += ((x >> k) & 0x01); }
	byte_pair_xor[i*256+j] = s;
      }
    }
  }

  inline int32_t sum_xor(uint8_t* x, uint8_t* y, int32_t nbytes) {
    int32_t sum = 0;
    for(int32_t i=0; i < nbytes; ++i) sum += byte_pair_xor[(x[i] << 8) + y[i]];
    return sum;
  }

  inline bool positive_xor(uint8_t* x, uint8_t* y, int32_t nbytes) {
    for(int32_t i=0; i < nbytes; ++i) {
      if ( byte_pair_xor[(x[i] << 8) + y[i]] ) return true;
    }
    return false;
  }  
};

// row-wise stored 1-bit matrix
class bitmatrix {
public:
  int32_t nrow;
  int32_t ncol;
  int nrow_alloc;
  int32_t nbytes_col; // internal values
  uint8_t* bytes;
  bitmatrix(int32_t _ncol, int32_t _nrow_alloc = 100);
  ~bitmatrix() { if ( bytes != NULL ) free(bytes); }
  int32_t reserve(int32_t new_nrow_alloc = 0);
  int32_t add_row_ints(int32_t* intarray);
  int32_t add_row_bytes(uint8_t* bytearray);
  int32_t add_row_bits(uint8_t* bitarray);
  bool transpose();
  inline uint8_t* get_row_bits(int32_t irow) { return bytes + irow*nbytes_col; }
  inline uint8_t get_byte_at(int32_t irow, int32_t icol) { return ( ( bytes[irow*nbytes_col + (icol >> 3)] >> (icol & 0x03) ) & 0x01 ); }
  //void head(int32_t r, int32_t c);
  void print(int32_t rbeg, int32_t rend, int32_t cbeg, int32_t cend);
};

extern byte_pair_operation byte_pair_op;

#endif // __COMPACT_MATRIX_H
