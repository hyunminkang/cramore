#include "compact_matrix.h"
#include "Error.h"

bitmatrix::bitmatrix(int32_t _ncol, int32_t _nrow_alloc) : nrow(0), ncol(_ncol), nrow_alloc(_nrow_alloc) {
  int32_t ncol_ceil8 = (ncol % 8 == 0) ? ncol : ( ncol + ( 8 - (ncol % 8) ) );  
  nbytes_col = ncol_ceil8 >> 3;
  bytes = (uint8_t*)calloc(nrow_alloc * nbytes_col, sizeof(uint8_t));
}

int32_t bitmatrix::reserve(int32_t new_nrow_alloc) {
  if ( new_nrow_alloc == 0 )
    new_nrow_alloc = 2 * nrow_alloc;

  //notice("reserve() called, from %d to %d", nrow_alloc, new_nrow_alloc);  

  uint8_t* new_bytes = (uint8_t*)calloc(new_nrow_alloc * nbytes_col, sizeof(uint8_t));
  memcpy(new_bytes, bytes, nrow * nbytes_col);
  free(bytes);
  bytes = new_bytes;
  nrow_alloc = new_nrow_alloc;
  
  return new_nrow_alloc;
}

// all nonzeros will be considered as 1
int32_t bitmatrix::add_row_ints(int32_t* intarray) {
  int32_t i, j;
  uint8_t byte;
  if ( nrow == nrow_alloc) reserve();
  uint8_t* rowbytes = get_row_bits(nrow);  
  for(i=0; i < ncol; i += 8) {
    for(j=0, byte = 0; (j < 8) && (j < ncol); ++j) {
      byte |= ((intarray[i+j] != 0) << (7-j));
    }
    rowbytes[i >> 3] = byte;
  }
  return ++nrow;
}

int32_t bitmatrix::add_row_bytes(uint8_t* bytearray) {
  int32_t i, j;
  uint8_t byte;
  if ( nrow == nrow_alloc) reserve();
  uint8_t* rowbytes = get_row_bits(nrow);  
  for(i=0; i < ncol; i += 8) {
    for(j=0, byte = 0; (j < 8) && (j < ncol); ++j) {
      byte |= ((bytearray[i+j] != 0) << (7-j));
    }
    rowbytes[i >> 3] = byte;
  }
  return ++nrow;  
}

int32_t bitmatrix::add_row_bits(uint8_t* bitarray) {
  if ( nrow == nrow_alloc) reserve();
  uint8_t* rowbytes = get_row_bits(nrow);
  memcpy(rowbytes, bitarray, nbytes_col);
  return ++nrow;
}

bool bitmatrix::transpose() {
  notice("Transposing the bit matrix...");

  int32_t nrow_ceil8 = (nrow % 8 == 0) ? nrow : ( nrow + ( 8 - (nrow % 8) ) );  
  int32_t nbytes_row = nrow_ceil8 >> 3;
  uint8_t* tbytes = (uint8_t*)calloc(ncol * nbytes_row, sizeof(uint8_t));
  for(int32_t i=0; i < nrow; ++i) {
    uint8_t* rowbytes = get_row_bits(i);
    for(int32_t j=0; j < ncol; ++j) {
      tbytes[j * nbytes_row + i/8] |= ( ( ( rowbytes[j/8] >> (7 - (j%8)) ) & 0x01 ) << ( 7 - i%8 ) );
    }
  }
  notice("Finished transposing the bit matrix...");

  nbytes_col = nbytes_row;
  int32_t tmp = nrow;
  nrow = ncol;
  ncol = tmp;
  nrow_alloc = ncol;

  free(bytes);
  bytes = tbytes;
  
  return true;
}

void bitmatrix::print(int32_t rbeg, int32_t rend, int32_t cbeg, int32_t cend) {
  for(int32_t i=rbeg; i < rend; ++i) {
    uint8_t* row = get_row_bits(i);
    for(int32_t j=cbeg; j < cend; ++j) {
      printf("%02x ",row[j]);
    }
    printf("\n");
  }
}

byte_pair_operation byte_pair_op;
