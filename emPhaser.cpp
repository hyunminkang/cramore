#include "emPhaser.h"
#include "Error.h"

#include <cmath>
#include <cstring>

hap8map_t hap8map;

#define BYTE_TO_BINARY_PATTERN "%c%c%c%c%c%c%c%c"
#define BYTE_TO_BINARY(byte)  \
  (byte & 0x80 ? '1' : '0'), \
  (byte & 0x40 ? '1' : '0'), \
  (byte & 0x20 ? '1' : '0'), \
  (byte & 0x10 ? '1' : '0'), \
  (byte & 0x08 ? '1' : '0'), \
  (byte & 0x04 ? '1' : '0'), \
  (byte & 0x02 ? '1' : '0'), \
  (byte & 0x01 ? '1' : '0')

// constructor
hap8map_t::hap8map_t() {
  // possible combinations of phase
  for(int32_t i=0; i < 256; ++i) {
    for(int32_t j=0; j <= i; ++j) {
      uint16_t hap2 = ((uint16_t)j << 8) | (uint16_t)i;
      uint16_t gen = 0;
      for(int32_t k=0; k < 8; ++k) {
	gen |= ((( (j >> k) & 0x01) + ((i >> k) & 0x01) + 1) << (k+k));
      }
      //fprintf(stderr,"%x\n",gen);
      //geno2hap2s[gen].push_back(hap2); // assume no missingness

      // the number of possible cases of missingness is 2^8 = 256
      for(int32_t k=0; k < 256; ++k) {
	uint16_t mgen = ( ( ( (k >> 0) & 0x01 ) ? ( 0x00 << 0 ) : ( ( ( gen >> 0 ) & 0x03 ) << 0 ) ) |
			  ( ( (k >> 1) & 0x01 ) ? ( 0x00 << 2 ) : ( ( ( gen >> 2 ) & 0x03 ) << 2 ) ) |
			  ( ( (k >> 2) & 0x01 ) ? ( 0x00 << 4 ) : ( ( ( gen >> 4 ) & 0x03 ) << 4 ) ) |
			  ( ( (k >> 3) & 0x01 ) ? ( 0x00 << 6 ) : ( ( ( gen >> 6 ) & 0x03 ) << 6 ) ) |
			  ( ( (k >> 4) & 0x01 ) ? ( 0x00 << 8 ) : ( ( ( gen >> 8 ) & 0x03 ) << 8 ) ) |
			  ( ( (k >> 5) & 0x01 ) ? ( 0x00 << 10 ) : ( ( ( gen >> 10 ) & 0x03 ) << 10 ) ) |
			  ( ( (k >> 6) & 0x01 ) ? ( 0x00 << 12 ) : ( ( ( gen >> 12 ) & 0x03 ) << 12 ) ) |
			  ( ( (k >> 7) & 0x01 ) ? ( 0x00 << 14 ) : ( ( ( gen >> 14 ) & 0x03 ) << 14 ) ) );
	geno2hap2s[mgen].push_back(hap2);
      }
    }
  }

  for(int32_t i=0; i < 65536; ++i) {
    for(int32_t k=0; k < 8; ++k) {
      geno8to1s[i][k] = ( ( i >> (k+k) ) & 0x03 );
    }
  }
}

// constructor
hap8phaser_t::hap8phaser_t(int32_t _ns, double _alpha) : nsamples(_ns), alpha(_alpha) {
  geno8s.resize(nsamples,0);
  //hap8s.resize(2*nsamples,0);
  //resolved.resize(nsamples,false);
  //acs.resize(8);
  //ans.resize(8);
}

// perform an EM iteration
double hap8phaser_t::emIteration() {
  int32_t nchanged = 0;
  double tmpUnresolved[256];
  memset(tmpUnresolved, 0, sizeof(double)*256);

  for(auto it = idxUnresolved.begin(); it != idxUnresolved.end(); ++it) {
    std::vector<uint16_t>& hap2 = hap8map.geno2hap2s[geno8s[*it]];
    int nh = hap2.size();
    std::vector<double> frqs(nh,0);
    double sumfreq = 0;
    for(int32_t j=0; j < nh; ++j) {
      uint8_t h1 = ( (hap2[j] >> 8) & 0xff ); // h1 and h2 cannot be the same	
      uint8_t h2 = (hap2[j] & 0xff);
      frqs[j] = hap8Freq[h1] * hap8Freq[h2];
      sumfreq += frqs[j];
    } // normalied frequency calculated

    for(int32_t j=0; j < nh; ++j) {
      uint8_t h1 = ( (hap2[j] >> 8) & 0xff ); // h1 and h2 cannot be the same	
      uint8_t h2 = (hap2[j] & 0xff);
      double r = frqs[j] / sumfreq;
      tmpUnresolved[h1] += r;
      tmpUnresolved[h2] += r;
    }
  }

  double diff = 0, tmp;
  for(int32_t i=0; i < 256; ++i) {
    tmp = ( hap8countsResolved[i] + tmpUnresolved[i] ) / ( nsamples + nsamples );
    diff += fabs(tmp - hap8Freq[i]);
    hap8Freq[i] = tmp;
  }
  return diff;
}

// initialize phasing
void hap8phaser_t::initPhase() {
  //std::fill(acs.begin(), acs.end(), 0);
  //std::fill(ans.begin(), ans.end(), 0);
  //std::fill(resolved.begin(), resolved.end(), false);
  memset(hap8countsResolved, 0, sizeof(uint32_t)*256); // initialize haplotype counts
  //memset(hap8countsUnresolved, 0, sizeof(double)*256);
  idxUnresolved.clear();                               // indices of unresolved haplotypes

  // calculate AC, AN, resolved haplotype frequencies
  for(int32_t i=0; i < nsamples; ++i) {
    uint16_t g = geno8s[i];
    
    /*
    for(int32_t k=0; k < 8; ++k) {
      switch(hap8map.geno8to1s[g][k]) {
      case 3: ++acs[k];
      case 2: ++acs[k];
      case 1: ++acs[k]; ans[k] +=2;
      }
    }
    */
    
    std::vector<uint16_t>& hap2 = hap8map.geno2hap2s[g]; // check possible haplotypes
    //int32_t nresolved = 0;
    if ( hap2.empty() ) { abort(); } // should not happen without missing genotypes
    else if ( hap2.size() == 1 ) {
      //resolved[i] = true;
      //++nresolved;
      uint8_t h1 = ( (hap2[0] >> 8) & 0xff );      
      uint8_t h2 = ( hap2[0] & 0xff);
      //hap8s[i+i] = h1;
      //hap8s[i+i+1] = h2;
      ++hap8countsResolved[h1];
      ++hap8countsResolved[h2];
    }
    else {
      idxUnresolved.push_back(i);
    }
  }

  // resolve unresolved haplotypes, partially trusting resolved counts
  int32_t nhResolved = 2 * (nsamples - idxUnresolved.size());
  memset(hap8Freq, 0, sizeof(double)*256);
  memset(hap8FreqU, 0, sizeof(double)*256);   
  for(auto it = idxUnresolved.begin(); it != idxUnresolved.end(); ++it) {
    std::vector<uint16_t>& hap2 = hap8map.geno2hap2s[geno8s[*it]];
    int nh = hap2.size();
    std::vector<double> frqs(nh,0);
    double sumfreq = 0;
    for(int32_t j=0; j < nh; ++j) {
      uint8_t h1 = ( (hap2[j] >> 8) & 0xff ); // h1 and h2 cannot be the same	
      uint8_t h2 = (hap2[j] & 0xff);
      frqs[j] = (1.0-alpha) * hap8countsResolved[h1] * hap8countsResolved[h2] / nhResolved / nhResolved + alpha / nh / nh; 
      sumfreq += frqs[j];
    } // normalzied frequency calculated

    //double p = (rand() + 0.5) / (RAND_MAX + 1.), q = 0, r = 0;
    for(int32_t j=0; j < nh; ++j) {
      uint8_t h1 = ( (hap2[j] >> 8) & 0xff ); // h1 and h2 cannot be the same	
      uint8_t h2 = (hap2[j] & 0xff);
      double r = frqs[j] / sumfreq;
      hap8Freq[h1] += r;
      hap8Freq[h2] += r;
      hap8FreqU[h1] += (1.0/nh);
      hap8FreqU[h2] += (1.0/nh);      
      //if ( ( q < p ) && ( q + r >= p ) ) { // randomly select haplotypes
      //hap8s[i+i]   = h1;
      //hap8s[i+i+1] = h2;
      //	break;
      //}
      //q += r;
      //}
    }
  }

  for(int32_t i=0; i < 256; ++i) {
    hap8Freq[i] = ( hap8Freq[i] + hap8countsResolved[i] ) / ( nsamples + nsamples );
    hap8FreqU[i] = ( hap8FreqU[i] + hap8countsResolved[i] ) / ( nsamples + nsamples );    
  }
}

int32_t hap8phaser_t::runEMPhase(int32_t maxiter, double thres) {
  // get the initial phasing
  notice("Initializing the phase...");
  initPhase();
  notice("Resolved = %d/%d (%.5lf)", (nsamples+nsamples-idxUnresolved.size()), nsamples+nsamples, 1.0 - (double)idxUnresolved.size()/(nsamples+nsamples));

  for(int32_t i=0; i < maxiter; ++i) {
    double diff = emIteration();
    notice("Performed EM iteration %d, diff = %.3le", i+1, diff);
    if ( diff < thres ) {
      notice("Difference is negligible.. Finishing early");
      return i+1;
    }
  }
  notice("Finished all %d iterations", maxiter);
  return maxiter;
}

void hap8phaser_t::printFreq() {
  //uint32_t hCnts[256];
  //memset(hCnts, 0, sizeof(uint32_t)*256);
  //for(int32_t i=0; i < nsamples+nsamples; ++i) 
  //  ++hCnts[hap8s[i]];

  for(int32_t i=0; i < 256; ++i) {
    uint32_t cntR = hap8countsResolved[i];
    double   freq = hap8Freq[i];
    double   freqU = hap8FreqU[i];    
    if ( ( cntR > 0 ) || ( freqU > 0 ) ) {
    notice("" BYTE_TO_BINARY_PATTERN "\t%u\t%.1le\t%.1le\t%.2lf", BYTE_TO_BINARY(i), cntR, freq, freqU, freqU*nsamples*2);
    }
  }
}

void hap8phaser_t::printHaplotypes() {
  for(int32_t i=0; i < nsamples; ++i) {
    std::vector<uint16_t>& hap2s = hap8map.geno2hap2s[geno8s[i]];
    if ( hap2s.size() == 1 ) {
      printf("%d\t%04x\t1\t%04x,1.0e+00\n", i, geno8s[i], hap2s[0]); 
    }
    else {
      double sumFreq = 0;
      int32_t nh = (int32_t) hap2s.size();
      std::vector<double> frqs(nh,0);      
      for(int32_t j=0; j < nh; ++j) {
	uint8_t h1 = (( hap2s[j] >> 8 ) & 0xff);
	uint8_t h2 = (hap2s[j] & 0xff);
	sumFreq += ( frqs[j] = hap8Freq[h1] * hap8Freq[h2] );
      }
      printf("%d\t%04x\t%d", i, geno8s[i], nh);
      for(int32_t j=0; j < nh; ++j) {
	printf("\t%04x,%.1le", hap2s[j], frqs[j]/sumFreq);
      }
      printf("\n");
    }
  }
}
