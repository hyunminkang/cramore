#include "cramore.h"
#include "tsv_reader.h"
#include <ctime>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>

/////////////////////////////////////////////////////////////////////////
// dge-shuffle : Shuffle digital expression matrix
////////////////////////////////////////////////////////////////////////
int32_t cmdDgeShuffle(int32_t argc, char** argv) {
  std::string mtxf("matrix.mtx");
  std::string outf;
  int32_t seed = 0;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("mtx",&mtxf,   "File containing sparse matrices of gene, barcode, UMI")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outf,"Output matrix file")
    LONG_INT_PARAM("seed",&seed,"Random seed")    
    LONG_INT_PARAM("verbose", &globalVerbosityThreshold, "Turn on verbose mode with specific verbosity threshold. 0: fully verbose, 100 : no verbose messages")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( mtxf.empty() || outf.empty() )
    error("Missing required option(s) : --mtx, --out");

  notice("Analysis Started");

  notice("Reading %s",mtxf.c_str());    
  tsv_reader tsv_mtxf(mtxf.c_str());
  int64_t n_bcds = -1, n_genes = -1, sum_umis = 0, num_nonzeros = -1;
  std::vector<int32_t> igenes;
  std::vector<int32_t> ibcds;
  int32_t iumi = 0;
  for(int32_t i=0; tsv_mtxf.read_line() > 0; ++i) {
    if ( i % 1000000 == 0 ) notice("Reading %d lines in %s",i,mtxf.c_str());
    if ( i < 2 ) {
      if ( tsv_mtxf.str_field_at(0)[0] != '%' )
	error("First two lines are expected to start with %%, but observed %s", tsv_mtxf.str_field_at(0));
    }
    else if ( i == 2 ) {
      n_genes = tsv_mtxf.int_field_at(0);
      n_bcds = tsv_mtxf.int_field_at(1);
      num_nonzeros = tsv_mtxf.int_field_at(2);
      sum_umis = num_nonzeros;
      igenes.resize(sum_umis);
      ibcds.resize(sum_umis);
    }
    else {
      int32_t numi = tsv_mtxf.int_field_at(2);
      for(int32_t j=0; j < numi; ++j) {
	igenes[iumi] = tsv_mtxf.int_field_at(0);
	ibcds[iumi]  = tsv_mtxf.int_field_at(1);
	++iumi;

	if ( iumi >= sum_umis ) {
	  sum_umis *= 2;            // resize if the array is full
	  igenes.resize(sum_umis);
	  ibcds.resize(sum_umis);	  
	}
	
      }
      //if ( iumi > sum_umis )
      //error("Exceeded the UMI count %d in the header", sum_umis);
    }
  }
  
  if ( sum_umis > iumi ) {
    sum_umis = iumi; // resize again to fit the size
    igenes.resize(sum_umis);
    ibcds.resize(sum_umis);
  }
  

  // set the random seed
  if ( seed == 0 ) srand(time(NULL));
  else srand((uint32_t) seed);

  notice("Randomizing UMIs with seed %d", seed);

  for(int32_t i=0; i < sum_umis-1; ++i) {
    int32_t r = i + ( rand() % (sum_umis-i) );
    if ( i != r ) { // swap
      int32_t tmp = igenes[i];
      igenes[i] = igenes[r];
      igenes[r] = tmp; 
    }
  }

  notice("Reconstructing digital expression matrix..");
  std::map< int32_t, std::map<int32_t, int32_t> > dge;
  for(int32_t i=0; i < sum_umis; ++i) {
    ++dge[ibcds[i]][igenes[i]];
  }

  notice("Writing new digital expression matrix to %s", outf.c_str());  
  htsFile* wf = hts_open(outf.c_str(), "w");
  hprintf(wf, "%%%%MatrixMarket matrix coordinate integer general\n%%\n%d %d %d\n", n_genes, n_bcds, num_nonzeros);
  for(std::map<int32_t,std::map<int32_t,int32_t> >::iterator it1 = dge.begin(); it1 != dge.end(); ++it1) {
    for(std::map<int32_t,int32_t>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2) {
      hprintf(wf, "%d %d %d\n", it2->first, it1->first, it2->second);
    }
  }
  hts_close(wf);

  notice("Analysis finished");

  return 0;
}
