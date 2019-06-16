#include "cramore.h"
#include "dsc_dge.h"
#include <ctime>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>

/////////////////////////////////////////////////////////////////////////
// dge-barcode-summary : Analyze barcodes from digital expression matrix
////////////////////////////////////////////////////////////////////////
int32_t cmdDgeBarcodeSummary(int32_t argc, char** argv) {
  std::string inDir;
  std::string bcdf("barcodes.tsv");
  std::string genef("genes.tsv");
  std::string mtxf("matrix.mtx");
  std::string outPrefix;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("in-dir",&inDir, "Input Directory that contains barcodes, genes, mtx file")
    LONG_STRING_PARAM("barcode",&bcdf, "File containing droplet barcodes")    
    LONG_STRING_PARAM("gene",&genef, "File containing gene ID and symbol")
    LONG_STRING_PARAM("mtx",&mtxf,   "File containing sparse matrices of gene, barcode, UMI")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outPrefix,"Output file prefix")
    LONG_INT_PARAM("verbose", &globalVerbosityThreshold, "Turn on verbose mode with specific verbosity threshold. 0: fully verbose, 100 : no verbose messages")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( inDir.empty() || outPrefix.empty() )
    error("Missing required option(s) : --in-dir, --out");

  notice("Analysis Started");

  std::vector<std::string> bcds;
  std::map<std::string,int> bcd2idx;
  std::vector<std::string> geneIDs;
  std::vector<std::string> geneSymbols;

  if ( inDir[inDir.size()-1] != '/' ) inDir += '/';
  char fbcd[65535];  
  sprintf(fbcd, "%s%s", inDir.c_str(), bcdf.c_str());
  char fgene[65535];  
  sprintf(fgene, "%s%s", inDir.c_str(), genef.c_str());
  char fmtx[65535];  
  sprintf(fmtx, "%s%s", inDir.c_str(), mtxf.c_str());

  dsc_dge(fbcd, fgene, fmtx);

  return 0;
}
