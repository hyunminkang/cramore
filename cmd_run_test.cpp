#include "cramore.h"
#include "tsv_reader.h"
#include <sys/stat.h>
#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>

//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0

void process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}

int32_t cmdRunTest(int32_t argc, char** argv) {
  std::string inputFile;
  
  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input files", NULL)
    LONG_STRING_PARAM("in",&inputFile, "Input file name")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( inputFile.empty() )
    error("Missing required parameter --in");

  double vm, rss;
  process_mem_usage(vm, rss);
  std::cout << "VM: " << vm << "; RSS: " << rss << std::endl;
  
  
  tsv_reader tr(inputFile.c_str());
  for(int32_t i=0;tr.read_line() > 0; ++i) {
    if ( (i+1) % 1000000 == 0 ) 
      notice("Processing %d lines from %s", i+1, inputFile.c_str());
  }
  notice("Finished Processing %d lines from %s", tr.nlines, inputFile.c_str());
  
  process_mem_usage(vm, rss);
  std::cout << "VM: " << vm << "; RSS: " << rss << std::endl;

  notice("Analysis finished");

  return 0;
}
