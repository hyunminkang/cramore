#include <stdint.h>
#include <cstdlib>
#include "louvain.h"


int main(int32_t argc, char** argv) {
  louvain l(16);

  srand(argc > 1 ? atoi(argv[1]) : 0);

  /*
  l.add_edge(0, 2, 1.0);
  l.add_edge(0, 3, 1.0);
  l.add_edge(0, 4, 1.0);  
  l.add_edge(0, 5, 1.0);
  l.add_edge(1, 2, 1.0);
  l.add_edge(1, 4, 1.0);
  l.add_edge(1, 7, 1.0);  
  l.add_edge(2, 4, 1.0);
  l.add_edge(2, 5, 1.0);
  l.add_edge(2, 6, 1.0);  
  l.add_edge(3, 7, 1.0);
  l.add_edge(4, 10, 1.0);
  l.add_edge(5, 7, 1.0);
  l.add_edge(5, 11, 1.0);
  l.add_edge(6, 7, 1.0);
  l.add_edge(6, 11, 1.0);
  l.add_edge(8, 9, 1.0);
  l.add_edge(8, 10, 1.0);
  l.add_edge(8, 11, 1.0);
  l.add_edge(8, 14, 1.0);
  l.add_edge(8, 15, 1.0);
  l.add_edge(9, 12, 1.0);
  l.add_edge(9, 14, 1.0);
  l.add_edge(10, 11, 1.0);
  l.add_edge(10, 12, 1.0);
  l.add_edge(10, 13, 1.0);
  l.add_edge(10, 14, 1.0);
  l.add_edge(11, 13, 1.0);
  */

  l.add_edge(0, 1, 1.0);
  l.add_edge(0, 2, 1.0);
  l.add_edge(0, 4, 1.0);
  //l.add_edge(0, 5, 1.0);
  l.add_edge(1, 2, 1.0);
  l.add_edge(1, 4, 1.0);
  l.add_edge(1, 5, 1.0);
  //l.add_edge(2, 4, 1.0);
  l.add_edge(2, 5, 1.0);
  l.add_edge(4, 5, 1.0);
  
  l.add_edge(3, 6, 1.0);
  l.add_edge(3, 7, 1.0);
  l.add_edge(6, 7, 1.0);
  
  l.add_edge(11, 13, 1.0);
  
  l.add_edge(8, 9, 1.0);
  l.add_edge(8, 10, 1.0);
  l.add_edge(8, 12, 1.0);
  l.add_edge(8, 14, 1.0);
  //l.add_edge(8, 15, 1.0);
  l.add_edge(9, 10, 1.0);
  l.add_edge(9, 12, 1.0);
  l.add_edge(9, 14, 1.0);
  l.add_edge(9, 15, 1.0);
  //l.add_edge(10, 12, 1.0);
  l.add_edge(10, 14, 1.0);
  //l.add_edge(10, 15, 1.0);
  //l.add_edge(12, 14, 1.0);
  l.add_edge(12, 15, 1.0);
  //l.add_edge(14, 15, 1.0);

  l.add_edge(0, 3, 1.0);
  l.add_edge(2, 6, 1.0);
  l.add_edge(10, 11, 1.0);  
  
  l.print();

  printf("%s\n\n", l.root->get_name().c_str()); 
  while( l.make_cluster_pass(true) ) {
    printf("%s\n\n", l.root->get_name().c_str());     
  }

  printf("%s\n\n", l.root->get_name().c_str());       
  
  l.print();  

  return 0;
}
