#include "aw_msp_ivf.h"

int main(int argc, char *argv[])
{
  if(argc!=2){
    printf("This program needs command line argument.\n");
    printf("Usage : %s 'output datafile name'\n",argv[0]);
    exit(0);
  }
  
  AMSP ms;
  
  read_data_amsp(&ms);
  print_data_amsp(&ms);
  setup_amsp(&ms);
  output_node_particles(argv[1],&ms);
  iterative_ops_amsp(&ms); 
  write_dat_amsp(argv[1],&ms);
  free_amsp(&ms);
  
  return 0;
}
