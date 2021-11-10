#include "aw_msp_ivf.h"

int main()
{
  AMSP ms;
  
  read_data_amsp(&ms);
  print_data_amsp(&ms);
  if(ms.n_sphr!=1) {
    printf("This program is available for single sphere. Exit\n");
    exit(1);
  }
  if(ms.aw.n_pw!=0 || ms.aw.n_bb!=1 || ms.aw.n_fb!=0){
    printf("This program is available for single Bessel beam. Exit\n");
    exit(1);
  }
  free_amsp(&ms);

  FILE *fp;
  double a,a_min,a_max,da,F[3],f_min,f_max,df,f;
  int ni,nj,i,j;
  
  // radius setting
  a_min= 1e-3;     // minimum radius [m]
  a_max=20.e-3;    // maximum radius [m]
  ni=150;          // sampling number
  da=(a_max-a_min)/(ni-1); // delta a
  
  // frequency setting
  f_min= 10e3;     // minimum frequency [Hz]
  f_max=200e3;     // maximum frequency [Hz]
  nj=150;          // sampling number 
  df=(f_max-f_min)/(nj-1); // delta f
  
  if((fp=fopen("radiation_force2.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp,"%s\n","# frequency radius F_z");
  
  for(j=0;j<nj;j++){
    f=f_min+(double)j*df;
    printf("\r%d/%d",j+1,nj);  fflush(stdout);
    for(i=0;i<ni;i++){
      read_data_amsp(&ms);
      a=a_min+(double)i*da;
      ms.sp[0].a=a;       // reset radius
      ms.aw.bd.bb[0].f=f; // reset frequency
      setup_amsp(&ms); // setup coefficients
      F_amsp(0,F,&ms); // calculate radiation force
      fprintf(fp,"%g %g %15.14e\n",f,a,F[2]); 
      free_amsp(&ms); // free allocated memory
    }
    fprintf(fp,"\n");
  }
  
  printf("\ndone\n");
  return 0;
}
