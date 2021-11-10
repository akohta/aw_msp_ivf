#include "aw_msp_ivf.h"

int main()
{
  AMSP ms;
  FILE *ff;
  double rang,dr,r[3],vf[3];
  int max,i,j;
  
  read_data_amsp(&ms);
  print_data_amsp(&ms);
  if(ms.n_sphr!=1) {
    printf("This program is available for single sphere only. Exit\n");
    exit(1);
  }
  setup_amsp(&ms);

  max=15;
  rang=2.0*ms.aw.lambda0;
  dr=rang*2/(double)(max-1);
  free_amsp(&ms);

  // x=0 plane 
  if((ff=fopen("force_yz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(ff,"%s\n","# y z Fx[N] Fy[N] Fz[N]");

  r[0]=0.0;
  for(i=0;i<max;i++){
    r[1]=-rang+(double)i*dr;
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      read_data_amsp(&ms);
      // reset shpere position
      ms.sp[0].cx=r[0];      ms.sp[0].cy=r[1];      ms.sp[0].cz=r[2];
      setup_amsp(&ms);
      F_amsp(0,vf,&ms);
      fprintf(ff,"%g %g %g %g %g\n",r[1],r[2],vf[0],vf[1],vf[2]);
      free_amsp(&ms);
     }
    fprintf(ff,"\n");
  }
  fclose(ff);
  printf("x=0 plane is finished\n"); fflush(stdout);
  

  // y=0 plane 
  if((ff=fopen("force_xz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(ff,"%s\n","# x z Fx[N] Fy[N] Fz[N]");

  r[1]=0.0;
  for(i=0;i<max;i++){
    r[0]=-rang+(double)i*dr;
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      read_data_amsp(&ms);
      // reset shpere position
      ms.sp[0].cx=r[0];      ms.sp[0].cy=r[1];      ms.sp[0].cz=r[2];
      setup_amsp(&ms);
      F_amsp(0,vf,&ms);
      fprintf(ff,"%g %g %g %g %g\n",r[0],r[2],vf[0],vf[1],vf[2]);
      free_amsp(&ms);
     }
    fprintf(ff,"\n");
  }
  fclose(ff);
  printf("y=0 plane is finished\n"); fflush(stdout);

  
  // z=0 plane
  if((ff=fopen("force_xy.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(ff,"%s\n","# x y Fx[N] Fy[N] Fz[N]");

  r[2]=0.0;
  for(i=0;i<max;i++){
    r[0]=-rang+(double)i*dr;
    for(j=0;j<max;j++){
      r[1]=-rang+(double)j*dr;
      read_data_amsp(&ms);
      // reset shpere position
      ms.sp[0].cx=r[0];      ms.sp[0].cy=r[1];      ms.sp[0].cz=r[2];
      setup_amsp(&ms);
      F_amsp(0,vf,&ms);
      fprintf(ff,"%g %g %g %g %g\n",r[0],r[1],vf[0],vf[1],vf[2]);
      free_amsp(&ms);
    }
    fprintf(ff,"\n");
  }
  fclose(ff);
  printf("z=0 plane is finished.\n");
  
  return 0;
}
