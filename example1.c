// calculation example using a datafile output by aw_msp_solver.
#include "aw_msp_ivf.h"

// verification function
void verification(AMSP *ms);

int main(int argc, char *argv[])
{
  double complex p,v[3];
  double r[3],f[3];
  int type,sid;
  AMSP ms;
  
  read_dat_amsp(argv[1],&ms); // read datafile 
  print_data_amsp(&ms);       // print data 
  
  r[0]=-0.030; // set x-coordinate
  r[1]= 0.010; // set y-coordinate
  r[2]= 0.002; // set z-coordinate
  type=1;      // select internal field at boundary
  total_pv_amsp(&p,v,r,type,&ms); // calculate total field
  
  printf("sound pressure and particle velocity at r=(%g, %g, %g)\n",r[0],r[1],r[2]);
  printf("p  =% 15.14e %+15.14e I [Pa]\n",creal(p),cimag(p));
  printf("vx =% 15.14e %+15.14e I [m/s]\n",creal(v[0]),cimag(v[0]));
  printf("vy =% 15.14e %+15.14e I [m/s]\n",creal(v[1]),cimag(v[1]));
  printf("vz =% 15.14e %+15.14e I [m/s]\n\n",creal(v[2]),cimag(v[2]));
  
  for(sid=0;sid<ms.n_sphr;sid++){ // sid : sphere id
    F_amsp(sid,f,&ms); // calculate radiation force
    printf("radiation force of sphere id %d\n",sid);
    printf("F=(% 15.14e, % 15.14e, % 15.14e) [N]\n",f[0],f[1],f[2]);
  }
  printf("\n\n");
  
  // verification
  verification(&ms);
  
  free_amsp(&ms);
  return 0;
}

void verification(AMSP *ms)
{
  double complex p1,v1[3],p2,v2[3],p3,v3[3],t1,t2;
  double a,theta,psi,d,r[3],n[3],rt[3],F1[3],N1[3],F2[3],N2[3],h,rc[3];
  int i,ret;
  
  // parameter for calculation point 
  theta=0.3;
  psi=0.5; 
  
  printf("-- verification of relation of sound pressure and particle velocity --\n");
  printf("relation : v = 1/k2 nabla p, k2 = i omega rho.\n");
  h=1.0e-6; 
  printf("calculate the derivative of p using the central difference. h=%g\n",h); 
  for(i=0;i<ms->n_sphr;i++){
    d=ms->sp[i].a*200.0; 
    r[0]=d*sin(theta)*cos(psi)+ms->sp[i].cx;
    r[1]=d*sin(theta)*sin(psi)+ms->sp[i].cy;
    r[2]=d*cos(theta)         +ms->sp[i].cz; // select a point outside sphere
    printf("sphere id %d, r=(%g, %g, %g), ",i,r[0],r[1],r[2]);
    ret=scattered_pv_amsp(&p1,v1,r,ms); 
    if(ret==-1) printf("outside the sphere\n");
    else printf("inside the sphere. please change the calculation point.\n");
    printf("  scattered_pv_amsp() : v_x =% 15.14e %+15.14e I\n",creal(v1[0]),cimag(v1[0]));
    rt[0]=r[0]+h;
    rt[1]=r[1];
    rt[2]=r[2];
    scattered_pv_amsp(&p2,v2,rt,ms);
    rt[0]=r[0]-h;
    scattered_pv_amsp(&p3,v3,rt,ms);
    t1=((p2-p3)/(2.0*h))/ms->aw.k2; // particle velocity using the center difference 
    printf("  central difference  : v_x =% 15.14e %+15.14e I\n",creal(t1),cimag(t1));
    printf("  scattered_pv_amsp() : v_y =% 15.14e %+15.14e I\n",creal(v1[1]),cimag(v1[1]));
    rt[0]=r[0];
    rt[1]=r[1]+h;
    rt[2]=r[2];
    scattered_pv_amsp(&p2,v2,rt,ms);
    rt[1]=r[1]-h;
    scattered_pv_amsp(&p3,v3,rt,ms);
    t1=((p2-p3)/(2.0*h))/ms->aw.k2;
    printf("  central difference  : v_y =% 15.14e %+15.14e I\n",creal(t1),cimag(t1));
    printf("  scattered_pv_amsp() : v_z =% 15.14e %+15.14e I\n",creal(v1[2]),cimag(v1[2]));
    rt[0]=r[0];
    rt[1]=r[1];
    rt[2]=r[2]+h;
    scattered_pv_amsp(&p2,v2,rt,ms);
    rt[2]=r[2]-h;
    scattered_pv_amsp(&p3,v3,rt,ms);
    t1=((p2-p3)/(2.0*h))/ms->aw.k2;
    printf("  central difference  : v_z =% 15.14e %+15.14e I\n",creal(t1),cimag(t1));
  }
 
  
  printf("\n-- verification of boundary conditions --\n");
  printf("internal field p_w, v_w, scattered field p_s, v_s, incident field p_i, v_i\n");
  printf("boundary condition of sound pressure    : p_i + p_s = p_w \n");
  printf("boundary condition of particle velocity : (v_i + v_s) dot n = v_w dot n\n");
  for(i=0;i<ms->n_sphr;i++){
    a=ms->sp[i].a;
    r[0]=a*sin(theta)*cos(psi)+ms->sp[i].cx;
    r[1]=a*sin(theta)*sin(psi)+ms->sp[i].cy;
    r[2]=a*cos(theta)         +ms->sp[i].cz; // select a point on the boundary 
    n[0]=sin(theta)*cos(psi);
    n[1]=sin(theta)*sin(psi);
    n[2]=cos(theta); // normal vector of the bounary
    printf("sphere id %d, r=(%g, %g, %g), n=(%g, %g, %g)\n",i,r[0],r[1],r[2],n[0],n[1],n[2]);
    internal_pv_amsp(&p1,v1,r,ms);
    scattered_pv_amsp(&p2,v2,r,ms);
    incident_pv_amsp(&p3,v3,r,ms);
    printf("  p_i + p_s =% 15.14e %+15.14e I\n",creal(p2+p3),cimag(p2+p3));
    printf("  p_w       =% 15.14e %+15.14e I\n",creal(p1),cimag(p1));
    t1=(v2[0]+v3[0])*n[0]+(v2[1]+v3[1])*n[1]+(v2[2]+v3[2])*n[2];
    t2=v1[0]*n[0]+v1[1]*n[1]+v1[2]*n[2];
    printf("  (v_i + v_s) dot n =% 15.14e %+15.14e I\n",creal(t1),cimag(t1));
    printf("  v_w dot n         =% 15.14e %+15.14e I\n",creal(t2),cimag(t2));
  }
  
  
  printf("\n-- verification of radiation force analysis --\n");
  printf("F_f :  far field method ( summation of field coefficients )\n");
  printf("F_n : near field method ( surface integral of radiation stress tensor )\n");
  rc[0]=0.0;
  rc[1]=0.0;
  rc[2]=0.0; // rotation center 
  for(i=0;i<ms->n_sphr;i++){
    printf("sphere id %d, rotation cetner rc=(%g, %g, %g)\n",i,rc[0],rc[1],rc[2]);
    FN_amsp(i,F1,N1,rc,ms);      // far field method 
    FN_sfint_amsp(i,F2,N2,rc,ms);//
    printf("  F_f=(% 15.14e, % 15.14e, % 15.14e)\n",F1[0],F1[1],F1[2]);
    printf("  F_n=(% 15.14e, % 15.14e, % 15.14e)\n",F2[0],F2[1],F2[2]);
    printf("  N_f=(% 15.14e, % 15.14e, % 15.14e)\n",N1[0],N1[1],N1[2]);    
    printf("  N_n=(% 15.14e, % 15.14e, % 15.14e)\n",N2[0],N2[1],N2[2]);
  }
}
