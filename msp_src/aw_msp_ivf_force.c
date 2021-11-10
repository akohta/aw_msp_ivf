#include "aw_msp_ivf.h"

void F_amsp(int id,double *f,AMSP *ms)
{
  void force_sphr(double *f,ASP *sp,Maw *aw,int spid);
  
  force_sphr(f,&(ms->sp[id]),&(ms->aw),id);
}

void FN_amsp(int id,double *f,double *n,double *rc,AMSP *ms)
{
  void force_sphr(double *f,ASP *sp,Maw *aw,int spid);
  
  force_sphr(f,&(ms->sp[id]),&(ms->aw),id);
  
  n[0]=f[1]*(rc[2]-ms->sp[id].cz)-f[2]*(rc[1]-ms->sp[id].cy);
  n[1]=f[2]*(rc[0]-ms->sp[id].cx)-f[0]*(rc[2]-ms->sp[id].cz);
  n[2]=f[0]*(rc[1]-ms->sp[id].cy)-f[1]*(rc[0]-ms->sp[id].cx);
}

void FN_sfint_amsp(int id,double *F,double *N,double *rc,AMSP *ms)
{
  double complex p,v[3];
  double ct,st,cp,sp,a,r[3],n[3],c1,c2,c3,t1,fx,fy,fz,nx,ny,nz,T[9];
  int i,j,ret,l,m;
  
  F[0]=0.0;  F[1]=0.0;  F[2]=0.0;
  N[0]=0.0;  N[1]=0.0;  N[2]=0.0;
  a=ms->sp[id].a+MEPS;
  c1=0.25/(ms->aw.c0*ms->aw.c0*ms->aw.rho0); // 1/(4*K0)
  c2=0.25*ms->aw.rho0; // rho0/4
  c3=0.50*ms->aw.rho0; // rho0/2
  
  for(i=0;i<ms->sp[id].dd.nt;i++){ // theta
    ct=ms->sp[id].dd.ct[i];
    st=ms->sp[id].dd.st[i];
    fx=0.0;    fy=0.0;    fz=0.0;
    nx=0.0;    ny=0.0;    nz=0.0;
    #pragma omp parallel for schedule(dynamic) reduction(+:fx,fy,fz,nx,ny,nz) private(cp,sp,n,r,ret,p,v,t1,l,m,T) // parallel
    for(j=0;j<ms->sp[id].dd.np;j++){ // psi
      cp=ms->sp[id].dd.cp[j];
      sp=ms->sp[id].dd.sp[j];
      
      n[0]=st*cp;
      n[1]=st*sp;
      n[2]=ct; // normal vector 
      r[0]=a*st*cp+ms->sp[id].cx;
      r[1]=a*st*sp+ms->sp[id].cy;
      r[2]=a*ct   +ms->sp[id].cz;
      ret=total_pv_amsp(&p,v,r,2,ms); 
      if(ret>=0){ // internal field returned
        printf("error in FN_sfint_amsp(). total_pv_amsp(). Exit...\n");
        exit(1);
      }
      // radiation stress tensor
      t1=c1*creal(p*conj(p))-c2*creal(v[0]*conj(v[0])+v[1]*conj(v[1])+v[2]*conj(v[2]));
      for(l=0;l<3;l++) for(m=0;m<3;m++) T[l*3+m]=c3*v[l]*conj(v[m]);
      T[0*3+0]+=t1;
      T[1*3+1]+=t1;
      T[2*3+2]+=t1; 
      // force
      fx+=(T[0*3+0]*n[0]+T[0*3+1]*n[1]+T[0*3+2]*n[2])*ms->sp[id].dd.wp[j];
      fy+=(T[1*3+0]*n[0]+T[1*3+1]*n[1]+T[1*3+2]*n[2])*ms->sp[id].dd.wp[j];
      fz+=(T[2*3+0]*n[0]+T[2*3+1]*n[1]+T[2*3+2]*n[2])*ms->sp[id].dd.wp[j]; 
      // torque
      nx+=( ((r[1]-rc[1])*T[2*3+0]-(r[2]-rc[2])*T[1*3+0])*n[0]
           +((r[1]-rc[1])*T[2*3+1]-(r[2]-rc[2])*T[1*3+1])*n[1]
           +((r[1]-rc[1])*T[2*3+2]-(r[2]-rc[2])*T[1*3+2])*n[2])*ms->sp[id].dd.wp[j];
      ny+=( ((r[2]-rc[2])*T[0*3+0]-(r[0]-rc[0])*T[2*3+0])*n[0]
           +((r[2]-rc[2])*T[0*3+1]-(r[0]-rc[0])*T[2*3+1])*n[1]
           +((r[2]-rc[2])*T[0*3+2]-(r[0]-rc[0])*T[2*3+2])*n[2])*ms->sp[id].dd.wp[j];
      nz+=( ((r[0]-rc[0])*T[1*3+0]-(r[1]-rc[1])*T[0*3+0])*n[0]
           +((r[0]-rc[0])*T[1*3+1]-(r[1]-rc[1])*T[0*3+1])*n[1]
           +((r[0]-rc[0])*T[1*3+2]-(r[1]-rc[1])*T[0*3+2])*n[2])*ms->sp[id].dd.wp[j];
    }
    F[0]+=fx*st*ms->sp[id].dd.wt[i];
    F[1]+=fy*st*ms->sp[id].dd.wt[i];
    F[2]+=fz*st*ms->sp[id].dd.wt[i];
    N[0]+=nx*st*ms->sp[id].dd.wt[i];
    N[1]+=ny*st*ms->sp[id].dd.wt[i];
    N[2]+=nz*st*ms->sp[id].dd.wt[i];
  }
  
  F[0]*=-a*a;
  F[1]*=-a*a;
  F[2]*=-a*a;
  N[0]*=-a*a;
  N[1]*=-a*a;
  N[2]*=-a*a;
}

///////////////////////////////////////////
void force_sphr(double *f,ASP *sp,Maw *aw,int spid)
{
  double complex pickup_ailm(int l,int m,ASP *sp); // aw_msp_ivf.c
  double coef_A_lm(int l,int m);
  double coef_B_lm(int l,int m);

  double complex fxy,fz,ail0m0,ailpm0,ailmm0,ailpmp,ailmmp,asl0m0,aslpm0,aslmm0,aslpmp,aslmmp,c1,c2,c3,c4;
  double A1,A2,B1,B2,rho0e,ifz,af;
  int l,m,ll,sm;
  static int st=0;
  
  ll=sp->l_limit;
  rho0e=aw->rho0;

  fxy=0.0;
  fz=0.0;
  
  m=0;
  for(l=0;l<ll;l++){
    sm=m;
    ail0m0=pickup_ailm(l+0,sm+0,sp);  asl0m0=sp->dd.cs[l+0]*ail0m0;
    ailpm0=pickup_ailm(l+1,sm+0,sp);  aslpm0=sp->dd.cs[l+1]*ailpm0;
    ailmm0=pickup_ailm(l-1,sm+0,sp);  aslmm0=sp->dd.cs[l-1]*ailmm0;
    ailpmp=pickup_ailm(l+1,sm+1,sp);  aslpmp=sp->dd.cs[l+1]*ailpmp;
    ailmmp=pickup_ailm(l-1,sm+1,sp);  aslmmp=sp->dd.cs[l-1]*ailmmp;
    c1= (ail0m0*conj(aslpm0)+asl0m0*conj(ailpm0)+2.0*asl0m0*conj(aslpm0));
    c2=-(ail0m0*conj(aslmm0)+asl0m0*conj(ailmm0)+2.0*asl0m0*conj(aslmm0));
    c3= (ail0m0*conj(aslpmp)+asl0m0*conj(ailpmp)+2.0*asl0m0*conj(aslpmp));
    c4=-(ail0m0*conj(aslmmp)+asl0m0*conj(ailmmp)+2.0*asl0m0*conj(aslmmp));
    A1=coef_A_lm(l+1,sm);    A2=coef_A_lm(l+0,sm);
    B1=coef_B_lm(l+1,sm);    B2=coef_B_lm(l+0,-sm-1);
    fz +=A1*c1+A2*c2;
    fxy+=B1*c3-B2*c4;
  }

  for(m=1;m<ll;m++){
    for(l=m;l<ll;l++){
      // m>0
      sm=m;
      ail0m0=pickup_ailm(l+0,sm+0,sp);  asl0m0=sp->dd.cs[l+0]*ail0m0;
      ailpm0=pickup_ailm(l+1,sm+0,sp);  aslpm0=sp->dd.cs[l+1]*ailpm0;
      ailmm0=pickup_ailm(l-1,sm+0,sp);  aslmm0=sp->dd.cs[l-1]*ailmm0;
      ailpmp=pickup_ailm(l+1,sm+1,sp);  aslpmp=sp->dd.cs[l+1]*ailpmp;
      ailmmp=pickup_ailm(l-1,sm+1,sp);  aslmmp=sp->dd.cs[l-1]*ailmmp;
      c1= (ail0m0*conj(aslpm0)+asl0m0*conj(ailpm0)+2.0*asl0m0*conj(aslpm0));
      c2=-(ail0m0*conj(aslmm0)+asl0m0*conj(ailmm0)+2.0*asl0m0*conj(aslmm0));
      c3= (ail0m0*conj(aslpmp)+asl0m0*conj(ailpmp)+2.0*asl0m0*conj(aslpmp));
      c4=-(ail0m0*conj(aslmmp)+asl0m0*conj(ailmmp)+2.0*asl0m0*conj(aslmmp));
      A1=coef_A_lm(l+1,sm);      A2=coef_A_lm(l+0,sm);
      B1=coef_B_lm(l+1,sm);      B2=coef_B_lm(l+0,-sm-1);
      fz +=A1*c1+A2*c2;
      fxy+=B1*c3-B2*c4;
      // m<0
      sm=-m;
      ail0m0=pickup_ailm(l+0,sm+0,sp);  asl0m0=sp->dd.cs[l+0]*ail0m0;
      ailpm0=pickup_ailm(l+1,sm+0,sp);  aslpm0=sp->dd.cs[l+1]*ailpm0;
      ailmm0=pickup_ailm(l-1,sm+0,sp);  aslmm0=sp->dd.cs[l-1]*ailmm0;
      ailpmp=pickup_ailm(l+1,sm+1,sp);  aslpmp=sp->dd.cs[l+1]*ailpmp;
      ailmmp=pickup_ailm(l-1,sm+1,sp);  aslmmp=sp->dd.cs[l-1]*ailmmp;
      c1= (ail0m0*conj(aslpm0)+asl0m0*conj(ailpm0)+2.0*asl0m0*conj(aslpm0));
      c2=-(ail0m0*conj(aslmm0)+asl0m0*conj(ailmm0)+2.0*asl0m0*conj(aslmm0));
      c3= (ail0m0*conj(aslpmp)+asl0m0*conj(ailpmp)+2.0*asl0m0*conj(aslpmp));
      c4=-(ail0m0*conj(aslmmp)+asl0m0*conj(ailmmp)+2.0*asl0m0*conj(aslmmp));
      A1=coef_A_lm(l+1,sm);      A2=coef_A_lm(l+0,sm);
      B1=coef_B_lm(l+1,sm);      B2=coef_B_lm(l+0,-sm-1);
      fz +=A1*c1+A2*c2;
      fxy+=B1*c3-B2*c4;
    }
  }

  fz *=-0.25*I*rho0e;
  fxy*= 0.25*I*rho0e;
  
  f[0]=creal(fxy);
  f[1]=cimag(fxy);
  f[2]=creal(fz);
  
  // error check and correction
  ifz=fabs(cimag(fz)); // theoretically 0
  if(ifz>fabs(f[2])*sqrt(MEPS) && st==0){
    st++;
    printf("sphere id %d, force_sphr()\n",spid);
    printf("the imaginary part of fz must be 0 theoretically, but it has a value that cannot be ignored.\n");
    printf("it is recommended to increase the limit of l and calculate again.\n");
    printf("fz=% 15.14e %+15.14e I\n",creal(fz),cimag(fz));
  }
  af=vabs_d(f); // defined in my_utlis.h
  if(fabs(f[0])<af*MEPS) f[0]=0.0;
  if(fabs(f[1])<af*MEPS) f[1]=0.0;
  if(fabs(f[2])<af*MEPS) f[2]=0.0;
}

double coef_A_lm(int l,int m)
{
  return sqrt((double)((l+m)*(l-m))/(double)((2*l+1)*(2*l-1)));
}

double coef_B_lm(int l,int m)
{
  return sqrt((double)((l+m)*(l+m+1))/(double)((2*l+1)*(2*l-1)));
}
