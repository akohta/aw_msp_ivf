#include "aw_msp_ivf.h"


int incident_pv_amsp(double complex *p,double complex *v,double *r,AMSP *ms)
{
  calc_maw_pv(p,v,r,&(ms->aw));
  
  return -1;
}

int internal_pv_amsp(double complex *p,double complex *v,double *r,AMSP *ms)
{
  int internal_pv_ss(double complex *p,double complex *v,double *xb,ASP *sp,Maw *aw);
  
  int id,ret;

  id=sphere_id_type1(r,ms);
  if(id>=0){
    ret=internal_pv_ss(p,v,r,&(ms->sp[id]),&(ms->aw));
    if(ret<0) {
      printf("error in internal_pv_amsp(), internal_pv_ss(). Exit...\n");
      exit(1);
    }
  }
  else { // outside the spheres
    *p=0.0;
    v[0]=0.0;  v[1]=0.0;  v[2]=0.0;
    ret=0;
  }
  return id; 
}
 
int scattered_pv_amsp(double complex *p,double complex *v,double *r,AMSP *ms)
{
  int scattered_pv_ss(double complex *p,double complex *v,double *xb,ASP *sp,Maw *aw);
  
  double complex tp,tv[3];
  int id,ret,i;
  
  id=sphere_id_type2(r,ms);
  if(id<0){ 
    *p=0.0;
    v[0]=0.0;  v[1]=0.0;  v[2]=0.0;
    for(i=0;i<ms->n_sphr;i++){
      ret=scattered_pv_ss(&tp,tv,r,&(ms->sp[i]),&(ms->aw));
      if(ret<0) {
        printf("error in scattered_pv_amsp(), scattered_pv_ss(). Exit...\n");
        exit(1);
      }
      *p+=tp;
      v[0]+=tv[0];  v[1]+=tv[1];  v[2]+=tv[2];
    }
  }
  else { // inside the spheres
   *p=0.0;
    v[0]=0.0;  v[1]=0.0;  v[2]=0.0;
    ret=0;
  }
  return id; 
}

int total_pv_amsp(double complex *p,double complex *v,double *r,int type,AMSP *ms)
{
  double complex pi,vi[3];
  int id;
  
  if(type==1) id=sphere_id_type1(r,ms);
  else id=sphere_id_type2(r,ms);
  
  if(id>=0) internal_pv_amsp(p,v,r,ms); 
  else { 
    scattered_pv_amsp(p,v,r,ms);
    incident_pv_amsp(&pi,vi,r,ms);
    *p+=pi;
    v[0]+=vi[0];
    v[1]+=vi[1];
    v[2]+=vi[2];
  }
  return id;
}

/////////////////////////////////////////////////////////////////////////////
int scattered_pv_ss(double complex *p,double complex *v,double *xb,ASP *sp,Maw *aw)
{
  double complex ph,dphr,dpht,dphp,Yp,Ym,dYp,dYm,dep,expi,aslmp,aslmm;
  double r,i_r,rxy,i_rxy,r2,cos_t,sin_t,cos_p,sin_p,ker,i_ker,ke,sig,x,y,z,cf;
  int i,l,m,tt,lm,ai;
  size_t ms,lmax=(size_t)sp->l_limit;
  ms=gsl_sf_legendre_array_n(lmax);
  lm=(int)lmax;
  
  double *sphPlm =(double *)m_alloc2(ms,sizeof(double),"scattered_pv_ss(),*sphPlm");
  double *dsphPlm=(double *)m_alloc2(ms,sizeof(double),"scattered_pv_ss(),*dsphPlm");
  double  *jl=(double *)m_alloc2(lm+2,sizeof(double),"scattered_pv_ss(),*jl");
  double  *yl=(double *)m_alloc2(lm+2,sizeof(double),"scattered_pv_ss(),*yl");
  double complex  *hl=(double complex *)m_alloc2(lm+2,sizeof(double complex),"scattered_pv_ss(),hl");
  double complex *dhl=(double complex *)m_alloc2(lm+1,sizeof(double complex),"scattered_sv_ss(),dhl");

  x=xb[0]-sp->cx;  y=xb[1]-sp->cy;  z=xb[2]-sp->cz;
  r2=x*x+y*y+z*z;    r=sqrt(r2);
  if(r<sp->a-MEPS) return -1; // inside the sphere 
  i_r=1.0/r;
  rxy=sqrt(x*x+y*y);  
  if(rxy==0.0){ // x==0,y==0
    x=z*0.7e-7;
    y=z*0.7e-7;
    r2=x*x+y*y+z*z;    r=sqrt(r2);  i_r=1.0/r;
    rxy=sqrt(x*x+y*y); 
  }
  i_rxy=1.0/rxy;
  cos_t=z*i_r;    sin_t=rxy*i_r;
  cos_p=x*i_rxy;  sin_p=y*i_rxy;

  ke =aw->k0;
  ker=ke*r;
  i_ker=1.0/ker;
  gsl_sf_bessel_jl_array(lm+1,ker,jl);
  gsl_sf_bessel_yl_array(lm+1,ker,yl);
  for(i=0;i<=lm+1;i++) hl[i]=jl[i]+I*yl[i];
  for(i=0;i<=lm;i++){
    cf=(double)i*i_ker;
    dhl[i]=cf*hl[i]-hl[i+1];
  }
  gsl_sf_legendre_deriv_alt_array_e(GSL_SF_LEGENDRE_SPHARM,lm,cos_t,-1,sphPlm,dsphPlm);
  dep=cos_p+I*sin_p; expi=1.0;
  sig=1.0;
  
  ph=0.0;
  dphr=0.0;
  dpht=0.0;
  dphp=0.0;
  tt=0;  m=0;
  for(l=0;l<=lm;l++){
    ai=gsl_sf_legendre_array_index(l,m);
    Yp = sphPlm[ai];
    dYp=dsphPlm[ai];
    aslmp=sp->dd.cs[l]*sp->dd.ailm[tt];
    ph  +=aslmp* hl[l]*Yp;
    dphr+=aslmp*dhl[l]*Yp;
    dpht+=aslmp* hl[l]*dYp;
    dphp+=aslmp* hl[l]*Yp*(double)m;
    tt++;
  }
  for(m=1;m<=lm;m++){
    expi*=dep;
    sig*=-1.0;
    for(l=m;l<=lm;l++){
      ai=gsl_sf_legendre_array_index(l,m);
      Yp = sphPlm[ai]*expi;
      dYp=dsphPlm[ai]*expi;
      Ym =sig*conj( Yp);
      dYm=sig*conj(dYp);
      aslmp=sp->dd.cs[l]*sp->dd.ailm[tt+0];
      aslmm=sp->dd.cs[l]*sp->dd.ailm[tt+1];
      
      ph  += hl[l]*(aslmp* Yp+aslmm* Ym);
      dphr+=dhl[l]*(aslmp* Yp+aslmm* Ym);
      dpht+= hl[l]*(aslmp*dYp+aslmm*dYm);
      dphp+= hl[l]*(aslmp* Yp-aslmm* Ym)*(double)m;
      tt+=2;
    }
  }
  dphr*=ke;
  dphp*=I;
  
  *p=-aw->k2*ph;
  v[0]=-(dphr*sin_t*cos_p+dpht*cos_t*cos_p*i_r-dphp*sin_p*i_rxy);
  v[1]=-(dphr*sin_t*sin_p+dpht*cos_t*sin_p*i_r+dphp*cos_p*i_rxy);
  v[2]=-(dphr*cos_t-dpht*sin_t*i_r);
  
  free(sphPlm);  free(dsphPlm);
  free(jl); free(yl);
  free(hl); free(dhl);
  
  return 0;
}

int scattered_phi_ss(double complex *phi,double *xb,ASP *sp,Maw *aw)
{
  double complex ph,Yp,Ym,dep,expi,aslmp,aslmm;
  double r,i_r,rxy,i_rxy,r2,cos_t,cos_p,sin_p,ker,ke,sig,x,y,z;
  int i,l,m,tt,lm,ai;
  size_t ms,lmax=(size_t)sp->l_limit;
  ms=gsl_sf_legendre_array_n(lmax);
  lm=(int)lmax;
  
  double *sphPlm =(double *)m_alloc2(ms,sizeof(double),"scattered_pv_ss(),*sphPlm");
  double  *jl=(double *)m_alloc2(lm+2,sizeof(double),"scattered_pv_ss(),*jl");
  double  *yl=(double *)m_alloc2(lm+2,sizeof(double),"scattered_pv_ss(),*yl");
  double complex  *hl=(double complex *)m_alloc2(lm+2,sizeof(double complex),"scattered_pv_ss(),hl");

  x=xb[0]-sp->cx;  y=xb[1]-sp->cy;  z=xb[2]-sp->cz;
  r2=x*x+y*y+z*z;    r=sqrt(r2);
  if(r<sp->a-MEPS) return -1; // inside the sphere 
  i_r=1.0/r;
  rxy=sqrt(x*x+y*y);  
  if(rxy==0.0){ // x==0,y==0
    x=z*0.7e-7;
    y=z*0.7e-7;
    r2=x*x+y*y+z*z;    r=sqrt(r2);  i_r=1.0/r;
    rxy=sqrt(x*x+y*y); 
  }
  i_rxy=1.0/rxy;
  cos_t=z*i_r;  
  cos_p=x*i_rxy;  sin_p=y*i_rxy;

  ke =aw->k0;
  ker=ke*r;
  gsl_sf_bessel_jl_array(lm+1,ker,jl);
  gsl_sf_bessel_yl_array(lm+1,ker,yl);
  for(i=0;i<=lm+1;i++) hl[i]=jl[i]+I*yl[i];

  gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM,lm,cos_t,-1,sphPlm);
  dep=cos_p+I*sin_p; expi=1.0;
  sig=1.0;
  
  ph=0.0;
  tt=0;  m=0;
  for(l=0;l<=lm;l++){
    ai=gsl_sf_legendre_array_index(l,m);
    Yp = sphPlm[ai];
    aslmp=sp->dd.cs[l]*sp->dd.ailm[tt];
    ph  +=aslmp* hl[l]*Yp;
    tt++;
  }
  for(m=1;m<=lm;m++){
    expi*=dep;
    sig*=-1.0;
    for(l=m;l<=lm;l++){
      ai=gsl_sf_legendre_array_index(l,m);
      Yp = sphPlm[ai]*expi;
      Ym =sig*conj( Yp);
      aslmp=sp->dd.cs[l]*sp->dd.ailm[tt+0];
      aslmm=sp->dd.cs[l]*sp->dd.ailm[tt+1];
      
      ph  += hl[l]*(aslmp* Yp+aslmm* Ym);
      tt+=2;
    }
  }
  
  *phi=ph;
  
  free(sphPlm);
  free(jl); free(yl);
  free(hl);
  
  return 0; 
  
}

int internal_pv_ss(double complex *p,double complex *v,double *xb,ASP *sp,Maw *aw)
{
  double complex ph,dphr,dpht,dphp,Yp,Ym,dYp,dYm,dep,expi,awlmp,awlmm;
  double r,i_r,rxy,i_rxy,r2,cos_t,sin_t,cos_p,sin_p,ki,kir,i_kir,sig,x,y,z,cf;
  int i,l,m,tt,lm,ai;
  size_t ms,lmax=(size_t)sp->l_limit;
  ms=gsl_sf_legendre_array_n(lmax);
  lm=(int)lmax;
  
  double *sphPlm =(double *)m_alloc2(ms,sizeof(double),"internal_pv_ss(),*sphPlm");
  double *dsphPlm=(double *)m_alloc2(ms,sizeof(double),"internal_pv_ss(),*dsphPlm");
  double  *jl=(double *)m_alloc2(lm+2,sizeof(double),"internal_pv_ss(),*jl");
  double *djl=(double *)m_alloc2(lm+1,sizeof(double),"internal_pv_ss(),*djl");

  x=xb[0]-sp->cx;  y=xb[1]-sp->cy;  z=xb[2]-sp->cz;
  r2=x*x+y*y+z*z;    r=sqrt(r2);
  if(r>sp->a+MEPS) return -1; // outside the sphere
  rxy=sqrt(x*x+y*y); 
  if(rxy==0.0){ // x==0,y==0
    if(z==0.0){
      z=( (sp->a <= 2.0*M_PI/sp->dd.k0)? sp->a : 2.0*M_PI/sp->dd.k0 )*sqrt(MEPS);
      x=z;
      y=z;
    }
    else {
      x=z*0.7e-7;
      y=z*0.7e-7;
    }
    r2=x*x+y*y+z*z;    r=sqrt(r2);
    rxy=sqrt(x*x+y*y); 
  }
  i_r=1.0/r;
  i_rxy=1.0/rxy;
  cos_t=z*i_r;    sin_t=rxy*i_r;
  cos_p=x*i_rxy;  sin_p=y*i_rxy;

  ki =sp->dd.k0;
  kir=ki*r;
  i_kir=1.0/kir;
  gsl_sf_bessel_jl_array(lm+1,kir,jl);
  for(i=0;i<=lm;i++){
    cf=(double)i*i_kir;
    djl[i]=cf*jl[i]-jl[i+1];
  }
  gsl_sf_legendre_deriv_alt_array_e(GSL_SF_LEGENDRE_SPHARM,lm,cos_t,-1,sphPlm,dsphPlm);
  dep=cos_p+I*sin_p; expi=1.0;
  sig=1.0;

  ph=0.0;
  dphr=0.0;
  dpht=0.0;
  dphp=0.0;
  tt=0;  m=0;
  for(l=0;l<=lm;l++){
    ai=gsl_sf_legendre_array_index(l,m);
    Yp = sphPlm[ai];
    dYp=dsphPlm[ai];
    awlmp=sp->dd.cw[l]*sp->dd.ailm[tt];

    ph  +=awlmp* jl[l]* Yp;
    dphr+=awlmp*djl[l]* Yp;
    dpht+=awlmp* jl[l]*dYp;
    dphp+=awlmp* jl[l]* Yp*(double)m;
    tt++;
  }
  for(m=1;m<=lm;m++){
    expi*=dep;
    sig*=-1.0;
    for(l=m;l<=lm;l++){
      ai=gsl_sf_legendre_array_index(l,m);
      Yp = sphPlm[ai]*expi;
      dYp=dsphPlm[ai]*expi;
      Ym =sig*conj( Yp);
      dYm=sig*conj(dYp);
      awlmp=sp->dd.cw[l]*sp->dd.ailm[tt+0];
      awlmm=sp->dd.cw[l]*sp->dd.ailm[tt+1];
      
      ph  += jl[l]*(awlmp* Yp+awlmm* Ym);
      dphr+=djl[l]*(awlmp* Yp+awlmm* Ym);
      dpht+= jl[l]*(awlmp*dYp+awlmm*dYm);
      dphp+= jl[l]*(awlmp* Yp-awlmm* Ym)*(double)m;
      tt+=2;
    }
  }
  dphr*=ki;
  dphp*=I;

  *p=-sp->dd.k2*ph;
  v[0]=-(dphr*sin_t*cos_p+dpht*cos_t*cos_p*i_r-dphp*sin_p*i_rxy);
  v[1]=-(dphr*sin_t*sin_p+dpht*cos_t*sin_p*i_r+dphp*cos_p*i_rxy);
  v[2]=-(dphr*cos_t-dpht*sin_t*i_r);

  free(sphPlm);  free(dsphPlm);
  free(jl);      free(djl);

  return 0;
}
