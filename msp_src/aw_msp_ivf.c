#include "aw_msp_ivf.h"

void read_data_amsp(AMSP *ms)
{
  FILE *fp;
  char buf[256]="";  int ti;  double td;
  int num,nc;
  
  if((fp=fopen(fn_sphr,"rt"))==NULL){    printf("Can not open the '%s' file. Exit...\n",fn_sphr);    exit(1);  }
  fgets(buf,256,fp);  fgets(buf,256,fp);

  fscanf(fp,"%d\n",&ti);   num=ti;
  if(num==0){    printf("No sphere defined. Exit...\n"); exit(1);  }
  fgets(buf,256,fp);

  ms->n_sphr=num;
  ms->sp=(ASP *)m_alloc2(num,sizeof(ASP),"read_data_amsp(), sm->sp");
  
  for(nc=0;nc<num;nc++){
    fscanf(fp,"%lf",&td);  ms->sp[nc].rho0=td;
    fscanf(fp,"%lf",&td);  ms->sp[nc].c0=td;
    fscanf(fp,"%lf",&td);  ms->sp[nc].a =td;
    fscanf(fp,"%lf",&td);  ms->sp[nc].cx=td;
    fscanf(fp,"%lf",&td);  ms->sp[nc].cy=td;
    fscanf(fp,"%lf",&td);  ms->sp[nc].cz=td;
    fscanf(fp,"%d",&ti);   ms->sp[nc].bsn=ti;
    fscanf(fp,"%d",&ti);   ms->sp[nc].bdv=ti;
    fscanf(fp,"%d\n",&ti);   ms->sp[nc].l_limit=ti;
  }
  
  fclose(fp);
  
  // multi_aw
  init_maw(&(ms->aw));      // initialize
  read_data_maw(&(ms->aw)); // search and read beam datafile
}

void print_data_amsp(AMSP *ms)
{
  int nc;
  
  print_data_maw(&(ms->aw));     // print beam data
  
  // print sphere data 
  printf("---- sphere data ( %s ) ----\n",fn_sphr);
  printf("number of spheres                             : %16d\n",ms->n_sphr);
  for(nc=0;nc<ms->n_sphr;nc++){
    printf("  Sphere ID %d\n",nc);
    printf("medium density                       [kg/m^3] : %16.15g\n",ms->sp[nc].rho0);
    printf("speed of sound in the medium            [m/s] : %16.15g\n",ms->sp[nc].c0);
    printf("radius of sphere                              : %16.15g\n",ms->sp[nc].a);
    printf("x-coordinate of sphere center                 : %16.15g\n",ms->sp[nc].cx);
    printf("y-coordinate of sphere center                 : %16.15g\n",ms->sp[nc].cy);
    printf("z-coordinate of sphere center                 : %16.15g\n",ms->sp[nc].cz);
    printf("basic sampling number on sphere surface       : %16d\n",ms->sp[nc].bsn);
    printf("division number for sphere surface    (per PI): %16d\n",ms->sp[nc].bdv);
    printf("limit of order number l                       : %16d\n",ms->sp[nc].l_limit);
  }
  printf("\n");
  
}

void setup_amsp(AMSP *ms)
{
  void check_position(AMSP *ms);
  void malloc_sp(ASP *sp);
  void setup_sp(ASP *sp,Maw *aw);
  void initialize_phii(ASP *sp,Maw *aw);
  void coefficient(ASP *sp);
  
  int i;
  
  check_position(ms);
  // multi_aw
  setup_maw(&(ms->aw));
  
  // spheres
  for(i=0;i<ms->n_sphr;i++){
    malloc_sp(&(ms->sp[i]));
    setup_sp(&(ms->sp[i]),&(ms->aw));
    initialize_phii(&(ms->sp[i]),&(ms->aw));
    coefficient(&(ms->sp[i]));
  } 
}

void free_amsp(AMSP *ms)
{
  void free_sp(ASP *sp);
  
  int i;
  
  // multi_aw
  free_maw(&(ms->aw));
  
  // spheres
  for(i=0;i<ms->n_sphr;i++){
    free_sp(&(ms->sp[i]));
  }
  free(ms->sp);  ms->n_sphr=0;
}

int sphere_id_type1(double *r,AMSP *ms)
{
  double r2;
  int i;
  
  for(i=0;i<ms->n_sphr;i++){
    r2=pow(r[0]-ms->sp[i].cx,2)+pow(r[1]-ms->sp[i].cy,2)+pow(r[2]-ms->sp[i].cz,2);
    if(ms->sp[i].a*ms->sp[i].a >= r2-MEPS) return i;
  }
  
  return -1;
}

int sphere_id_type2(double *r,AMSP *ms)
{
  double r2;
  int i;
  
  for(i=0;i<ms->n_sphr;i++){
    r2=pow(r[0]-ms->sp[i].cx,2)+pow(r[1]-ms->sp[i].cy,2)+pow(r[2]-ms->sp[i].cz,2);
    if(ms->sp[i].a*ms->sp[i].a > r2+MEPS) return i;
  }
  
  return -1;
}

void iterative_ops_amsp(AMSP *ms)
{
  void add_sfield(int src,int obj,AMSP *ms);
  void coefficient(ASP *sp); 
  
  int i,j,t,nn,sbc,num=ms->n_sphr;
  double vf[3],f1,*f0;
  int *bc;
  
  if(num<2) return;
  
  bc=(int *)m_alloc2(num,sizeof(int),"iterative_ops_amsp(),*bc");
  f0=(double *)m_alloc2(num,sizeof(double),"iterative_ops_amsp(),*f0");

  for(t=0;t<num;t++){
    F_amsp(t,vf,ms);
    f0[t]=vf[0]*vf[0]+vf[1]*vf[1]+vf[2]*vf[2];
    bc[t]=ito_breakcount;
  }

  printf("iterative operation start (convergence criterion : cv < %g)\n",ito_eps);
  for(nn=0;nn<ito_max;nn++){
    for(i=0;i<num;i++)
      for(j=0;j<num;j++) if(i!=j) add_sfield(i,j,ms);
    for(i=0;i<num;i++) coefficient(&(ms->sp[i]));
    printf("%3d, cv : ",nn);
    for(t=0;t<num;t++){
      F_amsp(t,vf,ms);
      f1=vf[0]*vf[0]+vf[1]*vf[1]+vf[2]*vf[2];
      if(fabs(f1/f0[t]-1.0)<ito_eps)  bc[t]--;
      printf("%g\t",fabs(f1/f0[t]-1.0));
      f0[t]=f1;
    }
    printf("\n");
    sbc=0;
    for(t=0;t<num;t++) if(bc[t]<=0) sbc++;
    if(sbc==num) break;
  }

  free(bc);  free(f0);
}

void output_node_particles(char *fname,AMSP *ms)
{
  FILE *fp;
  double a,st,ct,sp,cp,x,y,z;
  int s1,s2,oid,i,j;
  char *sd,fo[128]="",tmp[128]="";

  sd=strrchr(fname,'.');
  if(sd==NULL){ // no file extension
    sprintf(fo,"%s.particles",fname);
  }
  else {
    s1=strlen(fname);
    s2=strlen(sd);
    strncpy(tmp,fname,s1-s2);
    sprintf(fo,"%s.particles",tmp);
  }
  
  if((fp=fopen(fo,"wt"))==NULL){    printf("Can not open the %s file.\n",fo);    exit(1);  }
  fprintf(fp,"# x y z sphere_id\n");
  
  for(oid=0;oid<ms->n_sphr;oid++){
    a=ms->sp[oid].a;
    for(i=0;i<ms->sp[oid].dd.nt;i++){
      st=ms->sp[oid].dd.st[i];
      ct=ms->sp[oid].dd.ct[i];
      for(j=0;j<ms->sp[oid].dd.np;j++){
        sp=ms->sp[oid].dd.sp[j];
        cp=ms->sp[oid].dd.cp[j];

        x=a*st*cp+ms->sp[oid].cx;
        y=a*st*sp+ms->sp[oid].cy;
        z=a*ct   +ms->sp[oid].cz;
        fprintf(fp,"%15.14e %15.14e %15.14e %d\n",x,y,z,oid);
      }
    }
  }
  fclose(fp);
}

void write_dat_amsp(char *fn,AMSP *ms)
{
  FILE *fp;
  int i,mn,np,nt;
  
  if((fp=fopen(fn,"wb"))==NULL){    printf("write_dat_amsp(), Failed to open the file %s. Exit...\n",fn);    exit(1);  }
  
  fwrite(ms,sizeof(AMSP),1,fp);
  // beam data
  fwrite(ms->aw.bd.pw,sizeof(Apw),ms->aw.n_pw,fp);
  fwrite(ms->aw.bd.bb,sizeof(Abb),ms->aw.n_bb,fp);
  fwrite(ms->aw.bd.fb,sizeof(Afb),ms->aw.n_fb,fp);
  // sphere data  
  fwrite(ms->sp,sizeof(ASP),ms->n_sphr,fp);
  for(i=0;i<ms->n_sphr;i++){
    np=ms->sp[i].dd.np;
    nt=ms->sp[i].dd.nt;
    mn=ms->sp[i].l_limit;
    fwrite(ms->sp[i].dd.xt,sizeof(double),nt,fp);
    fwrite(ms->sp[i].dd.wt,sizeof(double),nt,fp);
    fwrite(ms->sp[i].dd.ct,sizeof(double),nt,fp);
    fwrite(ms->sp[i].dd.st,sizeof(double),nt,fp);
    fwrite(ms->sp[i].dd.xp,sizeof(double),np,fp);
    fwrite(ms->sp[i].dd.wp,sizeof(double),np,fp);
    fwrite(ms->sp[i].dd.cp,sizeof(double),np,fp);
    fwrite(ms->sp[i].dd.sp,sizeof(double),np,fp);
    fwrite(ms->sp[i].dd.cai,sizeof(double),mn+1,fp);
    fwrite(ms->sp[i].dd.cs,sizeof(double complex),mn+1,fp);
    fwrite(ms->sp[i].dd.cw,sizeof(double complex),mn+1,fp);
    fwrite(ms->sp[i].dd.ailm,sizeof(double complex),(mn+1)*(mn+1),fp);
  }
  fclose(fp);
}

void  read_dat_amsp(char *fn,AMSP *ms)
{
  void malloc_sp(ASP *sp);
  
  FILE *fp;
  int mn,i,nt,np;
  if((fp=fopen(fn,"rb"))==NULL){    printf("read_dat_amsp(), Failed to open the %s. Exit...\n",fn);    exit(1);  }
  
  fread(ms,sizeof(AMSP),1,fp);
  // beam data
  ms->aw.bd.pw=(Apw *)m_alloc2(ms->aw.n_pw,sizeof(Apw),"read_dat_amsp(),ms->aw.bd.pw");
  fread(ms->aw.bd.pw,sizeof(Apw),ms->aw.n_pw,fp);
  ms->aw.bd.bb=(Abb *)m_alloc2(ms->aw.n_bb,sizeof(Abb),"read_dat_amsp(),ms->aw.bd.bb");
  fread(ms->aw.bd.bb,sizeof(Abb),ms->aw.n_bb,fp);
  ms->aw.bd.fb=(Afb *)m_alloc2(ms->aw.n_fb,sizeof(Afb),"read_dat_amsp(),ms->aw.bd.fb");
  fread(ms->aw.bd.fb,sizeof(Afb),ms->aw.n_fb,fp);
  setup_maw(&(ms->aw));
  // sphere data
  ms->sp=(ASP *)m_alloc2(ms->n_sphr,sizeof(ASP),"read_dat_amsp(),ms->sp");
  fread(ms->sp,sizeof(ASP),ms->n_sphr,fp);
  for(i=0;i<ms->n_sphr;i++){
    malloc_sp(&(ms->sp[i]));
    np=ms->sp[i].dd.np;
    nt=ms->sp[i].dd.nt;
    mn=ms->sp[i].l_limit;
    fread(ms->sp[i].dd.xt,sizeof(double),nt,fp);
    fread(ms->sp[i].dd.wt,sizeof(double),nt,fp);
    fread(ms->sp[i].dd.ct,sizeof(double),nt,fp);
    fread(ms->sp[i].dd.st,sizeof(double),nt,fp);
    fread(ms->sp[i].dd.xp,sizeof(double),np,fp);
    fread(ms->sp[i].dd.wp,sizeof(double),np,fp);
    fread(ms->sp[i].dd.cp,sizeof(double),np,fp);
    fread(ms->sp[i].dd.sp,sizeof(double),np,fp);
    fread(ms->sp[i].dd.cai,sizeof(double),mn+1,fp);
    fread(ms->sp[i].dd.cs,sizeof(double complex),mn+1,fp);
    fread(ms->sp[i].dd.cw,sizeof(double complex),mn+1,fp);
    fread(ms->sp[i].dd.ailm,sizeof(double complex),(mn+1)*(mn+1),fp);
  }
  fclose(fp); 
}

//////////////////////////////////////////////////////////////////////
void check_position(AMSP *ms)
{
  double r,rs;
  int i,j;
  
  for(i=0;i<ms->n_sphr;i++){
    for(j=i+1;j<ms->n_sphr;j++){
      r =ms->sp[i].a+ms->sp[j].a;
      rs=sqrt(pow(ms->sp[j].cx-ms->sp[i].cx,2)+pow(ms->sp[j].cy-ms->sp[i].cy,2)+pow(ms->sp[j].cz-ms->sp[i].cz,2));
      if(rs<r){
        printf("Sphere Position Check Error! sphere id %d and sphere id %d is overlapped. Exit...\n",i,j);
        exit(1);
      }
    }
  }
}

void malloc_sp(ASP *sp)
{
  int ll;
  int nt=1*sp->bsn*sp->bdv;
  int np=2*sp->bsn*sp->bdv;
  
  // gauleg node and weight
  sp->dd.nt=nt;
  sp->dd.np=np;
  sp->dd.xt=(double *)m_alloc2(nt,sizeof(double),"malloc_sp(),sp->dd.xt");
  sp->dd.wt=(double *)m_alloc2(nt,sizeof(double),"malloc_sp(),sp->dd.wt");
  sp->dd.ct=(double *)m_alloc2(nt,sizeof(double),"malloc_sp(),sp->dd.ct");
  sp->dd.st=(double *)m_alloc2(nt,sizeof(double),"malloc_sp(),sp->dd.st");
  sp->dd.xp=(double *)m_alloc2(np,sizeof(double),"malloc_sp(),sp->dd.xp");
  sp->dd.wp=(double *)m_alloc2(np,sizeof(double),"malloc_sp(),sp->dd.wp");
  sp->dd.cp=(double *)m_alloc2(np,sizeof(double),"malloc_sp(),sp->dd.cp");
  sp->dd.sp=(double *)m_alloc2(np,sizeof(double),"malloc_sp(),sp->dd.sp");
  // incident field data on sphere surface
  sp->dd.phii=(double complex *)m_alloc2(nt*np,sizeof(double complex),"malloc_sp(), sp->dd.phii");
  sp->dd.phis=(double complex *)m_alloc2(nt*np,sizeof(double complex),"malloc_sp(), sp->dd.phis");
  // coefficient
  ll=sp->l_limit;
  sp->dd.cai=(double *)m_alloc2(ll+1,sizeof(double),"malloc_sp(),sp->dd.cai");
  sp->dd.cs=(double complex *)m_alloc2(ll+1,sizeof(double complex),"malloc_sp(),sp->dd.cs");
  sp->dd.cw=(double complex *)m_alloc2(ll+1,sizeof(double complex),"malloc_sp(),sp->dd.cw");
  sp->dd.ailm=(double complex *)m_alloc2((ll+1)*(ll+1),sizeof(double complex),"malloc_sp(),sp->dd.ailm");
}

void free_sp(ASP *sp)
{
  free(sp->dd.xt);  free(sp->dd.wt);
  free(sp->dd.ct);  free(sp->dd.st);
  free(sp->dd.xp);  free(sp->dd.wp);
  free(sp->dd.cp);  free(sp->dd.sp);
  
  free(sp->dd.phii);  free(sp->dd.phis);
  
  free(sp->dd.cai);
  free(sp->dd.cs);  free(sp->dd.cw);
  free(sp->dd.ailm);
}


void setup_sp(ASP *sp,Maw *aw)
{
  void gauleg_dv(double a,double b,double *x,double *w,int bn,int dv);
  
  double complex *h1l,*dh1l,tz;
  double kea,i_kea,kia,i_kia,bk,bk2,cfe,cfi,omega;
  double *jle,*djle,*yle,*jli,*djli;
  int ll,nt,np,i;

  ll=sp->l_limit;
  nt=sp->dd.nt;
  np=sp->dd.np;
  
  //gauleg node and weight
  gauleg_dv(0.0,    M_PI,sp->dd.xt,sp->dd.wt,sp->bsn,  sp->bdv);
  gauleg_dv(0.0,2.0*M_PI,sp->dd.xp,sp->dd.wp,sp->bsn,2*sp->bdv);
  for(i=0;i<nt;i++){
    sp->dd.st[i]=sin(sp->dd.xt[i]);
    sp->dd.ct[i]=cos(sp->dd.xt[i]);
  }
  for(i=0;i<np;i++){
    sp->dd.sp[i]=sin(sp->dd.xp[i]);
    sp->dd.cp[i]=cos(sp->dd.xp[i]);
  }

  // constant
  omega=2.0*M_PI*aw->f;
  sp->dd.K0=sp->rho0*sp->c0*sp->c0;
  sp->dd.k0=omega/sp->c0;
  sp->dd.k1=I*omega/sp->dd.K0;
  sp->dd.k2=I*omega*sp->rho0;
  
  // coefficient
  kea=aw->k0*sp->a;
  i_kea=1.0/kea;
  kia=sp->dd.k0*sp->a;
  i_kia=1.0/kia;
  bk=sp->dd.k0/aw->k0;
  bk2=sp->rho0/aw->rho0;
  jle=(double *)m_alloc2(ll+2,sizeof(double),"setup_sp(), jle");
  yle=(double *)m_alloc2(ll+2,sizeof(double),"setup_sp(), yle");
  jli=(double *)m_alloc2(ll+2,sizeof(double),"setup_sp(), jli");
  gsl_sf_bessel_jl_array(ll+1,kea,jle);
  gsl_sf_bessel_yl_array(ll+1,kea,yle);
  gsl_sf_bessel_jl_array(ll+1,kia,jli);
  h1l =(double complex *)m_alloc2(ll+2,sizeof(double complex),"setup_sp(), h1l");
  for(i=0;i<=ll+1;i++) h1l[i]=jle[i]+I*yle[i];
  djle=(double *)m_alloc2(ll+1,sizeof(double),"setup_sp(), djle");
  djli=(double *)m_alloc2(ll+1,sizeof(double),"setup_sp(), djli");
  dh1l=(double complex *)m_alloc2(ll+1,sizeof(double complex),"setup_sp(), dh1l");
  for(i=0;i<=ll;i++){
    cfe=(double)i*i_kea;
    cfi=(double)i*i_kia;
    djle[i]=cfe*jle[i]-jle[i+1];
    djli[i]=cfi*jli[i]-jli[i+1];
    dh1l[i]=cfe*h1l[i]-h1l[i+1];
  }
  for(i=0;i<=ll;i++){
    sp->dd.cai[i]=1.0/jle[i];
    tz=1.0/(bk2*jli[i]*dh1l[i]-bk*djli[i]*h1l[i]);
    sp->dd.cw[i]=(jle[i]*dh1l[i]-djle[i]*h1l[i])*tz;
    sp->dd.cs[i]=(bk*jle[i]*djli[i]-bk2*djle[i]*jli[i])*tz;
  }
  
  free(jle);  free(yle);  free(jli);  free(h1l);
  free(djle); free(djli); free(dh1l);
}

void gauleg_dv(double a,double b,double *x,double *w,int bn,int dv)
{
  double xt[bn],wt[bn];
  
  gauleg(-1.0, 1.0,xt,wt,bn);
  
  double h,dh,x0,x1,cx,cc;
  int d,i,j;
  h=b-a;
  dh=h/(double)dv;
  x1=a;
  j=0;
  for(d=0;d<dv;d++){
    x0=x1;
    x1=x0+dh;
    
    cx=0.5*(x1-x0);
    cc=0.5*(x1+x0);
    for(i=0;i<bn;i++){
      x[j]= cx*xt[i]+cc;
      w[j]= cx*wt[i];
      j++;
    }
  }
}

void initialize_phii(ASP *sp,Maw *aw)
{
  double complex p,v[3],cp;
  double r,x[3];
  int i,j,nt,np;
  
  nt=sp->dd.nt;
  np=sp->dd.np;
  r=sp->a;
  cp=-1.0/aw->k2;
  for(i=0;i<nt;i++){
    #pragma omp parallel for schedule(dynamic) private(x,p,v) // parallel for
    for(j=0;j<np;j++){
      x[0]=r*sp->dd.st[i]*sp->dd.cp[j]+sp->cx;
      x[1]=r*sp->dd.st[i]*sp->dd.sp[j]+sp->cy;
      x[2]=r*sp->dd.ct[i]+sp->cz;
      calc_maw_pv(&p,v,x,aw);
      sp->dd.phii[i*np+j]=cp*p;
    }
  }
}

void coefficient(ASP *sp)
{
  int ti,lm,np,nt,tt,l,m,i,j,t;
  size_t ms,lmax;
  
  lmax=(size_t)sp->l_limit;
  ms=gsl_sf_legendre_array_n(lmax);
  lm=(int)lmax;
  np=sp->dd.np;
  nt=sp->dd.nt;

  for(ti=0;ti<(lm+1)*(lm+1);ti++) sp->dd.ailm[ti]=0.0;
  
  #pragma omp parallel private(i,j,l,m,t,tt,ti) // parallel 
  {
  double complex Yp,Ym;
  double sin_t,cos_t,sig;
  
  double *sphPlm=(double *)m_alloc2(ms,sizeof(double),"coefficient(),*sphPlm");
  double complex *e_phim=(double complex *)m_alloc2(np,sizeof(double complex),"coefficient(),*e_phim");
  double complex *tlm=(double complex *)m_alloc2((lm+1)*(lm+1),sizeof(double complex),"coefficient(),*tlm");
  
  #pragma omp for schedule(dynamic) // parallel for
  for(i=0;i<nt;i++){ // theta
    sin_t=sp->dd.st[i];    
    cos_t=sp->dd.ct[i];
    tt=0;
    m=0;
    gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM,lmax,cos_t,-1,sphPlm);
    for(l=0;l<=lm;l++){
      Yp=sphPlm[gsl_sf_legendre_array_index(l,m)];
      for(j=0;j<np;j++){ // phi
        tlm[tt]+=(sp->dd.phii[i*np+j]+sp->dd.phis[i*np+j])*conj(Yp)*sp->dd.wp[j];
      }
      tt++;
    }
    sig=1.0;
    for(t=0;t<np;t++) e_phim[t]=1.0;
    for(m=1;m<=lm;m++){
      sig*=-1.0;
      for(t=0;t<np;t++) e_phim[t]*=sp->dd.cp[t]+I*sp->dd.sp[t];
      for(l=m;l<=lm;l++){
        for(j=0;j<np;j++){ // phi
          Yp=sphPlm[gsl_sf_legendre_array_index(l,m)]*e_phim[j];
          Ym=sig*conj(Yp);
          tlm[tt+0]+=(sp->dd.phii[i*np+j]+sp->dd.phis[i*np+j])*conj(Yp)*sp->dd.wp[j];
          tlm[tt+1]+=(sp->dd.phii[i*np+j]+sp->dd.phis[i*np+j])*conj(Ym)*sp->dd.wp[j];
        }
        tt+=2;
      }
    }
    for(ti=0;ti<(lm+1)*(lm+1);ti++){
      #pragma omp critical // omp critical
      sp->dd.ailm[ti]+=tlm[ti]*sin_t*sp->dd.wt[i];
      tlm[ti]=0.0;
    }
  }
  free(sphPlm);
  free(e_phim);
  free(tlm);
  } // end parallel
  
  tt=0;  m=0;
  for(l=0;l<=lm;l++){
    sp->dd.ailm[tt]*=sp->dd.cai[l]; 
    tt++;
  }
  for(m=1;m<=lm;m++){
    for(l=m;l<=lm;l++){
      sp->dd.ailm[tt  ]*=sp->dd.cai[l];
      sp->dd.ailm[tt+1]*=sp->dd.cai[l];
      tt+=2;
    }
  }
  for(i=0;i<nt;i++){
    for(j=0;j<np;j++){
      sp->dd.phis[i*np+j]=0.0;
      sp->dd.phis[i*np+j]=0.0;
    }
  } 
}

double complex pickup_ailm(int l,int m,ASP *sp)
{
  int lm=sp->l_limit;
  int am=abs(m);
  int ll;
  
  if(l<0) return 0.0;
  else if(l<am) return 0.0;
  else {
    if(m==0) return sp->dd.ailm[l];
    else {
      ll=(lm+1)*(2*am-1)-am*(am-1)+2*(l-am);
      if (m>0) return sp->dd.ailm[ll+0];
      else     return sp->dd.ailm[ll+1];
    }
  }
}

void add_sfield(int src,int obj,AMSP *ms)
{
  int scattered_phi_ss(double complex *phi,double *xb,ASP *sp,Maw *aw); // aw_msp_ivf_field.c
  
  double complex ps;
  double cos_t,sin_t,cos_p,sin_p;
  double r[3];
  double a=ms->sp[obj].a;
  int np=ms->sp[obj].dd.np;
  int nt=ms->sp[obj].dd.nt;
  int i,j;

  #pragma omp parallel for schedule(dynamic) private(j,cos_t,sin_t,cos_p,sin_p,r,ps)
  for(i=0;i<nt;i++){
    cos_t=ms->sp[obj].dd.ct[i];    sin_t=ms->sp[obj].dd.st[i];
    for(j=0;j<np;j++){
      cos_p=ms->sp[obj].dd.cp[j];      sin_p=ms->sp[obj].dd.sp[j];
      r[0]=a*sin_t*cos_p+ms->sp[obj].cx;
      r[1]=a*sin_t*sin_p+ms->sp[obj].cy;
      r[2]=a*cos_t      +ms->sp[obj].cz;
      scattered_phi_ss(&ps,r,&(ms->sp[src]),&(ms->aw));
      ms->sp[obj].dd.phis[i*np+j]+=ps;
    }
  }
}
