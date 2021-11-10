/* Analysis of sound wave scattering by a spherical object 
 * 
 * aw_msp_ivf.h
 *
 *  Created on: Nov 01, 2021
 *      Author: ohta
 */
#ifndef AW_MSP_IVF_H_
#define AW_MSP_IVF_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "my_utils.h"
#include "multi_aw.h"
#include "gsl/gsl_specfunc.h"

// machine epspilon
#define MEPS 1.0e-15

// default sphere datafile name
#define fn_sphr "msphr.txt"

// for iterative operations
#define ito_max 256
#define ito_eps 1.0e-12
#define ito_breakcount 2

// struct definition
typedef struct discrete_data{
  double K0;                     // bulk modulus [N/m^2], K0 = c0^2*rho0.
  double k0;                     // wavenumber, k0=omega/c0.
  double complex k1,k2;          // k1=I*omega/K_0, k2=I*omega*rho0, k0^2=-k1*k2.
  
  int    nt,np;                  // number of sampling points.
  double *xt,*wt,*xp,*wp;        // gauss-legender node and weight data.
  double *st,*ct,*sp,*cp;        // sin, cos value at the node.
  double complex *phii;          // velocity potential of incident field on sphere surface.  
  double complex *phis;          // velocity potential of scattered field on sphere surface (for iterative solution).
  
  double *cai;                   // coefficient for spherical harmonic expansion.
  double complex *cs,*cw;        // coefficient for scattered and internal field.
  double complex *ailm;          // coefficient of incident field.
}DDT;

typedef struct sphere_data{
  double rho0;          // density of the sphere objcet [kg/m^3].
  double c0;            // speed of sound in the sphere object [m/s].
  double a;             // radius of the sphere object [m].
  double cx,cy,cz;      // sphere center (cx,cy,cz) [m].
  int bsn;              // basic sampling number for Gauss-Legendre quadrature.
  int bdv;              // division number of sphere surface ( per M_PI ).
  int l_limit;          // limit of order number.
  
  DDT dd;               // discrete data.
}ASP;

typedef struct sphere_objects{
  int  n_sphr;          // number of spheres.
  ASP *sp;              // sphere data.
  Maw aw;               // incident field (multi_aw) data.
}AMSP;


// -- aw_msp_ivf.c - -
void read_data_amsp(AMSP *ms);           // search the sphere datafile (defined fn_sphr) and load.
void print_data_amsp(AMSP *ms);          // print loaded data.
void setup_amsp(AMSP *ms);               // allocate memory and setup coefficients.
void  free_amsp(AMSP *ms);               // free allocated memory.
int sphere_id_type1(double *r,AMSP *ms); // return sphere id at r. the boundary is determined to be inside the sphere. return -1 if outside the spheres. 
int sphere_id_type2(double *r,AMSP *ms); // return sphere id at r. the boundary is determined to be outside the sphere.
void iterative_ops_amsp(AMSP *ms);       // solve multiple scattering.
void output_node_particles(char *fname,AMSP *ms); // outputs the nodes for surface integral as point cloud data ( .particles file ).
void write_dat_amsp(char *fn,AMSP *ms);  // write datafile.
void  read_dat_amsp(char *fn,AMSP *ms);  // read datafile and allocate memory.

// -- aw_msp_ivf_force.c --
void  F_amsp(int id,double *f,AMSP *ms);                      // calculate radiation force. 
void FN_amsp(int id,double *f,double *n,double *rc,AMSP *ms); // calculate radiation force and torque.
void FN_sfint_amsp(int id,double *f,double *n,double *rc,AMSP *ms); // surface integral of radiation stress tensor on the boundary (test function).
// inputs 
// id : sphere id. 
// rc : coordinate of rotation center, rc[0]=x, rc[1]=y, rc[2]=z.
// ms : pointer of AMSP object 
// outputs
// f  : radiation force,  f[0]=f_x, f[1]=f_y, f[2]=f_z.
// n  : radiation torque, n[0]=n_x, n[1]=n_y, n[2]=n_z.


// -- aw_msp_ivf_field.c
int  incident_pv_amsp(double complex *p,double complex *v,double *r,AMSP *ms); // calculate incident field 
int  internal_pv_amsp(double complex *p,double complex *v,double *r,AMSP *ms); // calculate internal field 
int scattered_pv_amsp(double complex *p,double complex *v,double *r,AMSP *ms); // calculate scattered field
int     total_pv_amsp(double complex *p,double complex *v,double *r,int type,AMSP *ms); // calculate total field
// imputs
// r   : coordinate of calculation point, r[0]=x, r[1]=y, r[2]=z.
// type: select sphere_id_* function, type=1 : sphere_id_type1()(select inside field at the boundary), type!=1 : sphere_id_thpe2() (select outside field at the boundary). 
// ms  : pointer of AMSP object.
// outputs
// p   : sound pressure.
// v   : particle velocity, v[0]=v_x, v[1]=v_y, v[2]=v_z.
// return 
// sphere id. return -1 if outside the spheres.

#endif
