#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../c/include/tool.h"

static size_t idx3(int i,int j,int k,const int ex[3]){return (size_t)i+ (size_t)ex[0]*((size_t)j+(size_t)ex[1]*(size_t)k);} 

static void fill_grid(double *a,int n,double h){for(int i=0;i<n;++i)a[i]=i*h;}
static int test_fderivs(){
  int ex[3]={20,18,16}; size_t n=(size_t)ex[0]*ex[1]*ex[2];
  double *X=calloc(ex[0],8),*Y=calloc(ex[1],8),*Z=calloc(ex[2],8),*f=calloc(n,8),*fx=calloc(n,8),*fy=calloc(n,8),*fz=calloc(n,8);
  fill_grid(X,ex[0],0.05); fill_grid(Y,ex[1],0.04); fill_grid(Z,ex[2],0.03);
  for(int k=0;k<ex[2];++k)for(int j=0;j<ex[1];++j)for(int i=0;i<ex[0];++i){size_t p=idx3(i,j,k,ex); double x=X[i],y=Y[j],z=Z[k]; f[p]=sin(x)+cos(2*y)+z*z;}
  fderivs(ex,f,fx,fy,fz,X,Y,Z,1,1,1,0,1);
  double err=0; int cnt=0;
  for(int k=3;k<ex[2]-3;++k)for(int j=3;j<ex[1]-3;++j)for(int i=3;i<ex[0]-3;++i){size_t p=idx3(i,j,k,ex); double x=X[i],y=Y[j],z=Z[k];
    err=fmax(err,fabs(fx[p]-cos(x))); err=fmax(err,fabs(fy[p]-(-2*sin(2*y)))); err=fmax(err,fabs(fz[p]-(2*z))); cnt++;}
  printf("fderivs interior max err=%.3e (%d points)\n",err,cnt);
  return err<5e-4?0:1;
}

static int test_fdderivs(){
  int ex[3]={20,19,18}; size_t n=(size_t)ex[0]*ex[1]*ex[2];
  double *X=calloc(ex[0],8),*Y=calloc(ex[1],8),*Z=calloc(ex[2],8),*f=calloc(n,8);
  double *fxx=calloc(n,8),*fxy=calloc(n,8),*fxz=calloc(n,8),*fyy=calloc(n,8),*fyz=calloc(n,8),*fzz=calloc(n,8);
  fill_grid(X,ex[0],0.05); fill_grid(Y,ex[1],0.05); fill_grid(Z,ex[2],0.05);
  for(int k=0;k<ex[2];++k)for(int j=0;j<ex[1];++j)for(int i=0;i<ex[0];++i){size_t p=idx3(i,j,k,ex); double x=X[i],y=Y[j],z=Z[k]; f[p]=x*x+3*x*y+2*y*z+sin(z);}  
  fdderivs(ex,f,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,1,1,1,0,1);
  double err=0; int cnt=0;
  for(int k=3;k<ex[2]-3;++k)for(int j=3;j<ex[1]-3;++j)for(int i=3;i<ex[0]-3;++i){size_t p=idx3(i,j,k,ex); double z=Z[k];
    err=fmax(err,fabs(fxx[p]-2.0)); err=fmax(err,fabs(fxy[p]-3.0)); err=fmax(err,fabs(fxz[p]-0.0));
    err=fmax(err,fabs(fyy[p]-0.0)); err=fmax(err,fabs(fyz[p]-2.0)); err=fmax(err,fabs(fzz[p]-(-sin(z)))); cnt++;}
  printf("fdderivs interior max err=%.3e (%d points)\n",err,cnt);
  return err<8e-4?0:1;
}

static int test_kodis_constant(){
  int ex[3]={12,11,10}; size_t n=(size_t)ex[0]*ex[1]*ex[2];
  double *X=calloc(ex[0],8),*Y=calloc(ex[1],8),*Z=calloc(ex[2],8),*f=calloc(n,8),*rhs=calloc(n,8),*rhs0=calloc(n,8);
  fill_grid(X,ex[0],0.1); fill_grid(Y,ex[1],0.1); fill_grid(Z,ex[2],0.1);
  for(size_t i=0;i<n;++i){f[i]=2.5; rhs[i]=sin(i*0.01); rhs0[i]=rhs[i];}
  double soa[3]={1,1,1}; kodis(ex,X,Y,Z,f,rhs,soa,0,0.1);
  double err=0; for(size_t i=0;i<n;++i) err=fmax(err,fabs(rhs[i]-rhs0[i]));
  printf("kodis constant invariance err=%.3e\n",err);
  return err<1e-12?0:1;
}

static int test_lopsided_zero_shift(){
  int ex[3]={12,11,10}; size_t n=(size_t)ex[0]*ex[1]*ex[2];
  double *X=calloc(ex[0],8),*Y=calloc(ex[1],8),*Z=calloc(ex[2],8),*f=calloc(n,8),*rhs=calloc(n,8),*rhs0=calloc(n,8),*sfx=calloc(n,8),*sfy=calloc(n,8),*sfz=calloc(n,8);
  fill_grid(X,ex[0],0.1); fill_grid(Y,ex[1],0.1); fill_grid(Z,ex[2],0.1);
  for(size_t i=0;i<n;++i){f[i]=cos(i*0.02); rhs[i]=sin(i*0.03); rhs0[i]=rhs[i];}
  double soa[3]={1,1,1}; lopsided(ex,X,Y,Z,f,rhs,sfx,sfy,sfz,0,soa);
  double err=0; for(size_t i=0;i<n;++i) err=fmax(err,fabs(rhs[i]-rhs0[i]));
  printf("lopsided zero-shift invariance err=%.3e\n",err);
  return err<1e-12?0:1;
}

int main(){int rc=0; rc|=test_fderivs(); rc|=test_fdderivs(); rc|=test_kodis_constant(); rc|=test_lopsided_zero_shift(); return rc;}
