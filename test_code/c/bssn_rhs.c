#include "include/bssn_rhs_compute.h"
// 0-based i,j,k
// #define IDX_F(i,j,k,nx,ny) ((i) + (j)*(nx) + (k)*(nx)*(ny))
// ex(1)=nx, ex(2)=ny, ex(3)=nz

// 用法：a[ IDX_F(i,j,k,nx,ny) ]

// C function that calculates the right-hand side for BSSN equations
int f_compute_rhs_bssn(int *ex, double &T, 
                       double *X, double *Y, double *Z,
                       double *chi, double *trK,
                       double *dxx, double *gxy, double *gxz, double *dyy, double *gyz, double *dzz,
                       double *Axx, double *Axy, double *Axz, double *Ayy, double *Ayz, double *Azz,
                       double *Gamx, double *Gamy, double *Gamz,
                       double *Lap, double *betax, double *betay, double *betaz,
                       double *dtSfx, double *dtSfy, double *dtSfz,
                       double *chi_rhs, double *trK_rhs,
                       double *gxx_rhs, double *gxy_rhs, double *gxz_rhs, double *gyy_rhs, double *gyz_rhs, double *gzz_rhs,
                       double *Axx_rhs, double *Axy_rhs, double *Axz_rhs, double *Ayy_rhs, double *Ayz_rhs, double *Azz_rhs,
                       double *Gamx_rhs, double *Gamy_rhs, double *Gamz_rhs,
                       double *Lap_rhs, double *betax_rhs, double *betay_rhs, double *betaz_rhs,
                       double *dtSfx_rhs, double *dtSfy_rhs, double *dtSfz_rhs,
                       double *rho, double *Sx, double *Sy, double *Sz,
                       double *Sxx, double *Sxy, double *Sxz, double *Syy, double *Syz, double *Szz,
                       double *Gamxxx, double *Gamxxy, double *Gamxxz, double *Gamxyy, double *Gamxyz, double *Gamxzz,
                       double *Gamyxx, double *Gamyxy, double *Gamyxz, double *Gamyyy, double *Gamyyz, double *Gamyzz,
                       double *Gamzxx, double *Gamzxy, double *Gamzxz, double *Gamzyy, double *Gamzyz, double *Gamzzz,
                       double *Rxx, double *Rxy, double *Rxz, double *Ryy, double *Ryz, double *Rzz,
                       double *ham_Res, double *movx_Res, double *movy_Res, double *movz_Res,
                       double *Gmx_Res, double *Gmy_Res, double *Gmz_Res,
                       int &Symmetry, int &Lev, double &eps, int &co
                       )  // return gont
{   
    int nx = ex[0], ny = ex[1], nz=ex[2];
    int all = nx*ny*nz;

    // input
    // const int ex[3], Symmetry, Lev, co;
    // const double T;
    // const double X[ex[0]], Y[ex[1]], Z[ex[2]];
    // const double trK[all];
    // const double gxy[all], gxz[all], gyz[all];
    // const double Axx[all], Axy[all], Axz[all], Ayy[all], Ayz[all], Azz[all];
    // const double Gamx[all], Gamy[all], Gamz[all];
    // const double dtSfx[all], dtSfy[all], dtSfz[all];
    // const double rho[all], Sx[all], Sy[all], Sz[all];
    // const double Sxx[all], Sxy[all], Sxz[all], Syy[all], Syz[all], Szz[all];
    // const double eps;
    
    // double Lap[all], betax[all], betay[all], betaz[all];
    // double chi[all], dxx[all], dyy[all], dzz[all];
    // double chi_rhs[all], trK_rhs[all];
    // double gxx_rhs[all], gxy_rhs[all], gxz_rhs[all];
    // double gyy_rhs[all], gyz_rhs[all], gzz_rhs[all];
    // double Axx_rhs[all], Axy_rhs[all], Axz_rhs[all];
    // double Ayy_rhs[all], Ayz_rhs[all], Azz_rhs[all];
    // double Gamx_rhs[all], Gamy_rhs[all], Gamz_rhs[all];
    // double Lap_rhs[all], betax_rhs[all], betay_rhs[all], betaz_rhs[all];
    // double dtSfx_rhs[all], dtSfy_rhs[all], dtSfz_rhs[all];
    // double Gamxxx[all], Gamxxy[all], Gamxxz[all];
    // double Gamxyy[all], Gamxyz[all], Gamxzz[all];
    // double Gamyxx[all], Gamyxy[all], Gamyxz[all];
    // double Gamyyy[all], Gamyyz[all], Gamyzz[all];
    // double Gamzxx[all], Gamzxy[all], Gamzxz[all];
    // double Gamzyy[all], Gamzyz[all], Gamzzz[all];
    // double Rxx[all], Rxy[all], Rxz[all], Ryy[all], Ryz[all], Rzz[all];
    // double ham_Res[all], movx_Res[all], movy_Res[all], movz_Res[all];
    // double Gmx_Res[all], Gmy_Res[all], Gmz_Res[all];

    // temp variable
    double gxx[all],gyy[all],gzz[all];
    double chix[all],chiy[all],chiz[all];
    double gxxx[all],gxyx[all],gxzx[all],gyyx[all],gyzx[all],gzzx[all];
    double gxxy[all],gxyy[all],gxzy[all],gyyy[all],gyzy[all],gzzy[all];
    double gxxz[all],gxyz[all],gxzz[all],gyyz[all],gyzz[all],gzzz[all];
    double Lapx[all], Lapy[all], Lapz[all];
    double betaxx[all], betaxy[all], betaxz[all];
    double betayx[all], betayy[all], betayz[all];
    double betazx[all], betazy[all], betazz[all];
    double Gamxx[all],Gamxy[all],Gamxz[all];
    double Gamyx[all],Gamyy[all],Gamyz[all];
    double Gamzx[all],Gamzy[all],Gamzz[all];
    double Kx[all], Ky[all], Kz[all], div_beta[all], S[all];
    double f[all], fxx[all], fxy[all], fxz[all], fyy[all], fyz[all], fzz[all];
    double Gamxa[all], Gamya[all], Gamza[all], alpn1[all], chin1[all];
    double gupxx[all], gupxy[all], gupxz[all];
    double gupyy[all], gupyz[all], gupzz[all];
    double SSS[3],AAS[3],ASA[3],SAA[3],ASS[3],SAS[3],SSA[3];
    double dX, dY, dZ, PI;
    const double ZEO = 0, ONE = 1, TWO = 2, FOUR = 4;
    const double EIGHT = 8, HALF = 0.5, THR = 3;
    const double SYM = 1, ANTI = -1;
    const double FF = 0.75, eta = 2;
    const double F1o3 = 1/3, F2o3 = 2/3, F3o2 = 1.5, F1o6 = 1/6;
    const double F16 = 16, F8 = 8;
    #if (GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5)
    double reta[all];
    /* 使用时：reta[idx]，其中 idx = i + nx*(j + ny*k)  (Fortran列主序) */
    #endif
    PI = acos(-1.0); 
    dX = X[1] - X[0];
    dY = Y[1] - Y[0];
    dZ = Z[1] - Z[0];

    for(int i=0;i<all;i+=1){
        alpn1[i] = Lap[i] + 1;
    }

    for(int i=0;i<all;i+=1){
        chin1[i] = chi[i] + 1;
    }

    for(int i=0;i<all;i+=1){
        gxx[i] = dxx[i] + 1;
    }

    for(int i=0;i<all;i+=1){
        gyy[i] = dyy[i] + 1;
    }

    for(int i=0;i<all;i+=1){
        gzz[i] = dzz[i] + 1;
    }

    fderivs(ex,betax,betaxx,betaxy,betaxz,X,Y,Z,ANTI, SYM, SYM,Symmetry,Lev);
    fderivs(ex,betay,betayx,betayy,betayz,X,Y,Z, SYM,ANTI, SYM,Symmetry,Lev);
    fderivs(ex,betaz,betazx,betazy,betazz,X,Y,Z, SYM, SYM,ANTI,Symmetry,Lev);

    for(int i=0;i<all;i+=1){
        div_beta[i] = betaxx[i] + betayy[i] + betazz[i];
    }

    fderivs(ex,chi,chix,chiy,chiz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev);
    for(int i=0;i<all;i+=1){
        chi_rhs[i] = F2o3 * chin1[i] * (alpn1[i] * trK[i] - div_beta[i]); 
    }
    fderivs(ex,dxx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev);
    fderivs(ex,gxy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,Lev);
    fderivs(ex,gxz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,Lev);
    fderivs(ex,dyy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev);
    fderivs(ex,gyz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,Lev);
    fderivs(ex,dzz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev);
    
    for(int i=0;i<all;i+=1){
        gxx_rhs[i] = -TWO * alpn1[i] * Axx[i] - F2o3 * gxx[i] * div_beta[i] + 
                    TWO * (gxx[i] * betaxx[i] + gxy[i] * betayx[i] + gxz[i] * betazx[i]);
    }

    for(int i=0;i<all;i+=1){
        gyy_rhs[i] = -TWO * alpn1[i] * Ayy[i] - F2o3 * gyy[i] * div_beta[i] + 
                    TWO * (gxy[i] * betaxy[i] + gyy[i] * betayy[i] + gyz[i] * betazy[i]);
    }

    for(int i=0;i<all;i+=1){
        gzz_rhs[i] = -TWO * alpn1[i] * Azz[i] - F2o3 * gzz[i] * div_beta[i] + 
                    TWO * (gxz[i] * betaxz[i] + gyz[i] * betayz[i] + gzz[i] * betazz[i]);
    }

    for(int i=0;i<all;i+=1){
        gxy_rhs[i] = -TWO * alpn1[i] * Axy[i] + F1o3 * gxy[i] * div_beta[i] + 
                    gxx[i] * betaxy[i] + gxz[i] * betazy[i] + gyy[i] * betayx[i] 
                    + gyz[i] * betazx[i] - gxy[i] * betazz[i];
    }

    for(int i=0;i<all;i+=1){
        gyz_rhs[i] = -TWO * alpn1[i] * Ayz[i] + F1o3 * gyz[i] * div_beta[i] + 
                    gxy[i] * betaxz[i] + gyy[i] * betayz[i] + gxz[i] * betaxy[i] 
                    + gzz[i] * betazy[i] - gyz[i] * betaxx[i];
    }

    for(int i=0;i<all;i+=1){
        gxz_rhs[i] = -TWO * alpn1[i] * Axz[i] + F1o3 * gxz[i] * div_beta[i] + 
                    gxx[i] * betaxz[i] + gxy[i] * betayz[i] + gyz[i] * betayx[i] 
                    + gzz[i] * betazx[i] - gxz[i] * betayy[i];
    }

    for(int i=0;i<all;i+=1){
        gupzz[i] = gxx[i] * gyy[i] * gzz[i] + gxy[i] * gyz[i] * gxz[i] + gxz[i] * gxy[i] * gyz[i] - 
                    gxz[i] * gyy[i] * gxz[i] - gxy[i] * gxy[i] * gzz[i] - gxx[i] * gyz[i] * gyz[i];
    }
    for(int i=0;i<all;i+=1){
        gupxx[i] = (gyy[i] * gzz[i] - gyz[i] * gyz[i]) / gupzz[i];
        gupxy[i] = -(gxy[i] * gzz[i] - gyz[i] * gxz[i]) / gupzz[i];
        gupxz[i] = (gxy[i] * gyz[i] - gyy[i] * gxz[i]) / gupzz[i];
        gupyy[i] = (gxx[i] * gzz[i] - gxz[i] * gxz[i]) / gupzz[i];
        gupyz[i] = -(gxx[i] * gyz[i] - gxy[i] * gxz[i]) / gupzz[i];
        gupzz[i] = (gxx[i] * gyy[i] - gxy[i] * gxy[i]) / gupzz[i];
    }
    if(co==0){
        for (int i = 0; i < all; i += 1) {
            Gmx_Res[i] = Gamx[i] - (
                gupxx[i] * (gupxx[i]*gxxx[i] + gupxy[i]*gxyx[i] + gupxz[i]*gxzx[i]) +
                gupxy[i] * (gupxx[i]*gxyx[i] + gupxy[i]*gyyx[i] + gupxz[i]*gyzx[i]) +
                gupxz[i] * (gupxx[i]*gxzx[i] + gupxy[i]*gyzx[i] + gupxz[i]*gzzx[i]) +

                gupxx[i] * (gupxy[i]*gxxy[i] + gupyy[i]*gxyy[i] + gupyz[i]*gxzy[i]) +
                gupxy[i] * (gupxy[i]*gxyy[i] + gupyy[i]*gyyy[i] + gupyz[i]*gyzy[i]) +
                gupxz[i] * (gupxy[i]*gxzy[i] + gupyy[i]*gyzy[i] + gupyz[i]*gzzy[i]) +

                gupxx[i] * (gupxz[i]*gxxz[i] + gupyz[i]*gxyz[i] + gupzz[i]*gxzz[i]) +
                gupxy[i] * (gupxz[i]*gxyz[i] + gupyz[i]*gyyz[i] + gupzz[i]*gyzz[i]) +
                gupxz[i] * (gupxz[i]*gxzz[i] + gupyz[i]*gyzz[i] + gupzz[i]*gzzz[i])
            );

            Gmy_Res[i] = Gamy[i] - (
                gupxx[i] * (gupxy[i]*gxxx[i] + gupyy[i]*gxyx[i] + gupyz[i]*gxzx[i]) +
                gupxy[i] * (gupxy[i]*gxyx[i] + gupyy[i]*gyyx[i] + gupyz[i]*gyzx[i]) +
                gupxz[i] * (gupxy[i]*gxzx[i] + gupyy[i]*gyzx[i] + gupyz[i]*gzzx[i]) +

                gupxy[i] * (gupxy[i]*gxxy[i] + gupyy[i]*gxyy[i] + gupyz[i]*gxzy[i]) +
                gupyy[i] * (gupxy[i]*gxyy[i] + gupyy[i]*gyyy[i] + gupyz[i]*gyzy[i]) +
                gupyz[i] * (gupxy[i]*gxzy[i] + gupyy[i]*gyzy[i] + gupyz[i]*gzzy[i]) +

                gupxy[i] * (gupxz[i]*gxxz[i] + gupyz[i]*gxyz[i] + gupzz[i]*gxzz[i]) +
                gupyy[i] * (gupxz[i]*gxyz[i] + gupyz[i]*gyyz[i] + gupzz[i]*gyzz[i]) +
                gupyz[i] * (gupxz[i]*gxzz[i] + gupyz[i]*gyzz[i] + gupzz[i]*gzzz[i])
            );

            Gmz_Res[i] = Gamz[i] - (
                gupxx[i] * (gupxz[i]*gxxx[i] + gupyz[i]*gxyx[i] + gupzz[i]*gxzx[i]) +
                gupxy[i] * (gupxz[i]*gxyx[i] + gupyz[i]*gyyx[i] + gupzz[i]*gyzx[i]) +
                gupxz[i] * (gupxz[i]*gxzx[i] + gupyz[i]*gyzx[i] + gupzz[i]*gzzx[i]) +

                gupxy[i] * (gupxz[i]*gxxy[i] + gupyz[i]*gxyy[i] + gupzz[i]*gxzy[i]) +
                gupyy[i] * (gupxz[i]*gxyy[i] + gupyz[i]*gyyy[i] + gupzz[i]*gyzy[i]) +
                gupyz[i] * (gupxz[i]*gxzy[i] + gupyz[i]*gyzy[i] + gupzz[i]*gzzy[i]) +

                gupxz[i] * (gupxz[i]*gxxz[i] + gupyz[i]*gxyz[i] + gupzz[i]*gxzz[i]) +
                gupyz[i] * (gupxz[i]*gxyz[i] + gupyz[i]*gyyz[i] + gupzz[i]*gyzz[i]) +
                gupzz[i] * (gupxz[i]*gxzz[i] + gupyz[i]*gyzz[i] + gupzz[i]*gzzz[i])
            );
        }
    }
    for (int i = 0; i < all; i += 1) {

        Gamxxx[i] = HALF * ( gupxx[i]*gxxx[i]
                        + gupxy[i]*(TWO*gxyx[i] - gxxy[i])
                        + gupxz[i]*(TWO*gxzx[i] - gxxz[i]) );

        Gamyxx[i] = HALF * ( gupxy[i]*gxxx[i]
                        + gupyy[i]*(TWO*gxyx[i] - gxxy[i])
                        + gupyz[i]*(TWO*gxzx[i] - gxxz[i]) );

        Gamzxx[i] = HALF * ( gupxz[i]*gxxx[i]
                        + gupyz[i]*(TWO*gxyx[i] - gxxy[i])
                        + gupzz[i]*(TWO*gxzx[i] - gxxz[i]) );

        Gamxyy[i] = HALF * ( gupxx[i]*(TWO*gxyy[i] - gyyx[i])
                        + gupxy[i]*gyyy[i]
                        + gupxz[i]*(TWO*gyzy[i] - gyyz[i]) );

        Gamyyy[i] = HALF * ( gupxy[i]*(TWO*gxyy[i] - gyyx[i])
                        + gupyy[i]*gyyy[i]
                        + gupyz[i]*(TWO*gyzy[i] - gyyz[i]) );

        Gamzyy[i] = HALF * ( gupxz[i]*(TWO*gxyy[i] - gyyx[i])
                        + gupyz[i]*gyyy[i]
                        + gupzz[i]*(TWO*gyzy[i] - gyyz[i]) );

        Gamxzz[i] = HALF * ( gupxx[i]*(TWO*gxzz[i] - gzzx[i])
                        + gupxy[i]*(TWO*gyzz[i] - gzzy[i])
                        + gupxz[i]*gzzz[i] );

        Gamyzz[i] = HALF * ( gupxy[i]*(TWO*gxzz[i] - gzzx[i])
                        + gupyy[i]*(TWO*gyzz[i] - gzzy[i])
                        + gupyz[i]*gzzz[i] );

        Gamzzz[i] = HALF * ( gupxz[i]*(TWO*gxzz[i] - gzzx[i])
                        + gupyz[i]*(TWO*gyzz[i] - gzzy[i])
                        + gupzz[i]*gzzz[i] );

        Gamxxy[i] = HALF * ( gupxx[i]*gxxy[i]
                        + gupxy[i]*gyyx[i]
                        + gupxz[i]*(gxzy[i] + gyzx[i] - gxyz[i]) );

        Gamyxy[i] = HALF * ( gupxy[i]*gxxy[i]
                        + gupyy[i]*gyyx[i]
                        + gupyz[i]*(gxzy[i] + gyzx[i] - gxyz[i]) );

        Gamzxy[i] = HALF * ( gupxz[i]*gxxy[i]
                        + gupyz[i]*gyyx[i]
                        + gupzz[i]*(gxzy[i] + gyzx[i] - gxyz[i]) );

        Gamxxz[i] = HALF * ( gupxx[i]*gxxz[i]
                        + gupxy[i]*(gxyz[i] + gyzx[i] - gxzy[i])
                        + gupxz[i]*gzzx[i] );

        Gamyxz[i] = HALF * ( gupxy[i]*gxxz[i]
                        + gupyy[i]*(gxyz[i] + gyzx[i] - gxzy[i])
                        + gupyz[i]*gzzx[i] );

        Gamzxz[i] = HALF * ( gupxz[i]*gxxz[i]
                        + gupyz[i]*(gxyz[i] + gyzx[i] - gxzy[i])
                        + gupzz[i]*gzzx[i] );

        Gamxyz[i] = HALF * ( gupxx[i]*(gxyz[i] + gxzy[i] - gyzx[i])
                        + gupxy[i]*gyyz[i]
                        + gupxz[i]*gzzy[i] );

        Gamyyz[i] = HALF * ( gupxy[i]*(gxyz[i] + gxzy[i] - gyzx[i])
                        + gupyy[i]*gyyz[i]
                        + gupyz[i]*gzzy[i] );

        Gamzyz[i] = HALF * ( gupxz[i]*(gxyz[i] + gxzy[i] - gyzx[i])
                        + gupyz[i]*gyyz[i]
                        + gupzz[i]*gzzy[i] );

    }
    for (int i = 0; i < all; i += 1) {

        Rxx[i] = gupxx[i]*gupxx[i]*Axx[i]
            + gupxy[i]*gupxy[i]*Ayy[i]
            + gupxz[i]*gupxz[i]*Azz[i]
            + TWO * ( gupxx[i]*gupxy[i]*Axy[i]
                    + gupxx[i]*gupxz[i]*Axz[i]
                    + gupxy[i]*gupxz[i]*Ayz[i] );

        Ryy[i] = gupxy[i]*gupxy[i]*Axx[i]
            + gupyy[i]*gupyy[i]*Ayy[i]
            + gupyz[i]*gupyz[i]*Azz[i]
            + TWO * ( gupxy[i]*gupyy[i]*Axy[i]
                    + gupxy[i]*gupyz[i]*Axz[i]
                    + gupyy[i]*gupyz[i]*Ayz[i] );

        Rzz[i] = gupxz[i]*gupxz[i]*Axx[i]
            + gupyz[i]*gupyz[i]*Ayy[i]
            + gupzz[i]*gupzz[i]*Azz[i]
            + TWO * ( gupxz[i]*gupyz[i]*Axy[i]
                    + gupxz[i]*gupzz[i]*Axz[i]
                    + gupyz[i]*gupzz[i]*Ayz[i] );

        Rxy[i] = gupxx[i]*gupxy[i]*Axx[i]
            + gupxy[i]*gupyy[i]*Ayy[i]
            + gupxz[i]*gupyz[i]*Azz[i]
            + ( gupxx[i]*gupyy[i] + gupxy[i]*gupxy[i] ) * Axy[i]
            + ( gupxx[i]*gupyz[i] + gupxz[i]*gupxy[i] ) * Axz[i]
            + ( gupxy[i]*gupyz[i] + gupxz[i]*gupyy[i] ) * Ayz[i];

        Rxz[i] = gupxx[i]*gupxz[i]*Axx[i]
            + gupxy[i]*gupyz[i]*Ayy[i]
            + gupxz[i]*gupzz[i]*Azz[i]
            + ( gupxx[i]*gupyz[i] + gupxy[i]*gupxz[i] ) * Axy[i]
            + ( gupxx[i]*gupzz[i] + gupxz[i]*gupxz[i] ) * Axz[i]
            + ( gupxy[i]*gupzz[i] + gupxz[i]*gupyz[i] ) * Ayz[i];

        Ryz[i] = gupxy[i]*gupxz[i]*Axx[i]
            + gupyy[i]*gupyz[i]*Ayy[i]
            + gupyz[i]*gupzz[i]*Azz[i]
            + ( gupxy[i]*gupyz[i] + gupyy[i]*gupxz[i] ) * Axy[i]
            + ( gupxy[i]*gupzz[i] + gupyz[i]*gupxz[i] ) * Axz[i]
            + ( gupyy[i]*gupzz[i] + gupyz[i]*gupyz[i] ) * Ayz[i];
    }
    fderivs(ex,Lap,Lapx,Lapy,Lapz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev);
    fderivs(ex,trK,Kx,Ky,Kz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev);

    for(int i=0;i<all;i+=1){
        Gamx_rhs[i] = - TWO * (   Lapx[i] * Rxx[i] +   Lapy[i] * Rxy[i] +   Lapz[i] * Rxz[i] ) + 
            TWO * alpn1[i] * (                                                
            -F3o2/chin1[i] * (   chix[i] * Rxx[i] +   chiy[i] * Rxy[i] +   chiz[i] * Rxz[i] ) - 
                gupxx[i] * (   F2o3 * Kx[i]  +  EIGHT * PI * Sx[i]            ) - 
                gupxy[i] * (   F2o3 * Ky[i]  +  EIGHT * PI * Sy[i]            ) - 
                gupxz[i] * (   F2o3 * Kz[i]  +  EIGHT * PI * Sz[i]            ) + 
                            Gamxxx[i] * Rxx[i] + Gamxyy[i] * Ryy[i] + Gamxzz[i] * Rzz[i]   + 
                    TWO * ( Gamxxy[i] * Rxy[i] + Gamxxz[i] * Rxz[i] + Gamxyz[i] * Ryz[i] ) );

        Gamy_rhs[i] = -TWO * ( Lapx[i]*Rxy[i] + Lapy[i]*Ryy[i] + Lapz[i]*Ryz[i] )
                    +  TWO * alpn1[i] * (
                        -F3o2/chin1[i] * ( chix[i]*Rxy[i] + chiy[i]*Ryy[i] + chiz[i]*Ryz[i] )
                        - gupxy[i] * ( F2o3*Kx[i] + EIGHT*PI*Sx[i] )
                        - gupyy[i] * ( F2o3*Ky[i] + EIGHT*PI*Sy[i] )
                        - gupyz[i] * ( F2o3*Kz[i] + EIGHT*PI*Sz[i] )
                        + Gamyxx[i]*Rxx[i] + Gamyyy[i]*Ryy[i] + Gamyzz[i]*Rzz[i]
                        + TWO * ( Gamyxy[i]*Rxy[i] + Gamyxz[i]*Rxz[i] + Gamyyz[i]*Ryz[i] )
                    );

        Gamz_rhs[i] = -TWO * ( Lapx[i]*Rxz[i] + Lapy[i]*Ryz[i] + Lapz[i]*Rzz[i] )
                    +  TWO * alpn1[i] * (
                        -F3o2/chin1[i] * ( chix[i]*Rxz[i] + chiy[i]*Ryz[i] + chiz[i]*Rzz[i] )
                        - gupxz[i] * ( F2o3*Kx[i] + EIGHT*PI*Sx[i] )
                        - gupyz[i] * ( F2o3*Ky[i] + EIGHT*PI*Sy[i] )
                        - gupzz[i] * ( F2o3*Kz[i] + EIGHT*PI*Sz[i] )
                        + Gamzxx[i]*Rxx[i] + Gamzyy[i]*Ryy[i] + Gamzzz[i]*Rzz[i]
                        + TWO * ( Gamzxy[i]*Rxy[i] + Gamzxz[i]*Rxz[i] + Gamzyz[i]*Ryz[i] )
                    );
    }

    fdderivs(ex,betax,gxxx,gxyx,gxzx,gyyx,gyzx,gzzx,
                X,Y,Z,ANTI,SYM, SYM ,Symmetry,Lev);
    fdderivs(ex,betay,gxxy,gxyy,gxzy,gyyy,gyzy,gzzy,
                X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev);
    fdderivs(ex,betaz,gxxz,gxyz,gxzz,gyyz,gyzz,gzzz,
                X,Y,Z,SYM ,SYM, ANTI,Symmetry,Lev);
    
    for(int i=0;i<all;i+=1){
        fxx[i] = gxxx[i] + gxyy[i] + gxzz[i];
        fxy[i] = gxyx[i] + gyyy[i] + gyzz[i];
        fxz[i] = gxzx[i] + gyzy[i] + gzzz[i];
    }

    for(int i=0;i<all;i+=1){
        Gamxa[i] = gupxx[i]*Gamxxx[i] + gupyy[i]*Gamxyy[i] + gupzz[i]*Gamxzz[i]
                + TWO * ( gupxy[i]*Gamxxy[i] + gupxz[i]*Gamxxz[i] + gupyz[i]*Gamxyz[i] );

        Gamya[i] = gupxx[i]*Gamyxx[i] + gupyy[i]*Gamyyy[i] + gupzz[i]*Gamyzz[i]
                + TWO * ( gupxy[i]*Gamyxy[i] + gupxz[i]*Gamyxz[i] + gupyz[i]*Gamyyz[i] );

        Gamza[i] = gupxx[i]*Gamzxx[i] + gupyy[i]*Gamzyy[i] + gupzz[i]*Gamzzz[i]
                + TWO * ( gupxy[i]*Gamzxy[i] + gupxz[i]*Gamzxz[i] + gupyz[i]*Gamzyz[i] );
    }
    fderivs(ex,Gamx,Gamxx,Gamxy,Gamxz,X,Y,Z,ANTI,SYM ,SYM ,Symmetry,Lev);
    fderivs(ex,Gamy,Gamyx,Gamyy,Gamyz,X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev);
    fderivs(ex,Gamz,Gamzx,Gamzy,Gamzz,X,Y,Z,SYM ,SYM ,ANTI,Symmetry,Lev);
    for(int i=0;i<all;i+=1){
        Gamx_rhs[i] = Gamx_rhs[i]
                    + F2o3 * Gamxa[i] * div_beta[i]
                    - Gamxa[i] * betaxx[i] - Gamya[i] * betaxy[i] - Gamza[i] * betaxz[i]
                    + F1o3 * ( gupxx[i] * fxx[i] + gupxy[i] * fxy[i] + gupxz[i] * fxz[i] )
                    + gupxx[i] * gxxx[i] + gupyy[i] * gyyx[i] + gupzz[i] * gzzx[i]
                    + TWO * ( gupxy[i] * gxyx[i] + gupxz[i] * gxzx[i] + gupyz[i] * gyzx[i] );

        Gamy_rhs[i] = Gamy_rhs[i]
                    + F2o3 * Gamya[i] * div_beta[i]
                    - Gamxa[i] * betayx[i] - Gamya[i] * betayy[i] - Gamza[i] * betayz[i]
                    + F1o3 * ( gupxy[i] * fxx[i] + gupyy[i] * fxy[i] + gupyz[i] * fxz[i] )
                    + gupxx[i] * gxxy[i] + gupyy[i] * gyyy[i] + gupzz[i] * gzzy[i]
                    + TWO * ( gupxy[i] * gxyy[i] + gupxz[i] * gxzy[i] + gupyz[i] * gyzy[i] );

        Gamz_rhs[i] = Gamz_rhs[i]
                    + F2o3 * Gamza[i] * div_beta[i]
                    - Gamxa[i] * betazx[i] - Gamya[i] * betazy[i] - Gamza[i] * betazz[i]
                    + F1o3 * ( gupxz[i] * fxx[i] + gupyz[i] * fxy[i] + gupzz[i] * fxz[i] )
                    + gupxx[i] * gxxz[i] + gupyy[i] * gyyz[i] + gupzz[i] * gzzz[i]
                    + TWO * ( gupxy[i] * gxyz[i] + gupxz[i] * gxzz[i] + gupyz[i] * gyzz[i] );
    }
    for (int i = 0; i < all; i += 1) {

        gxxx[i] = gxx[i]*Gamxxx[i] + gxy[i]*Gamyxx[i] + gxz[i]*Gamzxx[i];
        gxyx[i] = gxx[i]*Gamxxy[i] + gxy[i]*Gamyxy[i] + gxz[i]*Gamzxy[i];
        gxzx[i] = gxx[i]*Gamxxz[i] + gxy[i]*Gamyxz[i] + gxz[i]*Gamzxz[i];
        gyyx[i] = gxx[i]*Gamxyy[i] + gxy[i]*Gamyyy[i] + gxz[i]*Gamzyy[i];
        gyzx[i] = gxx[i]*Gamxyz[i] + gxy[i]*Gamyyz[i] + gxz[i]*Gamzyz[i];
        gzzx[i] = gxx[i]*Gamxzz[i] + gxy[i]*Gamyzz[i] + gxz[i]*Gamzzz[i];

        gxxy[i] = gxy[i]*Gamxxx[i] + gyy[i]*Gamyxx[i] + gyz[i]*Gamzxx[i];
        gxyy[i] = gxy[i]*Gamxxy[i] + gyy[i]*Gamyxy[i] + gyz[i]*Gamzxy[i];
        gxzy[i] = gxy[i]*Gamxxz[i] + gyy[i]*Gamyxz[i] + gyz[i]*Gamzxz[i];
        gyyy[i] = gxy[i]*Gamxyy[i] + gyy[i]*Gamyyy[i] + gyz[i]*Gamzyy[i];
        gyzy[i] = gxy[i]*Gamxyz[i] + gyy[i]*Gamyyz[i] + gyz[i]*Gamzyz[i];
        gzzy[i] = gxy[i]*Gamxzz[i] + gyy[i]*Gamyzz[i] + gyz[i]*Gamzzz[i];

        gxxz[i] = gxz[i]*Gamxxx[i] + gyz[i]*Gamyxx[i] + gzz[i]*Gamzxx[i];
        gxyz[i] = gxz[i]*Gamxxy[i] + gyz[i]*Gamyxy[i] + gzz[i]*Gamzxy[i];
        gxzz[i] = gxz[i]*Gamxxz[i] + gyz[i]*Gamyxz[i] + gzz[i]*Gamzxz[i];
        gyyz[i] = gxz[i]*Gamxyy[i] + gyz[i]*Gamyyy[i] + gzz[i]*Gamzyy[i];
        gyzz[i] = gxz[i]*Gamxyz[i] + gyz[i]*Gamyyz[i] + gzz[i]*Gamzyz[i];
        gzzz[i] = gxz[i]*Gamxzz[i] + gyz[i]*Gamyzz[i] + gzz[i]*Gamzzz[i];
    }


    fdderivs(ex,dxx,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev);
    for (int i = 0; i < all; i += 1) {
        Rxx[i] =  gupxx[i] * fxx[i] + gupyy[i] * fyy[i] + gupzz[i] * fzz[i]
                + (gupxy[i] * fxy[i] + gupxz[i] * fxz[i] + gupyz[i] * fyz[i]) * TWO;
    }

    fdderivs(ex,dzz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev);
    for (int i = 0; i < all; i += 1) {
        Ryy[i] =  gupxx[i] * fxx[i] + gupyy[i] * fyy[i] + gupzz[i] * fzz[i]
                + (gupxy[i] * fxy[i] + gupxz[i] * fxz[i] + gupyz[i] * fyz[i]) * TWO;
    }

    fdderivs(ex,dyy,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev);
    for (int i = 0; i < all; i += 1) {
        Rzz[i] =  gupxx[i] * fxx[i] + gupyy[i] * fyy[i] + gupzz[i] * fzz[i]
                + (gupxy[i] * fxy[i] + gupxz[i] * fxz[i] + gupyz[i] * fyz[i]) * TWO;
    }

    fdderivs(ex,gxy,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,ANTI, ANTI,SYM ,Symmetry,Lev);
    for (int i = 0; i < all; i += 1) {
        Rxy[i] =  gupxx[i] * fxx[i] + gupyy[i] * fyy[i] + gupzz[i] * fzz[i]
                + (gupxy[i] * fxy[i] + gupxz[i] * fxz[i] + gupyz[i] * fyz[i]) * TWO;
    }

    fdderivs(ex,gxz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,ANTI ,SYM ,ANTI,Symmetry,Lev);
    for (int i = 0; i < all; i += 1) {
        Rxz[i] =  gupxx[i] * fxx[i] + gupyy[i] * fyy[i] + gupzz[i] * fzz[i]
                + (gupxy[i] * fxy[i] + gupxz[i] * fxz[i] + gupyz[i] * fyz[i]) * TWO;
    }

    fdderivs(ex,gyz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,ANTI ,ANTI,Symmetry,Lev);
    for (int i = 0; i < all; i += 1) {
        Ryz[i] =  gupxx[i] * fxx[i] + gupyy[i] * fyy[i] + gupzz[i] * fzz[i]
                + (gupxy[i] * fxy[i] + gupxz[i] * fxz[i] + gupyz[i] * fyz[i]) * TWO;
    }
/* 假设 all = ex1*ex2*ex3，所有量都是 length=all 的 double 数组（已按同一扁平化规则排布） */

    for (int i = 0; i < all; i += 1) {
        Rxx[i] =
            -HALF * Rxx[i]
            + gxx[i] * Gamxx[i] + gxy[i] * Gamyx[i] + gxz[i] * Gamzx[i]
            + Gamxa[i] * gxxx[i] + Gamya[i] * gxyx[i] + Gamza[i] * gxzx[i]
            + gupxx[i] * (
                TWO * (Gamxxx[i] * gxxx[i] + Gamyxx[i] * gxyx[i] + Gamzxx[i] * gxzx[i]) +
                    (Gamxxx[i] * gxxx[i] + Gamyxx[i] * gxxy[i] + Gamzxx[i] * gxxz[i])
            )
            + gupxy[i] * (
                TWO * (Gamxxx[i] * gxyx[i] + Gamyxx[i] * gyyx[i] + Gamzxx[i] * gyzx[i] +
                    Gamxxy[i] * gxxx[i] + Gamyxy[i] * gxyx[i] + Gamzxy[i] * gxzx[i]) +
                    (Gamxxy[i] * gxxx[i] + Gamyxy[i] * gxxy[i] + Gamzxy[i] * gxxz[i]) +
                    (Gamxxx[i] * gxyx[i] + Gamyxx[i] * gxyy[i] + Gamzxx[i] * gxyz[i])
            )
            + gupxz[i] * (
                TWO * (Gamxxx[i] * gxzx[i] + Gamyxx[i] * gyzx[i] + Gamzxx[i] * gzzx[i] +
                    Gamxxz[i] * gxxx[i] + Gamyxz[i] * gxyx[i] + Gamzxz[i] * gxzx[i]) +
                    (Gamxxz[i] * gxxx[i] + Gamyxz[i] * gxxy[i] + Gamzxz[i] * gxxz[i]) +
                    (Gamxxx[i] * gxzx[i] + Gamyxx[i] * gxzy[i] + Gamzxx[i] * gxzz[i])
            )
            + gupyy[i] * (
                TWO * (Gamxxy[i] * gxyx[i] + Gamyxy[i] * gyyx[i] + Gamzxy[i] * gyzx[i]) +
                    (Gamxxy[i] * gxyx[i] + Gamyxy[i] * gxyy[i] + Gamzxy[i] * gxyz[i])
            )
            + gupyz[i] * (
                TWO * (Gamxxy[i] * gxzx[i] + Gamyxy[i] * gyzx[i] + Gamzxy[i] * gzzx[i] +
                    Gamxxz[i] * gxyx[i] + Gamyxz[i] * gyyx[i] + Gamzxz[i] * gyzx[i]) +
                    (Gamxxz[i] * gxyx[i] + Gamyxz[i] * gxyy[i] + Gamzxz[i] * gxyz[i]) +
                    (Gamxxy[i] * gxzx[i] + Gamyxy[i] * gxzy[i] + Gamzxy[i] * gxzz[i])
            )
            + gupzz[i] * (
                TWO * (Gamxxz[i] * gxzx[i] + Gamyxz[i] * gyzx[i] + Gamzxz[i] * gzzx[i]) +
                    (Gamxxz[i] * gxzx[i] + Gamyxz[i] * gxzy[i] + Gamzxz[i] * gxzz[i])
            );

        Ryy[i] =
            -HALF * Ryy[i]
            + gxy[i] * Gamxy[i] + gyy[i] * Gamyy[i] + gyz[i] * Gamzy[i]
            + Gamxa[i] * gxyy[i] + Gamya[i] * gyyy[i] + Gamza[i] * gyzy[i]
            + gupxx[i] * (
                TWO * (Gamxxy[i] * gxxy[i] + Gamyxy[i] * gxyy[i] + Gamzxy[i] * gxzy[i]) +
                    (Gamxxy[i] * gxyx[i] + Gamyxy[i] * gxyy[i] + Gamzxy[i] * gxyz[i])
            )
            + gupxy[i] * (
                TWO * (Gamxxy[i] * gxyy[i] + Gamyxy[i] * gyyy[i] + Gamzxy[i] * gyzy[i] +
                    Gamxyy[i] * gxxy[i] + Gamyyy[i] * gxyy[i] + Gamzyy[i] * gxzy[i]) +
                    (Gamxyy[i] * gxyx[i] + Gamyyy[i] * gxyy[i] + Gamzyy[i] * gxyz[i]) +
                    (Gamxxy[i] * gyyx[i] + Gamyxy[i] * gyyy[i] + Gamzxy[i] * gyyz[i])
            )
            + gupxz[i] * (
                TWO * (Gamxxy[i] * gxzy[i] + Gamyxy[i] * gyzy[i] + Gamzxy[i] * gzzy[i] +
                    Gamxyz[i] * gxxy[i] + Gamyyz[i] * gxyy[i] + Gamzyz[i] * gxzy[i]) +
                    (Gamxyz[i] * gxyx[i] + Gamyyz[i] * gxyy[i] + Gamzyz[i] * gxyz[i]) +
                    (Gamxxy[i] * gyzx[i] + Gamyxy[i] * gyzy[i] + Gamzxy[i] * gyzz[i])
            )
            + gupyy[i] * (
                TWO * (Gamxyy[i] * gxyy[i] + Gamyyy[i] * gyyy[i] + Gamzyy[i] * gyzy[i]) +
                    (Gamxyy[i] * gyyx[i] + Gamyyy[i] * gyyy[i] + Gamzyy[i] * gyyz[i])
            )
            + gupyz[i] * (
                TWO * (Gamxyy[i] * gxzy[i] + Gamyyy[i] * gyzy[i] + Gamzyy[i] * gzzy[i] +
                    Gamxyz[i] * gxyy[i] + Gamyyz[i] * gyyy[i] + Gamzyz[i] * gyzy[i]) +
                    (Gamxyz[i] * gyyx[i] + Gamyyz[i] * gyyy[i] + Gamzyz[i] * gyyz[i]) +
                    (Gamxyy[i] * gyzx[i] + Gamyyy[i] * gyzy[i] + Gamzyy[i] * gyzz[i])
            )
            + gupzz[i] * (
                TWO * (Gamxyz[i] * gxzy[i] + Gamyyz[i] * gyzy[i] + Gamzyz[i] * gzzy[i]) +
                    (Gamxyz[i] * gyzx[i] + Gamyyz[i] * gyzy[i] + Gamzyz[i] * gyzz[i])
            );

        Rzz[i] =
            -HALF * Rzz[i]
            + gxz[i] * Gamxz[i] + gyz[i] * Gamyz[i] + gzz[i] * Gamzz[i]
            + Gamxa[i] * gxzz[i] + Gamya[i] * gyzz[i] + Gamza[i] * gzzz[i]
            + gupxx[i] * (
                TWO * (Gamxxz[i] * gxxz[i] + Gamyxz[i] * gxyz[i] + Gamzxz[i] * gxzz[i]) +
                    (Gamxxz[i] * gxzx[i] + Gamyxz[i] * gxzy[i] + Gamzxz[i] * gxzz[i])
            )
            + gupxy[i] * (
                TWO * (Gamxxz[i] * gxyz[i] + Gamyxz[i] * gyyz[i] + Gamzxz[i] * gyzz[i] +
                    Gamxyz[i] * gxxz[i] + Gamyyz[i] * gxyz[i] + Gamzyz[i] * gxzz[i]) +
                    (Gamxyz[i] * gxzx[i] + Gamyyz[i] * gxzy[i] + Gamzyz[i] * gxzz[i]) +
                    (Gamxxz[i] * gyzx[i] + Gamyxz[i] * gyzy[i] + Gamzxz[i] * gyzz[i])
            )
            + gupxz[i] * (
                TWO * (Gamxxz[i] * gxzz[i] + Gamyxz[i] * gyzz[i] + Gamzxz[i] * gzzz[i] +
                    Gamxzz[i] * gxxz[i] + Gamyzz[i] * gxyz[i] + Gamzzz[i] * gxzz[i]) +
                    (Gamxzz[i] * gxzx[i] + Gamyzz[i] * gxzy[i] + Gamzzz[i] * gxzz[i]) +
                    (Gamxxz[i] * gzzx[i] + Gamyxz[i] * gzzy[i] + Gamzxz[i] * gzzz[i])
            )
            + gupyy[i] * (
                TWO * (Gamxyz[i] * gxyz[i] + Gamyyz[i] * gyyz[i] + Gamzyz[i] * gyzz[i]) +
                    (Gamxyz[i] * gyzx[i] + Gamyyz[i] * gyzy[i] + Gamzyz[i] * gyzz[i])
            )
            + gupyz[i] * (
                TWO * (Gamxyz[i] * gxzz[i] + Gamyyz[i] * gyzz[i] + Gamzyz[i] * gzzz[i] +
                    Gamxzz[i] * gxyz[i] + Gamyzz[i] * gyyz[i] + Gamzzz[i] * gyzz[i]) +
                    (Gamxzz[i] * gyzx[i] + Gamyzz[i] * gyzy[i] + Gamzzz[i] * gyzz[i]) +
                    (Gamxyz[i] * gzzx[i] + Gamyyz[i] * gzzy[i] + Gamzyz[i] * gzzz[i])
            )
            + gupzz[i] * (
                TWO * (Gamxzz[i] * gxzz[i] + Gamyzz[i] * gyzz[i] + Gamzzz[i] * gzzz[i]) +
                    (Gamxzz[i] * gzzx[i] + Gamyzz[i] * gzzy[i] + Gamzzz[i] * gzzz[i])
            );

        Rxy[i] =
            HALF * (
                -Rxy[i]
                + gxx[i] * Gamxy[i] + gxy[i] * Gamyy[i] + gxz[i] * Gamzy[i]
                + gxy[i] * Gamxx[i] + gyy[i] * Gamyx[i] + gyz[i] * Gamzx[i]
                + Gamxa[i] * gxyx[i] + Gamya[i] * gyyx[i] + Gamza[i] * gyzx[i]
                + Gamxa[i] * gxxy[i] + Gamya[i] * gxyy[i] + Gamza[i] * gxzy[i]
            )
            + gupxx[i] * (
                Gamxxx[i] * gxxy[i] + Gamyxx[i] * gxyy[i] + Gamzxx[i] * gxzy[i]
                + Gamxxy[i] * gxxx[i] + Gamyxy[i] * gxyx[i] + Gamzxy[i] * gxzx[i]
                + Gamxxx[i] * gxyx[i] + Gamyxx[i] * gxyy[i] + Gamzxx[i] * gxyz[i]
            )
            + gupxy[i] * (
                Gamxxx[i] * gxyy[i] + Gamyxx[i] * gyyy[i] + Gamzxx[i] * gyzy[i]
                + Gamxxy[i] * gxyx[i] + Gamyxy[i] * gyyx[i] + Gamzxy[i] * gyzx[i]
                + Gamxxy[i] * gxyx[i] + Gamyxy[i] * gxyy[i] + Gamzxy[i] * gxyz[i]
                + Gamxxy[i] * gxxy[i] + Gamyxy[i] * gxyy[i] + Gamzxy[i] * gxzy[i]
                + Gamxyy[i] * gxxx[i] + Gamyyy[i] * gxyx[i] + Gamzyy[i] * gxzx[i]
                + Gamxxx[i] * gyyx[i] + Gamyxx[i] * gyyy[i] + Gamzxx[i] * gyyz[i]
            )
            + gupxz[i] * (
                Gamxxx[i] * gxzy[i] + Gamyxx[i] * gyzy[i] + Gamzxx[i] * gzzy[i]
                + Gamxxy[i] * gxzx[i] + Gamyxy[i] * gyzx[i] + Gamzxy[i] * gzzx[i]
                + Gamxxz[i] * gxyx[i] + Gamyxz[i] * gxyy[i] + Gamzxz[i] * gxyz[i]
                + Gamxxz[i] * gxxy[i] + Gamyxz[i] * gxyy[i] + Gamzxz[i] * gxzy[i]
                + Gamxyz[i] * gxxx[i] + Gamyyz[i] * gxyx[i] + Gamzyz[i] * gxzx[i]
                + Gamxxx[i] * gyzx[i] + Gamyxx[i] * gyzy[i] + Gamzxx[i] * gyzz[i]
            )
            + gupyy[i] * (
                Gamxxy[i] * gxyy[i] + Gamyxy[i] * gyyy[i] + Gamzxy[i] * gyzy[i]
                + Gamxyy[i] * gxyx[i] + Gamyyy[i] * gyyx[i] + Gamzyy[i] * gyzx[i]
                + Gamxxy[i] * gyyx[i] + Gamyxy[i] * gyyy[i] + Gamzxy[i] * gyyz[i]
            )
            + gupyz[i] * (
                Gamxxy[i] * gxzy[i] + Gamyxy[i] * gyzy[i] + Gamzxy[i] * gzzy[i]
                + Gamxyy[i] * gxzx[i] + Gamyyy[i] * gyzx[i] + Gamzyy[i] * gzzx[i]
                + Gamxxz[i] * gyyx[i] + Gamyxz[i] * gyyy[i] + Gamzxz[i] * gyyz[i]
                + Gamxxz[i] * gxyy[i] + Gamyxz[i] * gyyy[i] + Gamzxz[i] * gyzy[i]
                + Gamxyz[i] * gxyx[i] + Gamyyz[i] * gyyx[i] + Gamzyz[i] * gyzx[i]
                + Gamxxy[i] * gyzx[i] + Gamyxy[i] * gyzy[i] + Gamzxy[i] * gyzz[i]
            )
            + gupzz[i] * (
                Gamxxz[i] * gxzy[i] + Gamyxz[i] * gyzy[i] + Gamzxz[i] * gzzy[i]
                + Gamxyz[i] * gxzx[i] + Gamyyz[i] * gyzx[i] + Gamzyz[i] * gzzx[i]
                + Gamxxz[i] * gyzx[i] + Gamyxz[i] * gyzy[i] + Gamzxz[i] * gyzz[i]
            );

        Rxz[i] =
            HALF * (
                -Rxz[i]
                + gxx[i] * Gamxz[i] + gxy[i] * Gamyz[i] + gxz[i] * Gamzz[i]
                + gxz[i] * Gamxx[i] + gyz[i] * Gamyx[i] + gzz[i] * Gamzx[i]
                + Gamxa[i] * gxzx[i] + Gamya[i] * gyzx[i] + Gamza[i] * gzzx[i]
                + Gamxa[i] * gxxz[i] + Gamya[i] * gxyz[i] + Gamza[i] * gxzz[i]
            )
            + gupxx[i] * (
                Gamxxx[i] * gxxz[i] + Gamyxx[i] * gxyz[i] + Gamzxx[i] * gxzz[i]
                + Gamxxz[i] * gxxx[i] + Gamyxz[i] * gxyx[i] + Gamzxz[i] * gxzx[i]
                + Gamxxx[i] * gxzx[i] + Gamyxx[i] * gxzy[i] + Gamzxx[i] * gxzz[i]
            )
            + gupxy[i] * (
                Gamxxx[i] * gxyz[i] + Gamyxx[i] * gyyz[i] + Gamzxx[i] * gyzz[i]
                + Gamxxz[i] * gxyx[i] + Gamyxz[i] * gyyx[i] + Gamzxz[i] * gyzx[i]
                + Gamxxy[i] * gxzx[i] + Gamyxy[i] * gxzy[i] + Gamzxy[i] * gxzz[i]
                + Gamxxy[i] * gxxz[i] + Gamyxy[i] * gxyz[i] + Gamzxy[i] * gxzz[i]
                + Gamxyz[i] * gxxx[i] + Gamyyz[i] * gxyx[i] + Gamzyz[i] * gxzx[i]
                + Gamxxx[i] * gyzx[i] + Gamyxx[i] * gyzy[i] + Gamzxx[i] * gyzz[i]
            )
            + gupxz[i] * (
                Gamxxx[i] * gxzz[i] + Gamyxx[i] * gyzz[i] + Gamzxx[i] * gzzz[i]
                + Gamxxz[i] * gxzx[i] + Gamyxz[i] * gyzx[i] + Gamzxz[i] * gzzx[i]
                + Gamxxz[i] * gxzx[i] + Gamyxz[i] * gxzy[i] + Gamzxz[i] * gxzz[i]
                + Gamxxz[i] * gxxz[i] + Gamyxz[i] * gxyz[i] + Gamzxz[i] * gxzz[i]
                + Gamxzz[i] * gxxx[i] + Gamyzz[i] * gxyx[i] + Gamzzz[i] * gxzx[i]
                + Gamxxx[i] * gzzx[i] + Gamyxx[i] * gzzy[i] + Gamzxx[i] * gzzz[i]
            )
            + gupyy[i] * (
                Gamxxy[i] * gxyz[i] + Gamyxy[i] * gyyz[i] + Gamzxy[i] * gyzz[i]
                + Gamxyz[i] * gxyx[i] + Gamyyz[i] * gyyx[i] + Gamzyz[i] * gyzx[i]
                + Gamxxy[i] * gyzx[i] + Gamyxy[i] * gyzy[i] + Gamzxy[i] * gyzz[i]
            )
            + gupyz[i] * (
                Gamxxy[i] * gxzz[i] + Gamyxy[i] * gyzz[i] + Gamzxy[i] * gzzz[i]
                + Gamxyz[i] * gxzx[i] + Gamyyz[i] * gyzx[i] + Gamzyz[i] * gzzx[i]
                + Gamxxz[i] * gyzx[i] + Gamyxz[i] * gyzy[i] + Gamzxz[i] * gyzz[i]
                + Gamxxz[i] * gxyz[i] + Gamyxz[i] * gyyz[i] + Gamzxz[i] * gyzz[i]
                + Gamxzz[i] * gxyx[i] + Gamyzz[i] * gyyx[i] + Gamzzz[i] * gyzx[i]
                + Gamxxy[i] * gzzx[i] + Gamyxy[i] * gzzy[i] + Gamzxy[i] * gzzz[i]
            )
            + gupzz[i] * (
                Gamxxz[i] * gxzz[i] + Gamyxz[i] * gyzz[i] + Gamzxz[i] * gzzz[i]
                + Gamxzz[i] * gxzx[i] + Gamyzz[i] * gyzx[i] + Gamzzz[i] * gzzx[i]
                + Gamxxz[i] * gzzx[i] + Gamyxz[i] * gzzy[i] + Gamzxz[i] * gzzz[i]
            );

        Ryz[i] =
            HALF * (
                -Ryz[i]
                + gxy[i] * Gamxz[i] + gyy[i] * Gamyz[i] + gyz[i] * Gamzz[i]
                + gxz[i] * Gamxy[i] + gyz[i] * Gamyy[i] + gzz[i] * Gamzy[i]
                + Gamxa[i] * gxzy[i] + Gamya[i] * gyzy[i] + Gamza[i] * gzzy[i]
                + Gamxa[i] * gxyz[i] + Gamya[i] * gyyz[i] + Gamza[i] * gyzz[i]
            )
            + gupxx[i] * (
                Gamxxy[i] * gxxz[i] + Gamyxy[i] * gxyz[i] + Gamzxy[i] * gxzz[i]
                + Gamxxz[i] * gxxy[i] + Gamyxz[i] * gxyy[i] + Gamzxz[i] * gxzy[i]
                + Gamxxy[i] * gxzx[i] + Gamyxy[i] * gxzy[i] + Gamzxy[i] * gxzz[i]
            )
            + gupxy[i] * (
                Gamxxy[i] * gxyz[i] + Gamyxy[i] * gyyz[i] + Gamzxy[i] * gyzz[i]
                + Gamxxz[i] * gxyy[i] + Gamyxz[i] * gyyy[i] + Gamzxz[i] * gyzy[i]
                + Gamxyy[i] * gxzx[i] + Gamyyy[i] * gxzy[i] + Gamzyy[i] * gxzz[i]
                + Gamxyy[i] * gxxz[i] + Gamyyy[i] * gxyz[i] + Gamzyy[i] * gxzz[i]
                + Gamxyz[i] * gxxy[i] + Gamyyz[i] * gxyy[i] + Gamzyz[i] * gxzy[i]
                + Gamxxy[i] * gyzx[i] + Gamyxy[i] * gyzy[i] + Gamzxy[i] * gyzz[i]
            )
            + gupxz[i] * (
                Gamxxy[i] * gxzz[i] + Gamyxy[i] * gyzz[i] + Gamzxy[i] * gzzz[i]
                + Gamxxz[i] * gxzy[i] + Gamyxz[i] * gyzy[i] + Gamzxz[i] * gzzy[i]
                + Gamxyz[i] * gxzx[i] + Gamyyz[i] * gxzy[i] + Gamzyz[i] * gxzz[i]
                + Gamxyz[i] * gxxz[i] + Gamyyz[i] * gxyz[i] + Gamzyz[i] * gxzz[i]
                + Gamxzz[i] * gxxy[i] + Gamyzz[i] * gxyy[i] + Gamzzz[i] * gxzy[i]
                + Gamxxy[i] * gzzx[i] + Gamyxy[i] * gzzy[i] + Gamzxy[i] * gzzz[i]
            )
            + gupyy[i] * (
                Gamxyy[i] * gxyz[i] + Gamyyy[i] * gyyz[i] + Gamzyy[i] * gyzz[i]
                + Gamxyz[i] * gxyy[i] + Gamyyz[i] * gyyy[i] + Gamzyz[i] * gyzy[i]
                + Gamxyy[i] * gyzx[i] + Gamyyy[i] * gyzy[i] + Gamzyy[i] * gyzz[i]
            )
            + gupyz[i] * (
                Gamxyy[i] * gxzz[i] + Gamyyy[i] * gyzz[i] + Gamzyy[i] * gzzz[i]
                + Gamxyz[i] * gxzy[i] + Gamyyz[i] * gyzy[i] + Gamzyz[i] * gzzy[i]
                + Gamxyz[i] * gyzx[i] + Gamyyz[i] * gyzy[i] + Gamzyz[i] * gyzz[i]
                + Gamxyz[i] * gxyz[i] + Gamyyz[i] * gyyz[i] + Gamzyz[i] * gyzz[i]
                + Gamxzz[i] * gxyy[i] + Gamyzz[i] * gyyy[i] + Gamzzz[i] * gyzy[i]
                + Gamxyy[i] * gzzx[i] + Gamyyy[i] * gzzy[i] + Gamzyy[i] * gzzz[i]
            )
            + gupzz[i] * (
                Gamxyz[i] * gxzz[i] + Gamyyz[i] * gyzz[i] + Gamzyz[i] * gzzz[i]
                + Gamxzz[i] * gxzy[i] + Gamyzz[i] * gyzy[i] + Gamzzz[i] * gzzy[i]
                + Gamxyz[i] * gzzx[i] + Gamyyz[i] * gzzy[i] + Gamzyz[i] * gzzz[i]
            );
    }
    fdderivs(ex,chi,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev);
    for (int i = 0; i < all; i += 1) {
        fxx[i] = fxx[i] - Gamxxx[i] * chix[i] - Gamyxx[i] * chiy[i] - Gamzxx[i] * chiz[i];
        fxy[i] = fxy[i] - Gamxxy[i] * chix[i] - Gamyxy[i] * chiy[i] - Gamzxy[i] * chiz[i];
        fxz[i] = fxz[i] - Gamxxz[i] * chix[i] - Gamyxz[i] * chiy[i] - Gamzxz[i] * chiz[i];
        fyy[i] = fyy[i] - Gamxyy[i] * chix[i] - Gamyyy[i] * chiy[i] - Gamzyy[i] * chiz[i];
        fyz[i] = fyz[i] - Gamxyz[i] * chix[i] - Gamyyz[i] * chiy[i] - Gamzyz[i] * chiz[i];
        fzz[i] = fzz[i] - Gamxzz[i] * chix[i] - Gamyzz[i] * chiy[i] - Gamzzz[i] * chiz[i];
    }
    for (int i = 0; i < all; i += 1) {
        f[i] =
            gupxx[i] * (fxx[i] - (F3o2 / chin1[i]) * chix[i] * chix[i])
            + gupyy[i] * (fyy[i] - (F3o2 / chin1[i]) * chiy[i] * chiy[i])
            + gupzz[i] * (fzz[i] - (F3o2 / chin1[i]) * chiz[i] * chiz[i])
            + TWO * gupxy[i] * (fxy[i] - (F3o2 / chin1[i]) * chix[i] * chiy[i])
            + TWO * gupxz[i] * (fxz[i] - (F3o2 / chin1[i]) * chix[i] * chiz[i])
            + TWO * gupyz[i] * (fyz[i] - (F3o2 / chin1[i]) * chiy[i] * chiz[i]);
    }

    for (int i = 0; i < all; i += 1) {
        Rxx[i] = Rxx[i] + ( fxx[i] - (chix[i] * chix[i]) / (chin1[i] * TWO) + gxx[i] * f[i] ) / (chin1[i] * TWO);
        Ryy[i] = Ryy[i] + ( fyy[i] - (chiy[i] * chiy[i]) / (chin1[i] * TWO) + gyy[i] * f[i] ) / (chin1[i] * TWO);
        Rzz[i] = Rzz[i] + ( fzz[i] - (chiz[i] * chiz[i]) / (chin1[i] * TWO) + gzz[i] * f[i] ) / (chin1[i] * TWO);

        Rxy[i] = Rxy[i] + ( fxy[i] - (chix[i] * chiy[i]) / (chin1[i] * TWO) + gxy[i] * f[i] ) / (chin1[i] * TWO);
        Rxz[i] = Rxz[i] + ( fxz[i] - (chix[i] * chiz[i]) / (chin1[i] * TWO) + gxz[i] * f[i] ) / (chin1[i] * TWO);
        Ryz[i] = Ryz[i] + ( fyz[i] - (chiy[i] * chiz[i]) / (chin1[i] * TWO) + gyz[i] * f[i] ) / (chin1[i] * TWO);
    }

    fdderivs(ex,Lap,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev);

    for (int i = 0; i < all; i += 1) {
        /* gxxx,gxxy,gxxz (这里是“升指标后的chi导数/chi”那类量，你沿用原变量名即可) */
        gxxx[i] = (gupxx[i] * chix[i] + gupxy[i] * chiy[i] + gupxz[i] * chiz[i]) / chin1[i];
        gxxy[i] = (gupxy[i] * chix[i] + gupyy[i] * chiy[i] + gupyz[i] * chiz[i]) / chin1[i];
        gxxz[i] = (gupxz[i] * chix[i] + gupyz[i] * chiy[i] + gupzz[i] * chiz[i]) / chin1[i];

        /* Christoffel 修正项 */
        Gamxxx[i] = Gamxxx[i] - ( ((chix[i] + chix[i]) / chin1[i]) - gxx[i] * gxxx[i] ) * HALF;
        Gamyxx[i] = Gamyxx[i] - ( (0.0 / chin1[i])           - gxx[i] * gxxy[i] ) * HALF;  /* 原式只有 -gxx*gxxy */
        Gamzxx[i] = Gamzxx[i] - ( (0.0 / chin1[i])           - gxx[i] * gxxz[i] ) * HALF;

        Gamxyy[i] = Gamxyy[i] - ( (0.0 / chin1[i])           - gyy[i] * gxxx[i] ) * HALF;
        Gamyyy[i] = Gamyyy[i] - ( ((chiy[i] + chiy[i]) / chin1[i]) - gyy[i] * gxxy[i] ) * HALF;
        Gamzyy[i] = Gamzyy[i] - ( (0.0 / chin1[i])           - gyy[i] * gxxz[i] ) * HALF;

        Gamxzz[i] = Gamxzz[i] - ( (0.0 / chin1[i])           - gzz[i] * gxxx[i] ) * HALF;
        Gamyzz[i] = Gamyzz[i] - ( (0.0 / chin1[i])           - gzz[i] * gxxy[i] ) * HALF;
        Gamzzz[i] = Gamzzz[i] - ( ((chiz[i] + chiz[i]) / chin1[i]) - gzz[i] * gxxz[i] ) * HALF;

        Gamxxy[i] = Gamxxy[i] - ( ( chiy[i] / chin1[i])      - gxy[i] * gxxx[i] ) * HALF;
        Gamyxy[i] = Gamyxy[i] - ( ( chix[i] / chin1[i])      - gxy[i] * gxxy[i] ) * HALF;
        Gamzxy[i] = Gamzxy[i] - ( (0.0 / chin1[i])           - gxy[i] * gxxz[i] ) * HALF;

        Gamxxz[i] = Gamxxz[i] - ( ( chiz[i] / chin1[i])      - gxz[i] * gxxx[i] ) * HALF;
        Gamyxz[i] = Gamyxz[i] - ( (0.0 / chin1[i])           - gxz[i] * gxxy[i] ) * HALF;
        Gamzxz[i] = Gamzxz[i] - ( ( chix[i] / chin1[i])      - gxz[i] * gxxz[i] ) * HALF;

        Gamxyz[i] = Gamxyz[i] - ( (0.0 / chin1[i])           - gyz[i] * gxxx[i] ) * HALF;
        Gamyyz[i] = Gamyyz[i] - ( ( chiz[i] / chin1[i])      - gyz[i] * gxxy[i] ) * HALF;
        Gamzyz[i] = Gamzyz[i] - ( ( chiy[i] / chin1[i])      - gyz[i] * gxxz[i] ) * HALF;

        /* fxx..fyz 修正：减去 Γ * ∂Lap */
        fxx[i] = fxx[i] - Gamxxx[i] * Lapx[i] - Gamyxx[i] * Lapy[i] - Gamzxx[i] * Lapz[i];
        fyy[i] = fyy[i] - Gamxyy[i] * Lapx[i] - Gamyyy[i] * Lapy[i] - Gamzyy[i] * Lapz[i];
        fzz[i] = fzz[i] - Gamxzz[i] * Lapx[i] - Gamyzz[i] * Lapy[i] - Gamzzz[i] * Lapz[i];

        fxy[i] = fxy[i] - Gamxxy[i] * Lapx[i] - Gamyxy[i] * Lapy[i] - Gamzxy[i] * Lapz[i];
        fxz[i] = fxz[i] - Gamxxz[i] * Lapx[i] - Gamyxz[i] * Lapy[i] - Gamzxz[i] * Lapz[i];
        fyz[i] = fyz[i] - Gamxyz[i] * Lapx[i] - Gamyyz[i] * Lapy[i] - Gamzyz[i] * Lapz[i];
    }

    for (int i = 0; i < all; i += 1) {
        trK_rhs[i] =  gupxx[i] * fxx[i] + gupyy[i] * fyy[i] + gupzz[i] * fzz[i]
                +  TWO * ( gupxy[i] * fxy[i] + gupxz[i] * fxz[i] + gupyz[i] * fyz[i] );
    }

    for (int i = 0; i < all; i++) {

        S[i] = chin1[i] * (
                gupxx[i] * Sxx[i] + gupyy[i] * Syy[i] + gupzz[i] * Szz[i]
            + TWO * ( gupxy[i] * Sxy[i] + gupxz[i] * Sxz[i] + gupyz[i] * Syz[i] )
        );

        f[i] = F2o3 * trK[i] * trK[i]
            - (
                gupxx[i] * (
                    gupxx[i] * Axx[i] * Axx[i] + gupyy[i] * Axy[i] * Axy[i] + gupzz[i] * Axz[i] * Axz[i]
                    + TWO * ( gupxy[i] * Axx[i] * Axy[i] + gupxz[i] * Axx[i] * Axz[i] + gupyz[i] * Axy[i] * Axz[i] )
                )
                + gupyy[i] * (
                    gupxx[i] * Axy[i] * Axy[i] + gupyy[i] * Ayy[i] * Ayy[i] + gupzz[i] * Ayz[i] * Ayz[i]
                    + TWO * ( gupxy[i] * Axy[i] * Ayy[i] + gupxz[i] * Axy[i] * Ayz[i] + gupyz[i] * Ayy[i] * Ayz[i] )
                )
                + gupzz[i] * (
                    gupxx[i] * Axz[i] * Axz[i] + gupyy[i] * Ayz[i] * Ayz[i] + gupzz[i] * Azz[i] * Azz[i]
                    + TWO * ( gupxy[i] * Axz[i] * Ayz[i] + gupxz[i] * Axz[i] * Azz[i] + gupyz[i] * Ayz[i] * Azz[i] )
                )
                + TWO * (
                    gupxy[i] * (
                        gupxx[i] * Axx[i] * Axy[i] + gupyy[i] * Axy[i] * Ayy[i] + gupzz[i] * Axz[i] * Ayz[i]
                        + gupxy[i] * ( Axx[i] * Ayy[i] + Axy[i] * Axy[i] )
                        + gupxz[i] * ( Axx[i] * Ayz[i] + Axz[i] * Axy[i] )
                        + gupyz[i] * ( Axy[i] * Ayz[i] + Axz[i] * Ayy[i] )
                    )
                    + gupxz[i] * (
                        gupxx[i] * Axx[i] * Axz[i] + gupyy[i] * Axy[i] * Ayz[i] + gupzz[i] * Axz[i] * Azz[i]
                        + gupxy[i] * ( Axx[i] * Ayz[i] + Axy[i] * Axz[i] )
                        + gupxz[i] * ( Axx[i] * Azz[i] + Axz[i] * Axz[i] )
                        + gupyz[i] * ( Axy[i] * Azz[i] + Axz[i] * Ayz[i] )
                    )
                    + gupyz[i] * (
                        gupxx[i] * Axy[i] * Axz[i] + gupyy[i] * Ayy[i] * Ayz[i] + gupzz[i] * Ayz[i] * Azz[i]
                        + gupxy[i] * ( Axy[i] * Ayz[i] + Ayy[i] * Axz[i] )
                        + gupxz[i] * ( Axy[i] * Azz[i] + Ayz[i] * Axz[i] )
                        + gupyz[i] * ( Ayy[i] * Azz[i] + Ayz[i] * Ayz[i] )
                    )
                )
            )
            - F16 * PI * rho[i]
            + EIGHT * PI * S[i];

        f[i] = -F1o3 * (
                gupxx[i] * fxx[i] + gupyy[i] * fyy[i] + gupzz[i] * fzz[i]
            + TWO * ( gupxy[i] * fxy[i] + gupxz[i] * fxz[i] + gupyz[i] * fyz[i] )
            + (alpn1[i] / chin1[i]) * f[i]
        );

        fxx[i] = alpn1[i] * (Rxx[i] - EIGHT * PI * Sxx[i]) - fxx[i];
        fxy[i] = alpn1[i] * (Rxy[i] - EIGHT * PI * Sxy[i]) - fxy[i];
        fxz[i] = alpn1[i] * (Rxz[i] - EIGHT * PI * Sxz[i]) - fxz[i];
        fyy[i] = alpn1[i] * (Ryy[i] - EIGHT * PI * Syy[i]) - fyy[i];
        fyz[i] = alpn1[i] * (Ryz[i] - EIGHT * PI * Syz[i]) - fyz[i];
        fzz[i] = alpn1[i] * (Rzz[i] - EIGHT * PI * Szz[i]) - fzz[i];
    }

    for (int i = 0; i < all; i += 1) {

        /* Aij_rhs = fij - gij * f */
        Axx_rhs[i] = fxx[i] - gxx[i] * f[i];
        Ayy_rhs[i] = fyy[i] - gyy[i] * f[i];
        Azz_rhs[i] = fzz[i] - gzz[i] * f[i];
        Axy_rhs[i] = fxy[i] - gxy[i] * f[i];
        Axz_rhs[i] = fxz[i] - gxz[i] * f[i];
        Ayz_rhs[i] = fyz[i] - gyz[i] * f[i];

        /* Now: store A_il A^l_j into fij: */
        fxx[i] =
            gupxx[i] * Axx[i] * Axx[i]
            + gupyy[i] * Axy[i] * Axy[i]
            + gupzz[i] * Axz[i] * Axz[i]
            + TWO * ( gupxy[i] * Axx[i] * Axy[i]
                    + gupxz[i] * Axx[i] * Axz[i]
                    + gupyz[i] * Axy[i] * Axz[i] );

        fyy[i] =
            gupxx[i] * Axy[i] * Axy[i]
            + gupyy[i] * Ayy[i] * Ayy[i]
            + gupzz[i] * Ayz[i] * Ayz[i]
            + TWO * ( gupxy[i] * Axy[i] * Ayy[i]
                    + gupxz[i] * Axy[i] * Ayz[i]
                    + gupyz[i] * Ayy[i] * Ayz[i] );

        fzz[i] =
            gupxx[i] * Axz[i] * Axz[i]
            + gupyy[i] * Ayz[i] * Ayz[i]
            + gupzz[i] * Azz[i] * Azz[i]
            + TWO * ( gupxy[i] * Axz[i] * Ayz[i]
                    + gupxz[i] * Axz[i] * Azz[i]
                    + gupyz[i] * Ayz[i] * Azz[i] );

        fxy[i] =
            gupxx[i] * Axx[i] * Axy[i]
            + gupyy[i] * Axy[i] * Ayy[i]
            + gupzz[i] * Axz[i] * Ayz[i]
            + gupxy[i] * (Axx[i] * Ayy[i] + Axy[i] * Axy[i])
            + gupxz[i] * (Axx[i] * Ayz[i] + Axz[i] * Axy[i])
            + gupyz[i] * (Axy[i] * Ayz[i] + Axz[i] * Ayy[i]);

        fxz[i] =
            gupxx[i] * Axx[i] * Axz[i]
            + gupyy[i] * Axy[i] * Ayz[i]
            + gupzz[i] * Axz[i] * Azz[i]
            + gupxy[i] * (Axx[i] * Ayz[i] + Axy[i] * Axz[i])
            + gupxz[i] * (Axx[i] * Azz[i] + Axz[i] * Axz[i])
            + gupyz[i] * (Axy[i] * Azz[i] + Axz[i] * Ayz[i]);

        fyz[i] =
            gupxx[i] * Axy[i] * Axz[i]
            + gupyy[i] * Ayy[i] * Ayz[i]
            + gupzz[i] * Ayz[i] * Azz[i]
            + gupxy[i] * (Axy[i] * Ayz[i] + Ayy[i] * Axz[i])
            + gupxz[i] * (Axy[i] * Azz[i] + Ayz[i] * Axz[i])
            + gupyz[i] * (Ayy[i] * Azz[i] + Ayz[i] * Ayz[i]);

        /* f = chin1 */
        f[i] = chin1[i];

        /* store D^i D_i Lap in trK_rhs */
        trK_rhs[i] = f[i] * trK_rhs[i];

        /* rhs for Aij */
        Axx_rhs[i] =
            f[i] * Axx_rhs[i]
            + alpn1[i] * ( trK[i] * Axx[i] - TWO * fxx[i] )
            + TWO * ( Axx[i] * betaxx[i] + Axy[i] * betayx[i] + Axz[i] * betazx[i] )
            - F2o3 * Axx[i] * div_beta[i];

        Ayy_rhs[i] =
            f[i] * Ayy_rhs[i]
            + alpn1[i] * ( trK[i] * Ayy[i] - TWO * fyy[i] )
            + TWO * ( Axy[i] * betaxy[i] + Ayy[i] * betayy[i] + Ayz[i] * betazy[i] )
            - F2o3 * Ayy[i] * div_beta[i];

        Azz_rhs[i] =
            f[i] * Azz_rhs[i]
            + alpn1[i] * ( trK[i] * Azz[i] - TWO * fzz[i] )
            + TWO * ( Axz[i] * betaxz[i] + Ayz[i] * betayz[i] + Azz[i] * betazz[i] )
            - F2o3 * Azz[i] * div_beta[i];

        Axy_rhs[i] =
            f[i] * Axy_rhs[i]
            + alpn1[i] * ( trK[i] * Axy[i] - TWO * fxy[i] )
            + Axx[i] * betaxy[i]
            + Axz[i] * betazy[i]
            + Ayy[i] * betayx[i]
            + Ayz[i] * betazx[i]
            + F1o3 * Axy[i] * div_beta[i]
            - Axy[i] * betazz[i];

        Ayz_rhs[i] =
            f[i] * Ayz_rhs[i]
            + alpn1[i] * ( trK[i] * Ayz[i] - TWO * fyz[i] )
            + Axy[i] * betaxz[i]
            + Ayy[i] * betayz[i]
            + Axz[i] * betaxy[i]
            + Azz[i] * betazy[i]
            + F1o3 * Ayz[i] * div_beta[i]
            - Ayz[i] * betaxx[i];

        Axz_rhs[i] =
            f[i] * Axz_rhs[i]
            + alpn1[i] * ( trK[i] * Axz[i] - TWO * fxz[i] )
            + Axx[i] * betaxz[i]
            + Axy[i] * betayz[i]
            + Ayz[i] * betayx[i]
            + Azz[i] * betazx[i]
            + F1o3 * Axz[i] * div_beta[i]
            - Axz[i] * betayy[i];

        /* Compute trace of S_ij */
        S[i] =
            f[i] * ( gupxx[i] * Sxx[i] + gupyy[i] * Syy[i] + gupzz[i] * Szz[i]
                    + TWO * ( gupxy[i] * Sxy[i] + gupxz[i] * Sxz[i] + gupyz[i] * Syz[i] ) );

        /* rhs for trK */
        trK_rhs[i] =
            -trK_rhs[i]
            + alpn1[i] * ( F1o3 * trK[i] * trK[i]
                        + gupxx[i] * fxx[i] + gupyy[i] * fyy[i] + gupzz[i] * fzz[i]
                        + TWO * ( gupxy[i] * fxy[i] + gupxz[i] * fxz[i] + gupyz[i] * fyz[i] )
                        + FOUR * PI * ( rho[i] + S[i] ) );

        /* gauge variable part */
        Lap_rhs[i] = -TWO * alpn1[i] * trK[i];
    }
    for (int i = 0; i < all; i += 1) {

        betax_rhs[i] = FF * dtSfx[i];
        betay_rhs[i] = FF * dtSfy[i];
        betaz_rhs[i] = FF * dtSfz[i];
    }
    fderivs(ex,chi,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev);
    for(int i=0;i<all;i+=1){
        reta[i] =
            gupxx[i] * dtSfx_rhs[i] * dtSfx_rhs[i]
            + gupyy[i] * dtSfy_rhs[i] * dtSfy_rhs[i]
            + gupzz[i] * dtSfz_rhs[i] * dtSfz_rhs[i]
            + TWO * ( gupxy[i] * dtSfx_rhs[i] * dtSfy_rhs[i]
                    + gupxz[i] * dtSfx_rhs[i] * dtSfz_rhs[i]
                    + gupyz[i] * dtSfy_rhs[i] * dtSfz_rhs[i] );

        reta[i] = 1.31 / 2.0 * sqrt( reta[i] / chin1[i] ) / pow( (1.0 - sqrt(chin1[i])), 2.0 );

        dtSfx_rhs[i] = Gamx_rhs[i] - reta[i] * dtSfx[i];
        dtSfy_rhs[i] = Gamy_rhs[i] - reta[i] * dtSfy[i];
        dtSfz_rhs[i] = Gamz_rhs[i] - reta[i] * dtSfz[i];
    }

    SSS[0]=SYM;
    SSS[1]=SYM;
    SSS[2]=SYM;

    AAS[0] = ANTI;
    AAS[1] = ANTI;
    AAS[2] = SYM;

    ASA[0] = ANTI;
    ASA[1] = SYM;
    ASA[2] = ANTI;

    SAA[0] = SYM;
    SAA[1] = ANTI;
    SAA[2] = ANTI;

    ASS[0] = ANTI;
    ASS[1] = SYM;
    ASS[2] = SYM;

    SAS[0] = SYM;
    SAS[1] = ANTI;
    SAS[2] = SYM;

    SSA[0] = SYM;
    SSA[1] = SYM;
    SSA[2] = ANTI;

    lopsided(ex,X,Y,Z,gxx,gxx_rhs,betax,betay,betaz,Symmetry,SSS);
    lopsided(ex,X,Y,Z,gxy,gxy_rhs,betax,betay,betaz,Symmetry,AAS);
    lopsided(ex,X,Y,Z,gxz,gxz_rhs,betax,betay,betaz,Symmetry,ASA);
    lopsided(ex,X,Y,Z,gyy,gyy_rhs,betax,betay,betaz,Symmetry,SSS);
    lopsided(ex,X,Y,Z,gyz,gyz_rhs,betax,betay,betaz,Symmetry,SAA);
    lopsided(ex,X,Y,Z,gzz,gzz_rhs,betax,betay,betaz,Symmetry,SSS);

    lopsided(ex,X,Y,Z,Axx,Axx_rhs,betax,betay,betaz,Symmetry,SSS);
    lopsided(ex,X,Y,Z,Axy,Axy_rhs,betax,betay,betaz,Symmetry,AAS);
    lopsided(ex,X,Y,Z,Axz,Axz_rhs,betax,betay,betaz,Symmetry,ASA);
    lopsided(ex,X,Y,Z,Ayy,Ayy_rhs,betax,betay,betaz,Symmetry,SSS);
    lopsided(ex,X,Y,Z,Ayz,Ayz_rhs,betax,betay,betaz,Symmetry,SAA);
    lopsided(ex,X,Y,Z,Azz,Azz_rhs,betax,betay,betaz,Symmetry,SSS);

    lopsided(ex,X,Y,Z,chi,chi_rhs,betax,betay,betaz,Symmetry,SSS);
    lopsided(ex,X,Y,Z,trK,trK_rhs,betax,betay,betaz,Symmetry,SSS);

    lopsided(ex,X,Y,Z,Gamx,Gamx_rhs,betax,betay,betaz,Symmetry,ASS);
    lopsided(ex,X,Y,Z,Gamy,Gamy_rhs,betax,betay,betaz,Symmetry,SAS);
    lopsided(ex,X,Y,Z,Gamz,Gamz_rhs,betax,betay,betaz,Symmetry,SSA);
    lopsided(ex,X,Y,Z,Lap,Lap_rhs,betax,betay,betaz,Symmetry,SSS);
    lopsided(ex,X,Y,Z,betax,betax_rhs,betax,betay,betaz,Symmetry,ASS);
    lopsided(ex,X,Y,Z,betay,betay_rhs,betax,betay,betaz,Symmetry,SAS);
    lopsided(ex,X,Y,Z,betaz,betaz_rhs,betax,betay,betaz,Symmetry,SSA);
    lopsided(ex,X,Y,Z,dtSfx,dtSfx_rhs,betax,betay,betaz,Symmetry,ASS);
    lopsided(ex,X,Y,Z,dtSfy,dtSfy_rhs,betax,betay,betaz,Symmetry,SAS);
    lopsided(ex,X,Y,Z,dtSfz,dtSfz_rhs,betax,betay,betaz,Symmetry,SSA);

    if(eps>0){
        kodis(ex,X,Y,Z,chi,chi_rhs,SSS,Symmetry,eps);
        kodis(ex,X,Y,Z,trK,trK_rhs,SSS,Symmetry,eps);
        kodis(ex,X,Y,Z,dxx,gxx_rhs,SSS,Symmetry,eps);
        kodis(ex,X,Y,Z,gxy,gxy_rhs,AAS,Symmetry,eps);
        kodis(ex,X,Y,Z,gxz,gxz_rhs,ASA,Symmetry,eps);
        kodis(ex,X,Y,Z,dyy,gyy_rhs,SSS,Symmetry,eps);
        kodis(ex,X,Y,Z,gyz,gyz_rhs,SAA,Symmetry,eps);
        kodis(ex,X,Y,Z,dzz,gzz_rhs,SSS,Symmetry,eps);
        kodis(ex,X,Y,Z,Axx,Axx_rhs,SSS,Symmetry,eps);
        kodis(ex,X,Y,Z,Axy,Axy_rhs,AAS,Symmetry,eps);
        kodis(ex,X,Y,Z,Axz,Axz_rhs,ASA,Symmetry,eps);
        kodis(ex,X,Y,Z,Ayy,Ayy_rhs,SSS,Symmetry,eps);
        kodis(ex,X,Y,Z,Ayz,Ayz_rhs,SAA,Symmetry,eps);
        kodis(ex,X,Y,Z,Azz,Azz_rhs,SSS,Symmetry,eps);
        kodis(ex,X,Y,Z,Gamx,Gamx_rhs,ASS,Symmetry,eps);
        kodis(ex,X,Y,Z,Gamy,Gamy_rhs,SAS,Symmetry,eps);
        kodis(ex,X,Y,Z,Gamz,Gamz_rhs,SSA,Symmetry,eps);
        kodis(ex,X,Y,Z,Lap,Lap_rhs,SSS,Symmetry,eps);
        kodis(ex,X,Y,Z,betax,betax_rhs,ASS,Symmetry,eps);
        kodis(ex,X,Y,Z,betay,betay_rhs,SAS,Symmetry,eps);
        kodis(ex,X,Y,Z,betaz,betaz_rhs,SSA,Symmetry,eps);
        kodis(ex,X,Y,Z,dtSfx,dtSfx_rhs,ASS,Symmetry,eps);
        kodis(ex,X,Y,Z,dtSfy,dtSfy_rhs,SAS,Symmetry,eps);
        kodis(ex,X,Y,Z,dtSfz,dtSfz_rhs,SSA,Symmetry,eps);
    }
    if(co==0){
        for (int i = 0; i < all; i++) {
            ham_Res[i] =
                gupxx[i] * Rxx[i] + gupyy[i] * Ryy[i] + gupzz[i] * Rzz[i]
                + TWO * ( gupxy[i] * Rxy[i] + gupxz[i] * Rxz[i] + gupyz[i] * Ryz[i] );
            ham_Res[i] =
                chin1[i] * ham_Res[i]
                + F2o3 * trK[i] * trK[i]
                - (
                    gupxx[i] * (
                        gupxx[i] * Axx[i] * Axx[i] + gupyy[i] * Axy[i] * Axy[i] + gupzz[i] * Axz[i] * Axz[i]
                    + TWO * ( gupxy[i] * Axx[i] * Axy[i] + gupxz[i] * Axx[i] * Axz[i] + gupyz[i] * Axy[i] * Axz[i] )
                    )
                + gupyy[i] * (
                        gupxx[i] * Axy[i] * Axy[i] + gupyy[i] * Ayy[i] * Ayy[i] + gupzz[i] * Ayz[i] * Ayz[i]
                    + TWO * ( gupxy[i] * Axy[i] * Ayy[i] + gupxz[i] * Axy[i] * Ayz[i] + gupyz[i] * Ayy[i] * Ayz[i] )
                    )
                + gupzz[i] * (
                        gupxx[i] * Axz[i] * Axz[i] + gupyy[i] * Ayz[i] * Ayz[i] + gupzz[i] * Azz[i] * Azz[i]
                    + TWO * ( gupxy[i] * Axz[i] * Ayz[i] + gupxz[i] * Axz[i] * Azz[i] + gupyz[i] * Ayz[i] * Azz[i] )
                    )
                + TWO * (
                        gupxy[i] * (
                            gupxx[i] * Axx[i] * Axy[i] + gupyy[i] * Axy[i] * Ayy[i] + gupzz[i] * Axz[i] * Ayz[i]
                        + gupxy[i] * ( Axx[i] * Ayy[i] + Axy[i] * Axy[i] )
                        + gupxz[i] * ( Axx[i] * Ayz[i] + Axz[i] * Axy[i] )
                        + gupyz[i] * ( Axy[i] * Ayz[i] + Axz[i] * Ayy[i] )
                        )
                    + gupxz[i] * (
                            gupxx[i] * Axx[i] * Axz[i] + gupyy[i] * Axy[i] * Ayz[i] + gupzz[i] * Axz[i] * Azz[i]
                        + gupxy[i] * ( Axx[i] * Ayz[i] + Axy[i] * Axz[i] )
                        + gupxz[i] * ( Axx[i] * Azz[i] + Axz[i] * Axz[i] )
                        + gupyz[i] * ( Axy[i] * Azz[i] + Axz[i] * Ayz[i] )
                        )
                    + gupyz[i] * (
                            gupxx[i] * Axy[i] * Axz[i] + gupyy[i] * Ayy[i] * Ayz[i] + gupzz[i] * Ayz[i] * Azz[i]
                        + gupxy[i] * ( Axy[i] * Ayz[i] + Ayy[i] * Axz[i] )
                        + gupxz[i] * ( Axy[i] * Azz[i] + Ayz[i] * Axz[i] )
                        + gupyz[i] * ( Ayy[i] * Azz[i] + Ayz[i] * Ayz[i] )
                        )
                    )
                )
                - F16 * PI * rho[i];
        }
    }
    fderivs(ex,Axx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0);
    fderivs(ex,Axy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,0);
    fderivs(ex,Axz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,0);
    fderivs(ex,Ayy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0);
    fderivs(ex,Ayz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,0);
    fderivs(ex,Azz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0);

    return 0; // success
}
