#include "../include/tool.h"

/*
 * C 版 fderivs
 *
 * Fortran:
 * subroutine fderivs(ex,f,fx,fy,fz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff)
 *
 * 约定：
 *   f, fx, fy, fz: ex1*ex2*ex3，按 idx_ex 布局
 *   X: ex1, Y: ex2, Z: ex3
 */
void fderivs(const int ex[3],
             const double *f,
             double *fx, double *fy, double *fz,
             const double *X, const double *Y, const double *Z,
             double SYM1, double SYM2, double SYM3,
             int Symmetry, int onoff)
{
    (void)onoff; // Fortran 里没用到

    const double ZEO = 0.0, ONE = 1.0;
    const double TWO = 2.0, EIT = 8.0;
    const double F12 = 12.0;

    const int NO_SYMM = 0, EQ_SYMM = 1; // OCTANT=2 在本子程序里不直接用

    const int ex1 = ex[0], ex2 = ex[1], ex3 = ex[2];

    // dX = X(2)-X(1) -> C: X[1]-X[0]
    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

    // Fortran 1-based bounds
    const int imaxF = ex1;
    const int jmaxF = ex2;
    const int kmaxF = ex3;

    int iminF = 1, jminF = 1, kminF = 1;
    if (Symmetry > NO_SYMM && fabs(Z[0]) < dZ) kminF = -1;
    if (Symmetry > EQ_SYMM && fabs(X[0]) < dX) iminF = -1;
    if (Symmetry > EQ_SYMM && fabs(Y[0]) < dY) jminF = -1;

    // SoA(1:3) = SYM1,SYM2,SYM3
    const double SoA[3] = { SYM1, SYM2, SYM3 };

    // fh: (ex1+2)*(ex2+2)*(ex3+2) because ord=2
    const size_t nx = (size_t)ex1 + 2;
    const size_t ny = (size_t)ex2 + 2;
    const size_t nz = (size_t)ex3 + 2;
    const size_t fh_size = nx * ny * nz;

    double *fh = (double*)malloc(fh_size * sizeof(double));
    if (!fh) return;

    // call symmetry_bd(2,ex,f,fh,SoA)
    symmetry_bd(2, ex, f, fh, SoA);

    const double d12dx = ONE / F12 / dX;
    const double d12dy = ONE / F12 / dY;
    const double d12dz = ONE / F12 / dZ;

    const double d2dx  = ONE / TWO / dX;
    const double d2dy  = ONE / TWO / dY;
    const double d2dz  = ONE / TWO / dZ;

    // fx = fy = fz = 0
    const size_t all = (size_t)ex1 * (size_t)ex2 * (size_t)ex3;
    for (size_t p = 0; p < all; ++p) {
        fx[p] = ZEO;
        fy[p] = ZEO;
        fz[p] = ZEO;
    }

    /*
     * Fortran loops:
     * do k=1,ex3-1
     * do j=1,ex2-1
     * do i=1,ex1-1
     *
     * C: k0=0..ex3-2, j0=0..ex2-2, i0=0..ex1-2
     */
    for (int k0 = 0; k0 <= ex3 - 2; ++k0) {
        const int kF = k0 + 1;
        for (int j0 = 0; j0 <= ex2 - 2; ++j0) {
            const int jF = j0 + 1;
            for (int i0 = 0; i0 <= ex1 - 2; ++i0) {
                const int iF = i0 + 1;
                const size_t p = idx_ex(i0, j0, k0, ex);

                // if(i+2 <= imax .and. i-2 >= imin ... )  (全是 Fortran 索引)
                if ((iF + 2) <= imaxF && (iF - 2) >= iminF &&
                    (jF + 2) <= jmaxF && (jF - 2) >= jminF &&
                    (kF + 2) <= kmaxF && (kF - 2) >= kminF)
                {
                    fx[p] = d12dx * (
                        fh[idx_fh_F_ord2(iF - 2, jF,     kF,     ex)] -
                        EIT * fh[idx_fh_F_ord2(iF - 1, jF,     kF,     ex)] +
                        EIT * fh[idx_fh_F_ord2(iF + 1, jF,     kF,     ex)] -
                        fh[idx_fh_F_ord2(iF + 2, jF,     kF,     ex)]
                    );

                    fy[p] = d12dy * (
                        fh[idx_fh_F_ord2(iF,     jF - 2, kF,     ex)] -
                        EIT * fh[idx_fh_F_ord2(iF,     jF - 1, kF,     ex)] +
                        EIT * fh[idx_fh_F_ord2(iF,     jF + 1, kF,     ex)] -
                        fh[idx_fh_F_ord2(iF,     jF + 2, kF,     ex)]
                    );

                    fz[p] = d12dz * (
                        fh[idx_fh_F_ord2(iF,     jF,     kF - 2, ex)] -
                        EIT * fh[idx_fh_F_ord2(iF,     jF,     kF - 1, ex)] +
                        EIT * fh[idx_fh_F_ord2(iF,     jF,     kF + 1, ex)] -
                        fh[idx_fh_F_ord2(iF,     jF,     kF + 2, ex)]
                    );
                }
                // elseif(i+1 <= imax .and. i-1 >= imin ...)
                else if ((iF + 1) <= imaxF && (iF - 1) >= iminF &&
                         (jF + 1) <= jmaxF && (jF - 1) >= jminF &&
                         (kF + 1) <= kmaxF && (kF - 1) >= kminF)
                {
                    fx[p] = d2dx * (
                        -fh[idx_fh_F_ord2(iF - 1, jF,     kF,     ex)] +
                         fh[idx_fh_F_ord2(iF + 1, jF,     kF,     ex)]
                    );

                    fy[p] = d2dy * (
                        -fh[idx_fh_F_ord2(iF,     jF - 1, kF,     ex)] +
                         fh[idx_fh_F_ord2(iF,     jF + 1, kF,     ex)]
                    );

                    fz[p] = d2dz * (
                        -fh[idx_fh_F_ord2(iF,     jF,     kF - 1, ex)] +
                         fh[idx_fh_F_ord2(iF,     jF,     kF + 1, ex)]
                    );
                }
            }
        }
    }

    free(fh);
}