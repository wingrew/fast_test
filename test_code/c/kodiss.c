#include "include/tool.h"

/*
 * C 版 kodis
 *
 * Fortran signature:
 * subroutine kodis(ex,X,Y,Z,f,f_rhs,SoA,Symmetry,eps)
 *
 * 约定：
 *   X: ex1, Y: ex2, Z: ex3
 *   f, f_rhs: ex1*ex2*ex3 按 idx_ex 布局
 *   SoA[3]
 *   eps: double
 */
void kodis(const int ex[3],
           const double *X, const double *Y, const double *Z,
           const double *f, double *f_rhs,
           const double SoA[3],
           int Symmetry, double eps)
{
    const double ONE = 1.0, SIX = 6.0, FIT = 15.0, TWT = 20.0;
    const double cof = 64.0;             // 2^6
    const int NO_SYMM = 0, OCTANT = 2;

    const int ex1 = ex[0], ex2 = ex[1], ex3 = ex[2];

    // Fortran: dX = X(2)-X(1) -> C: X[1]-X[0]
    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];
    (void)ONE; // ONE 在原 Fortran 里只是参数，这里不一定用得上

    // Fortran: imax=ex(1) 等是 1-based 上界
    const int imaxF = ex1;
    const int jmaxF = ex2;
    const int kmaxF = ex3;

    // Fortran: imin=jmin=kmin=1，某些对称情况变 -2
    int iminF = 1, jminF = 1, kminF = 1;

    if (Symmetry > NO_SYMM && fabs(Z[0]) < dZ) kminF = -2;
    if (Symmetry == OCTANT && fabs(X[0]) < dX) iminF = -2;
    if (Symmetry == OCTANT && fabs(Y[0]) < dY) jminF = -2;

    // 分配 fh：大小 (ex1+3)*(ex2+3)*(ex3+3)，对应 ord=3
    const size_t nx = (size_t)ex1 + 3;
    const size_t ny = (size_t)ex2 + 3;
    const size_t nz = (size_t)ex3 + 3;
    const size_t fh_size = nx * ny * nz;

    double *fh = (double*)malloc(fh_size * sizeof(double));
    if (!fh) return;

    // Fortran: call symmetry_bd(3,ex,f,fh,SoA)
    symmetry_bd(3, ex, f, fh, SoA);

    /*
     * Fortran loops:
     * do k=1,ex3
     * do j=1,ex2
     * do i=1,ex1
     *
     * C: k0=0..ex3-1, j0=0..ex2-1, i0=0..ex1-1
     * 并定义 Fortran index: iF=i0+1, ...
     */
    for (int k0 = 0; k0 < ex3; ++k0) {
        const int kF = k0 + 1;
        for (int j0 = 0; j0 < ex2; ++j0) {
            const int jF = j0 + 1;
            for (int i0 = 0; i0 < ex1; ++i0) {
                const int iF = i0 + 1;

                // Fortran if 条件：
                // i-3 >= imin .and. i+3 <= imax  等（都是 Fortran 索引）
                if ((iF - 3) >= iminF && (iF + 3) <= imaxF &&
                    (jF - 3) >= jminF && (jF + 3) <= jmaxF &&
                    (kF - 3) >= kminF && (kF + 3) <= kmaxF)
                {
                    const size_t p = idx_ex(i0, j0, k0, ex);

                    // 三个方向各一份同型的 7 点组合（实际上是对称的 6th-order dissipation/filter 核）
                    const double Dx_term =
                        ( (fh[idx_fh_F(iF - 3, jF, kF, ex)] + fh[idx_fh_F(iF + 3, jF, kF, ex)]) -
                          SIX * (fh[idx_fh_F(iF - 2, jF, kF, ex)] + fh[idx_fh_F(iF + 2, jF, kF, ex)]) +
                          FIT * (fh[idx_fh_F(iF - 1, jF, kF, ex)] + fh[idx_fh_F(iF + 1, jF, kF, ex)]) -
                          TWT *  fh[idx_fh_F(iF    , jF, kF, ex)] ) / dX;

                    const double Dy_term =
                        ( (fh[idx_fh_F(iF, jF - 3, kF, ex)] + fh[idx_fh_F(iF, jF + 3, kF, ex)]) -
                          SIX * (fh[idx_fh_F(iF, jF - 2, kF, ex)] + fh[idx_fh_F(iF, jF + 2, kF, ex)]) +
                          FIT * (fh[idx_fh_F(iF, jF - 1, kF, ex)] + fh[idx_fh_F(iF, jF + 1, kF, ex)]) -
                          TWT *  fh[idx_fh_F(iF, jF    , kF, ex)] ) / dY;

                    const double Dz_term =
                        ( (fh[idx_fh_F(iF, jF, kF - 3, ex)] + fh[idx_fh_F(iF, jF, kF + 3, ex)]) -
                          SIX * (fh[idx_fh_F(iF, jF, kF - 2, ex)] + fh[idx_fh_F(iF, jF, kF + 2, ex)]) +
                          FIT * (fh[idx_fh_F(iF, jF, kF - 1, ex)] + fh[idx_fh_F(iF, jF, kF + 1, ex)]) -
                          TWT *  fh[idx_fh_F(iF, jF, kF    , ex)] ) / dZ;

                    // Fortran:
                    // f_rhs(i,j,k) = f_rhs(i,j,k) + eps/cof*(Dx_term + Dy_term + Dz_term)
                    f_rhs[p] += (eps / cof) * (Dx_term + Dy_term + Dz_term);
                }
            }
        }
    }

    free(fh);
}