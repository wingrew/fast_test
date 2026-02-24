#include "include/tool.h"
/*
 * 你需要提供 symmetry_bd 的 C 版本（或 Fortran 绑到 C 的接口）。
 * Fortran: call symmetry_bd(3,ex,f,fh,SoA)
 *
 * 约定：
 *   nghost = 3
 *   ex[3]  = {ex1,ex2,ex3}
 *   f      = 原始网格 (ex1*ex2*ex3)
 *   fh     = 扩展网格 ((ex1+3)*(ex2+3)*(ex3+3))，对应 Fortran 的 (-2:ex1, ...)
 *   SoA[3] = 输入参数
 */
void lopsided(const int ex[3],
              const double *X, const double *Y, const double *Z,
              const double *f, double *f_rhs,
              const double *Sfx, const double *Sfy, const double *Sfz,
              int Symmetry, const double SoA[3])
{
    const double ZEO = 0.0, ONE = 1.0, F3 = 3.0;
    const double TWO = 2.0, F6 = 6.0, F18 = 18.0;
    const double F12 = 12.0, F10 = 10.0, EIT = 8.0;

    const int NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2;
    (void)OCTANT; // 这里和 Fortran 一样只是定义了不用也没关系

    const int ex1 = ex[0], ex2 = ex[1], ex3 = ex[2];

    // 对应 Fortran: dX = X(2)-X(1)  （Fortran 1-based）
    // C: X[1]-X[0]
    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

    const double d12dx = ONE / F12 / dX;
    const double d12dy = ONE / F12 / dY;
    const double d12dz = ONE / F12 / dZ;

    // Fortran 里算了 d2dx/d2dy/d2dz 但本 subroutine 里没用到（保持一致也算出来）
    const double d2dx  = ONE / TWO / dX;
    const double d2dy  = ONE / TWO / dY;
    const double d2dz  = ONE / TWO / dZ;
    (void)d2dx; (void)d2dy; (void)d2dz;

    // Fortran:
    // imax = ex(1); jmax = ex(2); kmax = ex(3)
    const int imaxF = ex1;
    const int jmaxF = ex2;
    const int kmaxF = ex3;

    // Fortran:
    // imin=jmin=kmin=1; 若满足对称条件则设为 -2
    int iminF = 1, jminF = 1, kminF = 1;
    if (Symmetry > NO_SYMM && fabs(Z[0]) < dZ) kminF = -2;
    if (Symmetry > EQ_SYMM && fabs(X[0]) < dX) iminF = -2;
    if (Symmetry > EQ_SYMM && fabs(Y[0]) < dY) jminF = -2;

    // 分配 fh：大小 (ex1+3)*(ex2+3)*(ex3+3)
    const size_t nx = (size_t)ex1 + 3;
    const size_t ny = (size_t)ex2 + 3;
    const size_t nz = (size_t)ex3 + 3;
    const size_t fh_size = nx * ny * nz;

    double *fh = (double*)malloc(fh_size * sizeof(double));
    if (!fh) return; // 内存不足：直接返回（你也可以改成 abort/报错）

    // Fortran: call symmetry_bd(3,ex,f,fh,SoA)
    symmetry_bd(3, ex, f, fh, SoA);

    /*
     * Fortran 主循环：
     * do k=1,ex(3)-1
     * do j=1,ex(2)-1
     * do i=1,ex(1)-1
     *
     * 转成 C 0-based：
     * k0 = 0..ex3-2, j0 = 0..ex2-2, i0 = 0..ex1-2
     *
     * 并且 Fortran 里的 i/j/k 在 fh 访问时，仍然是 Fortran 索引值：
     * iF=i0+1, jF=j0+1, kF=k0+1
     */
    for (int k0 = 0; k0 <= ex3 - 2; ++k0) {
        const int kF = k0 + 1;
        for (int j0 = 0; j0 <= ex2 - 2; ++j0) {
            const int jF = j0 + 1;
            for (int i0 = 0; i0 <= ex1 - 2; ++i0) {
                const int iF = i0 + 1;

                const size_t p = idx_ex(i0, j0, k0, ex);

                // ---------------- x direction ----------------
                const double sfx = Sfx[p];
                if (sfx > ZEO) {
                    // Fortran: if(i+3 <= imax)
                    // iF+3 <= ex1  <=> i0+4 <= ex1 <=> i0 <= ex1-4
                    if (i0 <= ex1 - 4) {
                        f_rhs[p] += sfx * d12dx *
                            (-F3  * fh[idx_fh_F(iF - 1, jF, kF, ex)]
                             -F10 * fh[idx_fh_F(iF    , jF, kF, ex)]
                             +F18 * fh[idx_fh_F(iF + 1, jF, kF, ex)]
                             -F6  * fh[idx_fh_F(iF + 2, jF, kF, ex)]
                             +      fh[idx_fh_F(iF + 3, jF, kF, ex)]);
                    }
                    // elseif(i+2 <= imax)  <=> i0 <= ex1-3
                    else if (i0 <= ex1 - 3) {
                        f_rhs[p] += sfx * d12dx *
                            ( fh[idx_fh_F(iF - 2, jF, kF, ex)]
                             -EIT * fh[idx_fh_F(iF - 1, jF, kF, ex)]
                             +EIT * fh[idx_fh_F(iF + 1, jF, kF, ex)]
                             -      fh[idx_fh_F(iF + 2, jF, kF, ex)]);
                    }
                    // elseif(i+1 <= imax)  <=> i0 <= ex1-2（循环里总成立）
                    else if (i0 <= ex1 - 2) {
                        f_rhs[p] -= sfx * d12dx *
                            (-F3  * fh[idx_fh_F(iF + 1, jF, kF, ex)]
                             -F10 * fh[idx_fh_F(iF    , jF, kF, ex)]
                             +F18 * fh[idx_fh_F(iF - 1, jF, kF, ex)]
                             -F6  * fh[idx_fh_F(iF - 2, jF, kF, ex)]
                             +      fh[idx_fh_F(iF - 3, jF, kF, ex)]);
                    }
                } else if (sfx < ZEO) {
                    // Fortran: if(i-3 >= imin)
                    // (iF-3) >= iminF  <=> (i0-2) >= iminF
                    if ((i0 - 2) >= iminF) {
                        f_rhs[p] -= sfx * d12dx *
                            (-F3  * fh[idx_fh_F(iF + 1, jF, kF, ex)]
                             -F10 * fh[idx_fh_F(iF    , jF, kF, ex)]
                             +F18 * fh[idx_fh_F(iF - 1, jF, kF, ex)]
                             -F6  * fh[idx_fh_F(iF - 2, jF, kF, ex)]
                             +      fh[idx_fh_F(iF - 3, jF, kF, ex)]);
                    }
                    // elseif(i-2 >= imin) <=> (i0-1) >= iminF
                    else if ((i0 - 1) >= iminF) {
                        f_rhs[p] += sfx * d12dx *
                            ( fh[idx_fh_F(iF - 2, jF, kF, ex)]
                             -EIT * fh[idx_fh_F(iF - 1, jF, kF, ex)]
                             +EIT * fh[idx_fh_F(iF + 1, jF, kF, ex)]
                             -      fh[idx_fh_F(iF + 2, jF, kF, ex)]);
                    }
                    // elseif(i-1 >= imin) <=> i0 >= iminF
                    else if (i0 >= iminF) {
                        f_rhs[p] += sfx * d12dx *
                            (-F3  * fh[idx_fh_F(iF - 1, jF, kF, ex)]
                             -F10 * fh[idx_fh_F(iF    , jF, kF, ex)]
                             +F18 * fh[idx_fh_F(iF + 1, jF, kF, ex)]
                             -F6  * fh[idx_fh_F(iF + 2, jF, kF, ex)]
                             +      fh[idx_fh_F(iF + 3, jF, kF, ex)]);
                    }
                }

                // ---------------- y direction ----------------
                const double sfy = Sfy[p];
                if (sfy > ZEO) {
                    // jF+3 <= ex2 <=> j0+4 <= ex2 <=> j0 <= ex2-4
                    if (j0 <= ex2 - 4) {
                        f_rhs[p] += sfy * d12dy *
                            (-F3  * fh[idx_fh_F(iF, jF - 1, kF, ex)]
                             -F10 * fh[idx_fh_F(iF, jF    , kF, ex)]
                             +F18 * fh[idx_fh_F(iF, jF + 1, kF, ex)]
                             -F6  * fh[idx_fh_F(iF, jF + 2, kF, ex)]
                             +      fh[idx_fh_F(iF, jF + 3, kF, ex)]);
                    } else if (j0 <= ex2 - 3) {
                        f_rhs[p] += sfy * d12dy *
                            ( fh[idx_fh_F(iF, jF - 2, kF, ex)]
                             -EIT * fh[idx_fh_F(iF, jF - 1, kF, ex)]
                             +EIT * fh[idx_fh_F(iF, jF + 1, kF, ex)]
                             -      fh[idx_fh_F(iF, jF + 2, kF, ex)]);
                    } else if (j0 <= ex2 - 2) {
                        f_rhs[p] -= sfy * d12dy *
                            (-F3  * fh[idx_fh_F(iF, jF + 1, kF, ex)]
                             -F10 * fh[idx_fh_F(iF, jF    , kF, ex)]
                             +F18 * fh[idx_fh_F(iF, jF - 1, kF, ex)]
                             -F6  * fh[idx_fh_F(iF, jF - 2, kF, ex)]
                             +      fh[idx_fh_F(iF, jF - 3, kF, ex)]);
                    }
                } else if (sfy < ZEO) {
                    if ((j0 - 2) >= jminF) {
                        f_rhs[p] -= sfy * d12dy *
                            (-F3  * fh[idx_fh_F(iF, jF + 1, kF, ex)]
                             -F10 * fh[idx_fh_F(iF, jF    , kF, ex)]
                             +F18 * fh[idx_fh_F(iF, jF - 1, kF, ex)]
                             -F6  * fh[idx_fh_F(iF, jF - 2, kF, ex)]
                             +      fh[idx_fh_F(iF, jF - 3, kF, ex)]);
                    } else if ((j0 - 1) >= jminF) {
                        f_rhs[p] += sfy * d12dy *
                            ( fh[idx_fh_F(iF, jF - 2, kF, ex)]
                             -EIT * fh[idx_fh_F(iF, jF - 1, kF, ex)]
                             +EIT * fh[idx_fh_F(iF, jF + 1, kF, ex)]
                             -      fh[idx_fh_F(iF, jF + 2, kF, ex)]);
                    } else if (j0 >= jminF) {
                        f_rhs[p] += sfy * d12dy *
                            (-F3  * fh[idx_fh_F(iF, jF - 1, kF, ex)]
                             -F10 * fh[idx_fh_F(iF, jF    , kF, ex)]
                             +F18 * fh[idx_fh_F(iF, jF + 1, kF, ex)]
                             -F6  * fh[idx_fh_F(iF, jF + 2, kF, ex)]
                             +      fh[idx_fh_F(iF, jF + 3, kF, ex)]);
                    }
                }

                // ---------------- z direction ----------------
                const double sfz = Sfz[p];
                if (sfz > ZEO) {
                    if (k0 <= ex3 - 4) {
                        f_rhs[p] += sfz * d12dz *
                            (-F3  * fh[idx_fh_F(iF, jF, kF - 1, ex)]
                             -F10 * fh[idx_fh_F(iF, jF, kF    , ex)]
                             +F18 * fh[idx_fh_F(iF, jF, kF + 1, ex)]
                             -F6  * fh[idx_fh_F(iF, jF, kF + 2, ex)]
                             +      fh[idx_fh_F(iF, jF, kF + 3, ex)]);
                    } else if (k0 <= ex3 - 3) {
                        f_rhs[p] += sfz * d12dz *
                            ( fh[idx_fh_F(iF, jF, kF - 2, ex)]
                             -EIT * fh[idx_fh_F(iF, jF, kF - 1, ex)]
                             +EIT * fh[idx_fh_F(iF, jF, kF + 1, ex)]
                             -      fh[idx_fh_F(iF, jF, kF + 2, ex)]);
                    } else if (k0 <= ex3 - 2) {
                        f_rhs[p] -= sfz * d12dz *
                            (-F3  * fh[idx_fh_F(iF, jF, kF + 1, ex)]
                             -F10 * fh[idx_fh_F(iF, jF, kF    , ex)]
                             +F18 * fh[idx_fh_F(iF, jF, kF - 1, ex)]
                             -F6  * fh[idx_fh_F(iF, jF, kF - 2, ex)]
                             +      fh[idx_fh_F(iF, jF, kF - 3, ex)]);
                    }
                } else if (sfz < ZEO) {
                    if ((k0 - 2) >= kminF) {
                        f_rhs[p] -= sfz * d12dz *
                            (-F3  * fh[idx_fh_F(iF, jF, kF + 1, ex)]
                             -F10 * fh[idx_fh_F(iF, jF, kF    , ex)]
                             +F18 * fh[idx_fh_F(iF, jF, kF - 1, ex)]
                             -F6  * fh[idx_fh_F(iF, jF, kF - 2, ex)]
                             +      fh[idx_fh_F(iF, jF, kF - 3, ex)]);
                    } else if ((k0 - 1) >= kminF) {
                        f_rhs[p] += sfz * d12dz *
                            ( fh[idx_fh_F(iF, jF, kF - 2, ex)]
                             -EIT * fh[idx_fh_F(iF, jF, kF - 1, ex)]
                             +EIT * fh[idx_fh_F(iF, jF, kF + 1, ex)]
                             -      fh[idx_fh_F(iF, jF, kF + 2, ex)]);
                    } else if (k0 >= kminF) {
                        f_rhs[p] += sfz * d12dz *
                            (-F3  * fh[idx_fh_F(iF, jF, kF - 1, ex)]
                             -F10 * fh[idx_fh_F(iF, jF, kF    , ex)]
                             +F18 * fh[idx_fh_F(iF, jF, kF + 1, ex)]
                             -F6  * fh[idx_fh_F(iF, jF, kF + 2, ex)]
                             +      fh[idx_fh_F(iF, jF, kF + 3, ex)]);
                    }
                }
            }
        }
    }
    free(fh);
}





