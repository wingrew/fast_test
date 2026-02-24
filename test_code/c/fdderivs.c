#include "../include/tool.h"
void fdderivs(const int ex[3],
              const double *f,
              double *fxx, double *fxy, double *fxz,
              double *fyy, double *fyz, double *fzz,
              const double *X, const double *Y, const double *Z,
              double SYM1, double SYM2, double SYM3,
              int Symmetry, int onoff)
{
    (void)onoff;

    const int NO_SYMM = 0, EQ_SYMM = 1;
    const double ZEO = 0.0, ONE = 1.0, TWO = 2.0;
    const double F1o4   = 2.5e-1;          // 1/4
    const double F8     = 8.0;
    const double F16    = 16.0;
    const double F30    = 30.0;
    const double F1o12  = ONE / 12.0;
    const double F1o144 = ONE / 144.0;

    const int ex1 = ex[0], ex2 = ex[1], ex3 = ex[2];

    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

    const int imaxF = ex1;
    const int jmaxF = ex2;
    const int kmaxF = ex3;

    int iminF = 1, jminF = 1, kminF = 1;
    if (Symmetry > NO_SYMM && fabs(Z[0]) < dZ) kminF = -1;
    if (Symmetry > EQ_SYMM && fabs(X[0]) < dX) iminF = -1;
    if (Symmetry > EQ_SYMM && fabs(Y[0]) < dY) jminF = -1;

    const double SoA[3] = { SYM1, SYM2, SYM3 };

    /* fh: (ex1+2)*(ex2+2)*(ex3+2) because ord=2 */
    const size_t nx = (size_t)ex1 + 2;
    const size_t ny = (size_t)ex2 + 2;
    const size_t nz = (size_t)ex3 + 2;
    const size_t fh_size = nx * ny * nz;

    double *fh = (double*)malloc(fh_size * sizeof(double));
    if (!fh) return;

    symmetry_bd(2, ex, f, fh, SoA);

    /* 系数：按 Fortran 原式 */
    const double Sdxdx = ONE / (dX * dX);
    const double Sdydy = ONE / (dY * dY);
    const double Sdzdz = ONE / (dZ * dZ);

    const double Fdxdx = F1o12 / (dX * dX);
    const double Fdydy = F1o12 / (dY * dY);
    const double Fdzdz = F1o12 / (dZ * dZ);

    const double Sdxdy = F1o4 / (dX * dY);
    const double Sdxdz = F1o4 / (dX * dZ);
    const double Sdydz = F1o4 / (dY * dZ);

    const double Fdxdy = F1o144 / (dX * dY);
    const double Fdxdz = F1o144 / (dX * dZ);
    const double Fdydz = F1o144 / (dY * dZ);

    /* 输出清零：fxx,fyy,fzz,fxy,fxz,fyz = 0 */
    const size_t all = (size_t)ex1 * (size_t)ex2 * (size_t)ex3;
    for (size_t p = 0; p < all; ++p) {
        fxx[p] = ZEO; fyy[p] = ZEO; fzz[p] = ZEO;
        fxy[p] = ZEO; fxz[p] = ZEO; fyz[p] = ZEO;
    }

    /*
     * Fortran:
     * do k=1,ex3-1
     * do j=1,ex2-1
     * do i=1,ex1-1
     */
    for (int k0 = 0; k0 <= ex3 - 2; ++k0) {
        const int kF = k0 + 1;
        for (int j0 = 0; j0 <= ex2 - 2; ++j0) {
            const int jF = j0 + 1;
            for (int i0 = 0; i0 <= ex1 - 2; ++i0) {
                const int iF = i0 + 1;
                const size_t p = idx_ex(i0, j0, k0, ex);

                /* 高阶分支：i±2,j±2,k±2 都在范围内 */
                if ((iF + 2) <= imaxF && (iF - 2) >= iminF &&
                    (jF + 2) <= jmaxF && (jF - 2) >= jminF &&
                    (kF + 2) <= kmaxF && (kF - 2) >= kminF)
                {
                    fxx[p] = Fdxdx * (
                        -fh[idx_fh_F_ord2(iF - 2, jF,     kF,     ex)] +
                        F16 * fh[idx_fh_F_ord2(iF - 1, jF,     kF,     ex)] -
                        F30 * fh[idx_fh_F_ord2(iF,     jF,     kF,     ex)] -
                        fh[idx_fh_F_ord2(iF + 2, jF,     kF,     ex)] +
                        F16 * fh[idx_fh_F_ord2(iF + 1, jF,     kF,     ex)]
                    );

                    fyy[p] = Fdydy * (
                        -fh[idx_fh_F_ord2(iF,     jF - 2, kF,     ex)] +
                        F16 * fh[idx_fh_F_ord2(iF,     jF - 1, kF,     ex)] -
                        F30 * fh[idx_fh_F_ord2(iF,     jF,     kF,     ex)] -
                        fh[idx_fh_F_ord2(iF,     jF + 2, kF,     ex)] +
                        F16 * fh[idx_fh_F_ord2(iF,     jF + 1, kF,     ex)]
                    );

                    fzz[p] = Fdzdz * (
                        -fh[idx_fh_F_ord2(iF,     jF,     kF - 2, ex)] +
                        F16 * fh[idx_fh_F_ord2(iF,     jF,     kF - 1, ex)] -
                        F30 * fh[idx_fh_F_ord2(iF,     jF,     kF,     ex)] -
                        fh[idx_fh_F_ord2(iF,     jF,     kF + 2, ex)] +
                        F16 * fh[idx_fh_F_ord2(iF,     jF,     kF + 1, ex)]
                    );

                    /* fxy 高阶：完全照搬 Fortran 的括号结构 */
                    {
                        const double t_jm2 =
                            ( fh[idx_fh_F_ord2(iF - 2, jF - 2, kF, ex)]
                             -F8*fh[idx_fh_F_ord2(iF - 1, jF - 2, kF, ex)]
                             +F8*fh[idx_fh_F_ord2(iF + 1, jF - 2, kF, ex)]
                             -    fh[idx_fh_F_ord2(iF + 2, jF - 2, kF, ex)] );

                        const double t_jm1 =
                            ( fh[idx_fh_F_ord2(iF - 2, jF - 1, kF, ex)]
                             -F8*fh[idx_fh_F_ord2(iF - 1, jF - 1, kF, ex)]
                             +F8*fh[idx_fh_F_ord2(iF + 1, jF - 1, kF, ex)]
                             -    fh[idx_fh_F_ord2(iF + 2, jF - 1, kF, ex)] );

                        const double t_jp1 =
                            ( fh[idx_fh_F_ord2(iF - 2, jF + 1, kF, ex)]
                             -F8*fh[idx_fh_F_ord2(iF - 1, jF + 1, kF, ex)]
                             +F8*fh[idx_fh_F_ord2(iF + 1, jF + 1, kF, ex)]
                             -    fh[idx_fh_F_ord2(iF + 2, jF + 1, kF, ex)] );

                        const double t_jp2 =
                            ( fh[idx_fh_F_ord2(iF - 2, jF + 2, kF, ex)]
                             -F8*fh[idx_fh_F_ord2(iF - 1, jF + 2, kF, ex)]
                             +F8*fh[idx_fh_F_ord2(iF + 1, jF + 2, kF, ex)]
                             -    fh[idx_fh_F_ord2(iF + 2, jF + 2, kF, ex)] );

                        fxy[p] = Fdxdy * ( t_jm2 - F8 * t_jm1 + F8 * t_jp1 - t_jp2 );
                    }

                    /* fxz 高阶 */
                    {
                        const double t_km2 =
                            ( fh[idx_fh_F_ord2(iF - 2, jF, kF - 2, ex)]
                             -F8*fh[idx_fh_F_ord2(iF - 1, jF, kF - 2, ex)]
                             +F8*fh[idx_fh_F_ord2(iF + 1, jF, kF - 2, ex)]
                             -    fh[idx_fh_F_ord2(iF + 2, jF, kF - 2, ex)] );

                        const double t_km1 =
                            ( fh[idx_fh_F_ord2(iF - 2, jF, kF - 1, ex)]
                             -F8*fh[idx_fh_F_ord2(iF - 1, jF, kF - 1, ex)]
                             +F8*fh[idx_fh_F_ord2(iF + 1, jF, kF - 1, ex)]
                             -    fh[idx_fh_F_ord2(iF + 2, jF, kF - 1, ex)] );

                        const double t_kp1 =
                            ( fh[idx_fh_F_ord2(iF - 2, jF, kF + 1, ex)]
                             -F8*fh[idx_fh_F_ord2(iF - 1, jF, kF + 1, ex)]
                             +F8*fh[idx_fh_F_ord2(iF + 1, jF, kF + 1, ex)]
                             -    fh[idx_fh_F_ord2(iF + 2, jF, kF + 1, ex)] );

                        const double t_kp2 =
                            ( fh[idx_fh_F_ord2(iF - 2, jF, kF + 2, ex)]
                             -F8*fh[idx_fh_F_ord2(iF - 1, jF, kF + 2, ex)]
                             +F8*fh[idx_fh_F_ord2(iF + 1, jF, kF + 2, ex)]
                             -    fh[idx_fh_F_ord2(iF + 2, jF, kF + 2, ex)] );

                        fxz[p] = Fdxdz * ( t_km2 - F8 * t_km1 + F8 * t_kp1 - t_kp2 );
                    }

                    /* fyz 高阶 */
                    {
                        const double t_km2 =
                            ( fh[idx_fh_F_ord2(iF, jF - 2, kF - 2, ex)]
                             -F8*fh[idx_fh_F_ord2(iF, jF - 1, kF - 2, ex)]
                             +F8*fh[idx_fh_F_ord2(iF, jF + 1, kF - 2, ex)]
                             -    fh[idx_fh_F_ord2(iF, jF + 2, kF - 2, ex)] );

                        const double t_km1 =
                            ( fh[idx_fh_F_ord2(iF, jF - 2, kF - 1, ex)]
                             -F8*fh[idx_fh_F_ord2(iF, jF - 1, kF - 1, ex)]
                             +F8*fh[idx_fh_F_ord2(iF, jF + 1, kF - 1, ex)]
                             -    fh[idx_fh_F_ord2(iF, jF + 2, kF - 1, ex)] );

                        const double t_kp1 =
                            ( fh[idx_fh_F_ord2(iF, jF - 2, kF + 1, ex)]
                             -F8*fh[idx_fh_F_ord2(iF, jF - 1, kF + 1, ex)]
                             +F8*fh[idx_fh_F_ord2(iF, jF + 1, kF + 1, ex)]
                             -    fh[idx_fh_F_ord2(iF, jF + 2, kF + 1, ex)] );

                        const double t_kp2 =
                            ( fh[idx_fh_F_ord2(iF, jF - 2, kF + 2, ex)]
                             -F8*fh[idx_fh_F_ord2(iF, jF - 1, kF + 2, ex)]
                             +F8*fh[idx_fh_F_ord2(iF, jF + 1, kF + 2, ex)]
                             -    fh[idx_fh_F_ord2(iF, jF + 2, kF + 2, ex)] );

                        fyz[p] = Fdydz * ( t_km2 - F8 * t_km1 + F8 * t_kp1 - t_kp2 );
                    }
                }
                /* 二阶分支：i±1,j±1,k±1 在范围内 */
                else if ((iF + 1) <= imaxF && (iF - 1) >= iminF &&
                         (jF + 1) <= jmaxF && (jF - 1) >= jminF &&
                         (kF + 1) <= kmaxF && (kF - 1) >= kminF)
                {
                    fxx[p] = Sdxdx * (
                        fh[idx_fh_F_ord2(iF - 1, jF,     kF,     ex)] -
                        TWO * fh[idx_fh_F_ord2(iF,     jF,     kF,     ex)] +
                        fh[idx_fh_F_ord2(iF + 1, jF,     kF,     ex)]
                    );

                    fyy[p] = Sdydy * (
                        fh[idx_fh_F_ord2(iF,     jF - 1, kF,     ex)] -
                        TWO * fh[idx_fh_F_ord2(iF,     jF,     kF,     ex)] +
                        fh[idx_fh_F_ord2(iF,     jF + 1, kF,     ex)]
                    );

                    fzz[p] = Sdzdz * (
                        fh[idx_fh_F_ord2(iF,     jF,     kF - 1, ex)] -
                        TWO * fh[idx_fh_F_ord2(iF,     jF,     kF,     ex)] +
                        fh[idx_fh_F_ord2(iF,     jF,     kF + 1, ex)]
                    );

                    fxy[p] = Sdxdy * (
                        fh[idx_fh_F_ord2(iF - 1, jF - 1, kF, ex)] -
                        fh[idx_fh_F_ord2(iF + 1, jF - 1, kF, ex)] -
                        fh[idx_fh_F_ord2(iF - 1, jF + 1, kF, ex)] +
                        fh[idx_fh_F_ord2(iF + 1, jF + 1, kF, ex)]
                    );

                    fxz[p] = Sdxdz * (
                        fh[idx_fh_F_ord2(iF - 1, jF, kF - 1, ex)] -
                        fh[idx_fh_F_ord2(iF + 1, jF, kF - 1, ex)] -
                        fh[idx_fh_F_ord2(iF - 1, jF, kF + 1, ex)] +
                        fh[idx_fh_F_ord2(iF + 1, jF, kF + 1, ex)]
                    );

                    fyz[p] = Sdydz * (
                        fh[idx_fh_F_ord2(iF, jF - 1, kF - 1, ex)] -
                        fh[idx_fh_F_ord2(iF, jF + 1, kF - 1, ex)] -
                        fh[idx_fh_F_ord2(iF, jF - 1, kF + 1, ex)] +
                        fh[idx_fh_F_ord2(iF, jF + 1, kF + 1, ex)]
                    );
                }
            }
        }
    }

    free(fh);
}