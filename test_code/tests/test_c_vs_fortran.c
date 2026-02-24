#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../c/include/tool.h"

static size_t idx3(int i, int j, int k, const int ex[3]) {
    return (size_t)i + (size_t)ex[0] * ((size_t)j + (size_t)ex[1] * (size_t)k);
}

static void fill_grid(double *X, int n, double h) {
    for (int i = 0; i < n; ++i) X[i] = i * h;
}

static void fill_field(double *f, const int ex[3]) {
    for (int k = 0; k < ex[2]; ++k)
        for (int j = 0; j < ex[1]; ++j)
            for (int i = 0; i < ex[0]; ++i) {
                size_t p = idx3(i, j, k, ex);
                f[p] = sin(0.31 * i) + cos(0.27 * j) + 0.13 * k + 0.01 * i * j;
            }
}

static void fill_shift(double *s, const int ex[3], double phase) {
    for (int k = 0; k < ex[2]; ++k)
        for (int j = 0; j < ex[1]; ++j)
            for (int i = 0; i < ex[0]; ++i) {
                size_t p = idx3(i, j, k, ex);
                s[p] = sin(phase + 0.2 * i - 0.1 * j + 0.05 * k);
            }
}

static double max_abs_diff(const double *a, const double *b, size_t n) {
    double m = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double d = fabs(a[i] - b[i]);
        if (d > m) m = d;
    }
    return m;
}

#ifdef TEST_FDERIVS
extern void fderivs_(const int *ex, const double *f, double *fx, double *fy, double *fz,
                     const double *X, const double *Y, const double *Z,
                     const double *SYM1, const double *SYM2, const double *SYM3,
                     const int *Symmetry, const int *onoff);

int main(void) {
    const int ex[3] = {10, 9, 8};
    const size_t n = (size_t)ex[0] * ex[1] * ex[2];
    double *X = calloc(ex[0], sizeof(double)), *Y = calloc(ex[1], sizeof(double)), *Z = calloc(ex[2], sizeof(double));
    double *f = calloc(n, sizeof(double));
    double *fx_c = calloc(n, sizeof(double)), *fy_c = calloc(n, sizeof(double)), *fz_c = calloc(n, sizeof(double));
    double *fx_f = calloc(n, sizeof(double)), *fy_f = calloc(n, sizeof(double)), *fz_f = calloc(n, sizeof(double));
    fill_grid(X, ex[0], 0.1); fill_grid(Y, ex[1], 0.2); fill_grid(Z, ex[2], 0.15); fill_field(f, ex);
    double s1 = 1.0, s2 = -1.0, s3 = 1.0; int sym = 2, onoff = 1;
    fderivs(ex, f, fx_c, fy_c, fz_c, X, Y, Z, s1, s2, s3, sym, onoff);
    fderivs_(ex, f, fx_f, fy_f, fz_f, X, Y, Z, &s1, &s2, &s3, &sym, &onoff);
    double d = fmax(max_abs_diff(fx_c, fx_f, n), fmax(max_abs_diff(fy_c, fy_f, n), max_abs_diff(fz_c, fz_f, n)));
    printf("fderivs max_abs_diff = %.3e\n", d);
    return d < 1e-12 ? 0 : 1;
}
#endif

#ifdef TEST_FDDERIVS
extern void fdderivs_(const int *ex, const double *f, double *fxx, double *fxy, double *fxz,
                      double *fyy, double *fyz, double *fzz,
                      const double *X, const double *Y, const double *Z,
                      const double *SYM1, const double *SYM2, const double *SYM3,
                      const int *Symmetry, const int *onoff);

int main(void) {
    const int ex[3] = {11, 10, 9};
    const size_t n = (size_t)ex[0] * ex[1] * ex[2];
    double *X = calloc(ex[0], sizeof(double)), *Y = calloc(ex[1], sizeof(double)), *Z = calloc(ex[2], sizeof(double));
    double *f = calloc(n, sizeof(double));
    double *c[6], *ff[6];
    for (int i = 0; i < 6; ++i) { c[i] = calloc(n, sizeof(double)); ff[i] = calloc(n, sizeof(double)); }
    fill_grid(X, ex[0], 0.1); fill_grid(Y, ex[1], 0.2); fill_grid(Z, ex[2], 0.12); fill_field(f, ex);
    double s1 = 1.0, s2 = 1.0, s3 = -1.0; int sym = 1, onoff = 1;
    fdderivs(ex, f, c[0], c[1], c[2], c[3], c[4], c[5], X, Y, Z, s1, s2, s3, sym, onoff);
    fdderivs_(ex, f, ff[0], ff[1], ff[2], ff[3], ff[4], ff[5], X, Y, Z, &s1, &s2, &s3, &sym, &onoff);
    double d = 0.0;
    for (int i = 0; i < 6; ++i) d = fmax(d, max_abs_diff(c[i], ff[i], n));
    printf("fdderivs max_abs_diff = %.3e\n", d);
    return d < 1e-12 ? 0 : 1;
}
#endif

#ifdef TEST_KODIS
extern void kodis_(const int *ex, const double *X, const double *Y, const double *Z,
                   const double *f, double *f_rhs, const double *SoA,
                   const int *Symmetry, const double *eps);

int main(void) {
    const int ex[3] = {12, 11, 10};
    const size_t n = (size_t)ex[0] * ex[1] * ex[2];
    double *X = calloc(ex[0], sizeof(double)), *Y = calloc(ex[1], sizeof(double)), *Z = calloc(ex[2], sizeof(double));
    double *f = calloc(n, sizeof(double)), *rhs_c = calloc(n, sizeof(double)), *rhs_f = calloc(n, sizeof(double));
    fill_grid(X, ex[0], 0.08); fill_grid(Y, ex[1], 0.09); fill_grid(Z, ex[2], 0.11); fill_field(f, ex);
    for (size_t i = 0; i < n; ++i) rhs_c[i] = rhs_f[i] = 0.5 * sin((double)i * 0.01);
    const double soa[3] = {1.0, -1.0, 1.0}; int sym = 2; double eps = 0.03;
    kodis(ex, X, Y, Z, f, rhs_c, soa, sym, eps);
    kodis_(ex, X, Y, Z, f, rhs_f, soa, &sym, &eps);
    double d = max_abs_diff(rhs_c, rhs_f, n);
    printf("kodis max_abs_diff = %.3e\n", d);
    return d < 1e-12 ? 0 : 1;
}
#endif

#ifdef TEST_LOPSIDED
extern void lopsided_(const int *ex, const double *X, const double *Y, const double *Z,
                      const double *f, double *f_rhs,
                      const double *Sfx, const double *Sfy, const double *Sfz,
                      const int *Symmetry, const double *SoA);

int main(void) {
    const int ex[3] = {12, 10, 9};
    const size_t n = (size_t)ex[0] * ex[1] * ex[2];
    double *X = calloc(ex[0], sizeof(double)), *Y = calloc(ex[1], sizeof(double)), *Z = calloc(ex[2], sizeof(double));
    double *f = calloc(n, sizeof(double)), *rhs_c = calloc(n, sizeof(double)), *rhs_f = calloc(n, sizeof(double));
    double *sfx = calloc(n, sizeof(double)), *sfy = calloc(n, sizeof(double)), *sfz = calloc(n, sizeof(double));
    fill_grid(X, ex[0], 0.1); fill_grid(Y, ex[1], 0.1); fill_grid(Z, ex[2], 0.1); fill_field(f, ex);
    fill_shift(sfx, ex, 0.1); fill_shift(sfy, ex, 0.7); fill_shift(sfz, ex, -0.5);
    for (size_t i = 0; i < n; ++i) rhs_c[i] = rhs_f[i] = 0.25 * cos((double)i * 0.02);
    const double soa[3] = {1.0, 1.0, -1.0}; int sym = 2;
    lopsided(ex, X, Y, Z, f, rhs_c, sfx, sfy, sfz, sym, soa);
    lopsided_(ex, X, Y, Z, f, rhs_f, sfx, sfy, sfz, &sym, soa);
    double d = max_abs_diff(rhs_c, rhs_f, n);
    printf("lopsided max_abs_diff = %.3e\n", d);
    return d < 1e-12 ? 0 : 1;
}
#endif
