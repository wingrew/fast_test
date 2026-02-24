#include "share_func.h"
void fdderivs(const int ex[3],
              const double *f,
              double *fxx, double *fxy, double *fxz,
              double *fyy, double *fyz, double *fzz,
              const double *X, const double *Y, const double *Z,
              double SYM1, double SYM2, double SYM3,
              int Symmetry, int onoff);

void fderivs(const int ex[3],
             const double *f,
             double *fx, double *fy, double *fz,
             const double *X, const double *Y, const double *Z,
             double SYM1, double SYM2, double SYM3,
             int Symmetry, int onoff);

void kodis(const int ex[3],
           const double *X, const double *Y, const double *Z,
           const double *f, double *f_rhs,
           const double SoA[3],
           int Symmetry, double eps);

void lopsided(const int ex[3],
              const double *X, const double *Y, const double *Z,
              const double *f, double *f_rhs,
              const double *Sfx, const double *Sfy, const double *Sfz,
              int Symmetry, const double SoA[3]);