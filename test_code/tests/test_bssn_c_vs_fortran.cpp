#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <string>

#include "../c/include/bssn_rhs_compute.h"

extern "C" int compute_rhs_bssn_(int *ex, double *T,
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
                       int *Symmetry, int *Lev, double *eps, int *co);

#define ALL_FIELDS(X) \
X(chi) X(trK) X(dxx) X(gxy) X(gxz) X(dyy) X(gyz) X(dzz) \
X(Axx) X(Axy) X(Axz) X(Ayy) X(Ayz) X(Azz) X(Gamx) X(Gamy) X(Gamz) \
X(Lap) X(betax) X(betay) X(betaz) X(dtSfx) X(dtSfy) X(dtSfz) \
X(chi_rhs) X(trK_rhs) X(gxx_rhs) X(gxy_rhs) X(gxz_rhs) X(gyy_rhs) X(gyz_rhs) X(gzz_rhs) \
X(Axx_rhs) X(Axy_rhs) X(Axz_rhs) X(Ayy_rhs) X(Ayz_rhs) X(Azz_rhs) \
X(Gamx_rhs) X(Gamy_rhs) X(Gamz_rhs) X(Lap_rhs) X(betax_rhs) X(betay_rhs) X(betaz_rhs) \
X(dtSfx_rhs) X(dtSfy_rhs) X(dtSfz_rhs) X(rho) X(Sx) X(Sy) X(Sz) \
X(Sxx) X(Sxy) X(Sxz) X(Syy) X(Syz) X(Szz) \
X(Gamxxx) X(Gamxxy) X(Gamxxz) X(Gamxyy) X(Gamxyz) X(Gamxzz) \
X(Gamyxx) X(Gamyxy) X(Gamyxz) X(Gamyyy) X(Gamyyz) X(Gamyzz) \
X(Gamzxx) X(Gamzxy) X(Gamzxz) X(Gamzyy) X(Gamzyz) X(Gamzzz) \
X(Rxx) X(Rxy) X(Rxz) X(Ryy) X(Ryz) X(Rzz) X(ham_Res) X(movx_Res) X(movy_Res) X(movz_Res) X(Gmx_Res) X(Gmy_Res) X(Gmz_Res)


struct DiffItem {
    const char *name;
    double diff;
};

struct Fields {
#define DECL(name) std::vector<double> name;
    ALL_FIELDS(DECL)
#undef DECL

    explicit Fields(size_t n) {
#define INIT(name) name.assign(n, 0.0);
        ALL_FIELDS(INIT)
#undef INIT
    }
};

static void fill_grid(std::vector<double> &g, double h) {
    for (size_t i = 0; i < g.size(); ++i) g[i] = h * static_cast<double>(i);
}

static void fill_fields(Fields &f) {
#define FILL(name) for (size_t i=0;i<f.name.size();++i) f.name[i] = std::sin(0.01*i + 0.1) + std::cos(0.02*i + 0.3);
    ALL_FIELDS(FILL)
#undef FILL
}

static double max_abs_diff(const std::vector<double> &a, const std::vector<double> &b) {
    double m = 0.0;
    for (size_t i = 0; i < a.size(); ++i) m = std::max(m, std::fabs(a[i] - b[i]));
    return m;
}

int main() {
    int ex[3] = {8, 7, 6};
    const size_t n = static_cast<size_t>(ex[0]) * ex[1] * ex[2];

    std::vector<double> X(ex[0]), Y(ex[1]), Z(ex[2]);
    fill_grid(X, 0.1); fill_grid(Y, 0.11); fill_grid(Z, 0.12);

    Fields c(n), f(n);
    fill_fields(c);
    f = c;

    int Symmetry = 0, Lev = 1, co = 0;
    double eps = 0.01, T = 0.0;

    int gont_c = f_compute_rhs_bssn(ex, T,
        X.data(), Y.data(), Z.data(),
        c.chi.data(), c.trK.data(),
        c.dxx.data(), c.gxy.data(), c.gxz.data(), c.dyy.data(), c.gyz.data(), c.dzz.data(),
        c.Axx.data(), c.Axy.data(), c.Axz.data(), c.Ayy.data(), c.Ayz.data(), c.Azz.data(),
        c.Gamx.data(), c.Gamy.data(), c.Gamz.data(),
        c.Lap.data(), c.betax.data(), c.betay.data(), c.betaz.data(),
        c.dtSfx.data(), c.dtSfy.data(), c.dtSfz.data(),
        c.chi_rhs.data(), c.trK_rhs.data(),
        c.gxx_rhs.data(), c.gxy_rhs.data(), c.gxz_rhs.data(), c.gyy_rhs.data(), c.gyz_rhs.data(), c.gzz_rhs.data(),
        c.Axx_rhs.data(), c.Axy_rhs.data(), c.Axz_rhs.data(), c.Ayy_rhs.data(), c.Ayz_rhs.data(), c.Azz_rhs.data(),
        c.Gamx_rhs.data(), c.Gamy_rhs.data(), c.Gamz_rhs.data(),
        c.Lap_rhs.data(), c.betax_rhs.data(), c.betay_rhs.data(), c.betaz_rhs.data(),
        c.dtSfx_rhs.data(), c.dtSfy_rhs.data(), c.dtSfz_rhs.data(),
        c.rho.data(), c.Sx.data(), c.Sy.data(), c.Sz.data(),
        c.Sxx.data(), c.Sxy.data(), c.Sxz.data(), c.Syy.data(), c.Syz.data(), c.Szz.data(),
        c.Gamxxx.data(), c.Gamxxy.data(), c.Gamxxz.data(), c.Gamxyy.data(), c.Gamxyz.data(), c.Gamxzz.data(),
        c.Gamyxx.data(), c.Gamyxy.data(), c.Gamyxz.data(), c.Gamyyy.data(), c.Gamyyz.data(), c.Gamyzz.data(),
        c.Gamzxx.data(), c.Gamzxy.data(), c.Gamzxz.data(), c.Gamzyy.data(), c.Gamzyz.data(), c.Gamzzz.data(),
        c.Rxx.data(), c.Rxy.data(), c.Rxz.data(), c.Ryy.data(), c.Ryz.data(), c.Rzz.data(),
        c.ham_Res.data(), c.movx_Res.data(), c.movy_Res.data(), c.movz_Res.data(),
        c.Gmx_Res.data(), c.Gmy_Res.data(), c.Gmz_Res.data(),
        Symmetry, Lev, eps, co);

    int gont_f = compute_rhs_bssn_(ex, &T,
        X.data(), Y.data(), Z.data(),
        f.chi.data(), f.trK.data(),
        f.dxx.data(), f.gxy.data(), f.gxz.data(), f.dyy.data(), f.gyz.data(), f.dzz.data(),
        f.Axx.data(), f.Axy.data(), f.Axz.data(), f.Ayy.data(), f.Ayz.data(), f.Azz.data(),
        f.Gamx.data(), f.Gamy.data(), f.Gamz.data(),
        f.Lap.data(), f.betax.data(), f.betay.data(), f.betaz.data(),
        f.dtSfx.data(), f.dtSfy.data(), f.dtSfz.data(),
        f.chi_rhs.data(), f.trK_rhs.data(),
        f.gxx_rhs.data(), f.gxy_rhs.data(), f.gxz_rhs.data(), f.gyy_rhs.data(), f.gyz_rhs.data(), f.gzz_rhs.data(),
        f.Axx_rhs.data(), f.Axy_rhs.data(), f.Axz_rhs.data(), f.Ayy_rhs.data(), f.Ayz_rhs.data(), f.Azz_rhs.data(),
        f.Gamx_rhs.data(), f.Gamy_rhs.data(), f.Gamz_rhs.data(),
        f.Lap_rhs.data(), f.betax_rhs.data(), f.betay_rhs.data(), f.betaz_rhs.data(),
        f.dtSfx_rhs.data(), f.dtSfy_rhs.data(), f.dtSfz_rhs.data(),
        f.rho.data(), f.Sx.data(), f.Sy.data(), f.Sz.data(),
        f.Sxx.data(), f.Sxy.data(), f.Sxz.data(), f.Syy.data(), f.Syz.data(), f.Szz.data(),
        f.Gamxxx.data(), f.Gamxxy.data(), f.Gamxxz.data(), f.Gamxyy.data(), f.Gamxyz.data(), f.Gamxzz.data(),
        f.Gamyxx.data(), f.Gamyxy.data(), f.Gamyxz.data(), f.Gamyyy.data(), f.Gamyyz.data(), f.Gamyzz.data(),
        f.Gamzxx.data(), f.Gamzxy.data(), f.Gamzxz.data(), f.Gamzyy.data(), f.Gamzyz.data(), f.Gamzzz.data(),
        f.Rxx.data(), f.Rxy.data(), f.Rxz.data(), f.Ryy.data(), f.Ryz.data(), f.Rzz.data(),
        f.ham_Res.data(), f.movx_Res.data(), f.movy_Res.data(), f.movz_Res.data(),
        f.Gmx_Res.data(), f.Gmy_Res.data(), f.Gmz_Res.data(),
        &Symmetry, &Lev, &eps, &co);

    std::vector<DiffItem> diffs;
    diffs.reserve(128);
    double maxd = 0.0;
#define CMP(name) { double d = max_abs_diff(c.name, f.name); maxd = std::max(maxd, d); diffs.push_back({#name, d}); }
    ALL_FIELDS(CMP)
#undef CMP

    std::sort(diffs.begin(), diffs.end(), [](const DiffItem &a, const DiffItem &b) { return a.diff > b.diff; });

    std::printf("compute_rhs_bssn gont: C=%d Fortran=%d\n", gont_c, gont_f);
    std::printf("compute_rhs_bssn max_abs_diff = %.3e\n", maxd);
    std::printf("Top 12 field diffs:\n");
    for (size_t i = 0; i < std::min<size_t>(12, diffs.size()); ++i) {
        std::printf("  %2zu) %-10s %.6e\n", i + 1, diffs[i].name, diffs[i].diff);
    }
    return (gont_c == gont_f && maxd < 1e-10) ? 0 : 1;
}
