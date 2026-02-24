#ifndef SHARE_FUNC_H
#define SHARE_FUNC_H

#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <stdio.h>

/* 主网格：0-based -> 1D */
static inline size_t idx_ex(int i0, int j0, int k0, const int ex[3]) {
    const int ex1 = ex[0], ex2 = ex[1];
    return (size_t)i0 + (size_t)j0 * (size_t)ex1 + (size_t)k0 * (size_t)ex1 * (size_t)ex2;
}

/*
 * fh 对应 Fortran: fh(-1:ex1, -1:ex2, -1:ex3)
 * ord=2 => shift=1
 * iF/jF/kF 为 Fortran 索引（可为 -1,0,1..ex）
 */
static inline size_t idx_fh_F_ord2(int iF, int jF, int kF, const int ex[3]) {
    const int shift = 1;
    const int nx = ex[0] + 2;      // ex1 + ord
    const int ny = ex[1] + 2;

    const int ii = iF + shift;     // 0..ex1+1
    const int jj = jF + shift;     // 0..ex2+1
    const int kk = kF + shift;     // 0..ex3+1

    return (size_t)ii + (size_t)jj * (size_t)nx + (size_t)kk * (size_t)nx * (size_t)ny;
}

/*
 * fh 对应 Fortran: fh(-2:ex1, -2:ex2, -2:ex3)
 * ord=3 => shift=2
 * iF/jF/kF 是 Fortran 索引（可为负）
 */
static inline size_t idx_fh_F(int iF, int jF, int kF, const int ex[3]) {
    const int shift = 2;                 // ord=3 -> -2..ex
    const int nx = ex[0] + 3;            // ex1 + ord
    const int ny = ex[1] + 3;

    const int ii = iF + shift;           // 0..ex1+2
    const int jj = jF + shift;           // 0..ex2+2
    const int kk = kF + shift;           // 0..ex3+2

    return (size_t)ii + (size_t)jj * (size_t)nx + (size_t)kk * (size_t)nx * (size_t)ny;
}

/*
 * func:  (1..extc1, 1..extc2, 1..extc3)   1-based in Fortran
 * funcc: (-ord+1..extc1, -ord+1..extc2, -ord+1..extc3) in Fortran
 *
 * C 里我们把：
 *   func  视为 0-based: i0=0..extc1-1, j0=0..extc2-1, k0=0..extc3-1
 *   funcc 用“平移下标”存为一维数组：
 *     iF in [-ord+1..extc1]  -> ii = iF + (ord-1)  in [0..extc1+ord-1]
 *     总长度 nx = extc1 + ord
 *     同理 ny = extc2 + ord, nz = extc3 + ord
 */

static inline size_t idx_func0(int i0, int j0, int k0, const int extc[3]) {
    const int nx = extc[0], ny = extc[1];
    return (size_t)i0 + (size_t)j0 * (size_t)nx + (size_t)k0 * (size_t)nx * (size_t)ny;
}

static inline size_t idx_funcc_F(int iF, int jF, int kF, int ord, const int extc[3]) {
    const int shift = ord - 1;          // iF = -shift .. extc1
    const int nx = extc[0] + ord;       // [-shift..extc1] 共 extc1+ord 个
    const int ny = extc[1] + ord;

    const int ii = iF + shift;          // 0..extc1+shift
    const int jj = jF + shift;          // 0..extc2+shift
    const int kk = kF + shift;          // 0..extc3+shift

    return (size_t)ii + (size_t)jj * (size_t)nx + (size_t)kk * (size_t)nx * (size_t)ny;
}

/*
 * 等价于 Fortran:
 * funcc(1:extc1,1:extc2,1:extc3)=func
 * do i=0,ord-1
 *   funcc(-i,1:extc2,1:extc3) = funcc(i+1,1:extc2,1:extc3)*SoA(1)
 * enddo
 * do i=0,ord-1
 *   funcc(:,-i,1:extc3) = funcc(:,i+1,1:extc3)*SoA(2)
 * enddo
 * do i=0,ord-1
 *   funcc(:,:,-i) = funcc(:,:,i+1)*SoA(3)
 * enddo
 */
static inline void symmetry_bd(int ord,
                 const int extc[3],
                 const double *func,
                 double *funcc,
                 const double SoA[3])
{
    const int extc1 = extc[0], extc2 = extc[1], extc3 = extc[2];

    // 1) funcc(1:extc1,1:extc2,1:extc3) = func
    // Fortran 的 (iF=1..extc1) 对应 C 的 func(i0=0..extc1-1)
    for (int k0 = 0; k0 < extc3; ++k0) {
        for (int j0 = 0; j0 < extc2; ++j0) {
            for (int i0 = 0; i0 < extc1; ++i0) {
                const int iF = i0 + 1, jF = j0 + 1, kF = k0 + 1;
                funcc[idx_funcc_F(iF, jF, kF, ord, extc)] = func[idx_func0(i0, j0, k0, extc)];
            }
        }
    }

    // 2) do i=0..ord-1: funcc(-i, 1:extc2, 1:extc3) = funcc(i+1, ...)*SoA(1)
    for (int ii = 0; ii <= ord - 1; ++ii) {
        const int iF_dst = -ii;       // 0, -1, -2, ...
        const int iF_src = ii + 1;    // 1, 2, 3, ...
        for (int kF = 1; kF <= extc3; ++kF) {
            for (int jF = 1; jF <= extc2; ++jF) {
                funcc[idx_funcc_F(iF_dst, jF, kF, ord, extc)] =
                    funcc[idx_funcc_F(iF_src, jF, kF, ord, extc)] * SoA[0];
            }
        }
    }

    // 3) do i=0..ord-1: funcc(:,-i, 1:extc3) = funcc(:, i+1, 1:extc3)*SoA(2)
    // 注意 Fortran 这里的 ":" 表示 iF 从 (-ord+1..extc1) 全覆盖
    for (int jj = 0; jj <= ord - 1; ++jj) {
        const int jF_dst = -jj;
        const int jF_src = jj + 1;
        for (int kF = 1; kF <= extc3; ++kF) {
            for (int iF = -ord + 1; iF <= extc1; ++iF) {
                funcc[idx_funcc_F(iF, jF_dst, kF, ord, extc)] =
                    funcc[idx_funcc_F(iF, jF_src, kF, ord, extc)] * SoA[1];
            }
        }
    }

    // 4) do i=0..ord-1: funcc(:,:,-i) = funcc(:,:, i+1)*SoA(3)
    for (int kk = 0; kk <= ord - 1; ++kk) {
        const int kF_dst = -kk;
        const int kF_src = kk + 1;
        for (int jF = -ord + 1; jF <= extc2; ++jF) {
            for (int iF = -ord + 1; iF <= extc1; ++iF) {
                funcc[idx_funcc_F(iF, jF, kF_dst, ord, extc)] =
                    funcc[idx_funcc_F(iF, jF, kF_src, ord, extc)] * SoA[2];
            }
        }
    }
}
#endif
