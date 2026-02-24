#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD_DIR="$ROOT/tests/build"
mkdir -p "$BUILD_DIR"

if ! command -v gfortran >/dev/null 2>&1; then
  echo "ERROR: gfortran not found. Please install gfortran (provides f951) before running this script." >&2
  exit 2
fi

if ! gfortran -v >/dev/null 2>&1; then
  echo "ERROR: gfortran is present but not functional (f951 frontend may be missing)." >&2
  exit 2
fi

compile_and_run () {
  local name="$1"
  local macro="$2"
  local c_src="$3"
  local f_src="$4"
  local exe="$BUILD_DIR/$name"

  gcc -std=c11 -O2 -D"$macro" "$ROOT/tests/test_c_vs_fortran.c" "$ROOT/c/$c_src"       "$ROOT/fortran/$f_src" "$ROOT/fortran/symmetry_bd.f" -lgfortran -lm -o "$exe"
  "$exe"
}

compile_and_run_bssn () {
  local exe="$BUILD_DIR/bssn"
  local cxxflags="-std=c++17 -O2"

  g++ $cxxflags -c "$ROOT/tests/test_bssn_c_vs_fortran.cpp" -o "$BUILD_DIR/test_bssn_c_vs_fortran.o"
  g++ $cxxflags -c "$ROOT/c/bssn_rhs.c" -o "$BUILD_DIR/bssn_rhs.o"
  g++ $cxxflags -c "$ROOT/c/fderivs.c" -o "$BUILD_DIR/fderivs_cpp.o"
  g++ $cxxflags -c "$ROOT/c/fdderivs.c" -o "$BUILD_DIR/fdderivs_cpp.o"
  g++ $cxxflags -c "$ROOT/c/kodiss.c" -o "$BUILD_DIR/kodiss_cpp.o"
  g++ $cxxflags -c "$ROOT/c/lopsided.c" -o "$BUILD_DIR/lopsided_cpp.o"

  gfortran -O2 -cpp -c "$ROOT/fortran/compute_rhs_bssn.f" -o "$BUILD_DIR/compute_rhs_bssn.o"
  gfortran -O2 -cpp -c "$ROOT/fortran/fderivs.f" -o "$BUILD_DIR/fderivs_f.o"
  gfortran -O2 -cpp -c "$ROOT/fortran/fdderivs.f" -o "$BUILD_DIR/fdderivs_f.o"
  gfortran -O2 -cpp -c "$ROOT/fortran/kodiss.f" -o "$BUILD_DIR/kodiss_f.o"
  gfortran -O2 -cpp -c "$ROOT/fortran/lopsided.f" -o "$BUILD_DIR/lopsided_f.o"
  gfortran -O2 -cpp -c "$ROOT/fortran/symmetry_bd.f" -o "$BUILD_DIR/symmetry_bd_f.o"

  g++ $cxxflags     "$BUILD_DIR/test_bssn_c_vs_fortran.o"     "$BUILD_DIR/bssn_rhs.o" "$BUILD_DIR/fderivs_cpp.o" "$BUILD_DIR/fdderivs_cpp.o" "$BUILD_DIR/kodiss_cpp.o" "$BUILD_DIR/lopsided_cpp.o"     "$BUILD_DIR/compute_rhs_bssn.o" "$BUILD_DIR/fderivs_f.o" "$BUILD_DIR/fdderivs_f.o" "$BUILD_DIR/kodiss_f.o" "$BUILD_DIR/lopsided_f.o" "$BUILD_DIR/symmetry_bd_f.o"     -lgfortran -lm -o "$exe"
  "$exe"
}

compile_and_run fderivs TEST_FDERIVS fderivs.c fderivs.f
compile_and_run fdderivs TEST_FDDERIVS fdderivs.c fdderivs.f
compile_and_run kodis TEST_KODIS kodiss.c kodiss.f
compile_and_run lopsided TEST_LOPSIDED lopsided.c lopsided.f

compile_and_run_bssn

echo "All C vs Fortran comparisons passed."
