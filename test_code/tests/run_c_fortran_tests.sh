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

compile_and_run fderivs TEST_FDERIVS fderivs.c fderivs.f
compile_and_run fdderivs TEST_FDDERIVS fdderivs.c fdderivs.f
compile_and_run kodis TEST_KODIS kodiss.c kodiss.f
compile_and_run lopsided TEST_LOPSIDED lopsided.c lopsided.f

echo "All C vs Fortran comparisons passed."
