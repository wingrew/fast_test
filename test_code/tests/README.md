# C/Fortran 数值一致性测试

## 依赖
- `gcc`
- `gfortran`（包含 `f951` 前端）

## 运行
```bash
./test_code/tests/run_c_fortran_tests.sh
```

脚本会分别编译并运行：
- `fderivs`
- `fdderivs`
- `kodis`
- `lopsided`
- `compute_rhs_bssn`（C: `f_compute_rhs_bssn` vs Fortran: `compute_rhs_bssn`）

若缺少 Fortran 前端，脚本会在启动时给出明确提示并退出。
