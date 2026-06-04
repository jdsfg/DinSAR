# 上游整合说明 (DONGYUSEN/gmtsar)

**整合日期:** 2026-01-30  
**上游仓库:** https://github.com/DONGYUSEN/gmtsar  
**整合范围:** 2026-01-21 及之后在 `origin/master` 上提交并更新的代码

## 整合的提交 (7 个)

| 提交 | 日期 | 说明 |
|------|------|------|
| 61a2382 | 2026-01-30 | Add files via upload – `preproc/LT1_preproc/src/make_slc_lt1.c` |
| 94a369c | 2026-01-30 | Add files via upload – `gmtsar/fit_coefficients.c`, `gmtsar/resamp_lt1.c` |
| af964a6 | 2026-01-26 | Add files via upload – `preproc/LT1_preproc/src/make_slc_lt1.c` |
| e915d5d | 2026-01-25 | Refactor xcorr2.c for improved performance and clarity |
| dd1dffa | 2026-01-24 | Add files via upload – `gmtsar/SAT_llt2rat2.c` |
| c965527 | 2026-01-22 | Add files via upload – `gmtsar/csh/offset_2_plot.sh` |
| 53885f2 | 2026-01-21 | Add files via upload – `gmtsar/csh/grd2tiff.csh` |

## 涉及文件 (7 个路径)

- **已覆盖（本地原有）：**  
  `preproc/LT1_preproc/src/make_slc_lt1.c`  
  `gmtsar/xcorr2.c`  
  `gmtsar/SAT_llt2rat2.c`  
  `gmtsar/csh/grd2tiff.csh`

- **新增：**  
  `gmtsar/fit_coefficients.c`  
  `gmtsar/resamp_lt1.c`  
  `gmtsar/csh/offset_2_plot.sh`

## 备份

整合前被覆盖的 4 个文件已备份到：

- `.integration_backup_20260130/preproc/LT1_preproc/src/make_slc_lt1.c`
- `.integration_backup_20260130/gmtsar/xcorr2.c`
- `.integration_backup_20260130/gmtsar/SAT_llt2rat2.c`
- `.integration_backup_20260130/gmtsar/csh/grd2tiff.csh`

如需恢复本地版本，可从上述目录复制回对应路径。

## 构建说明

- `fit_coefficients.c`、`resamp_lt1.c` 当前未加入顶层 `gmtsar/Makefile` 的 `PROGS_C` 或 `LIB_C`；若上游后续将其纳入构建，需再同步 Makefile。
- 其余已整合文件（如 `xcorr2.c`、`SAT_llt2rat2.c`）已在现有 Makefile 中，整合后直接 `make` 即可参与编译。
