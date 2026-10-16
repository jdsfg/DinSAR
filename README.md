# PowerChina-1 (Spacety) DInSAR Processing System

本项目是电建一号（涪城一号 / 天仪 BC3）X波段雷达卫星的高精度 DInSAR 数据自动化处理系统（V2.3）。基于深度定制的 GMTSAR，旨在实现从原始数据到高精度形变监测结果的一键式自动化处理。

## 核心特性
- **支持硬件及数据:** 针对电建一号 (DJ1) 与天仪系列 (BC3/LT1) X 波段数据定制优化。
- **自动化流水线 (`run_dinsar.sh`):** 全流程自动化处理，包含影像对齐、干涉图生成、自适应 Goldstein 滤波与解缠控制。
- **APS 大气相位校正系统:** 支持三种修正模式 (`--atm <none|empirical|gacos|st_filter>`) 以减少大气引起的相位延迟。
- **高稳健性:** 支持 snaphu OOM 自动分块重试和意外中断解缠恢复，具备生成 ArcGIS 兼容的 GeoTIFF 及自动地理编码能力。

## 新版电建一号数据 (XML + HDF5) 与多尺度配准
新版产品的 XML 结构和 SLC 格式（h5 单精度复数）与旧版不同，旧的 `make_slc_dj1` 不能处理。相关工具（都带单元测试）：
- `gmtsar-master/preproc/DJ1_preproc/src_swath/dj1_pair_select.py ROOT [--json pairs.json]`：只读 XML 的像对筛选。
- `make_slc_dj1h5.py PRODUCT.XML PRODUCT.h5 20260904`（需要 numpy 和 h5py）：转成 GMTSAR 的 `.SLC/.PRM/.LED`，并写转换说明 `.convert.txt`。实时轨道只覆盖成像时段，转换时会外推到场景前后；PRM 的 `fd1` 固定为 0（零多普勒产品，非零值会让 `SAT_llt2rat` 的方位坐标错开上千行）。PRM 里没有 `SC_vel/earth_radius/SC_height`，按 GMTSAR 惯例用 `calc_dop_orb` 追加（从影像沿用主影像的 `earth_radius`）。
- `gmtsar-master/gmtsar/csh/ms_register.py register M.PRM S.PRM ...`：多尺度配准（不依赖 PRM 里的初始偏移；估计、`resamp` 重采样、独立点位复核一步完成，未通过复核时返回非 0 且不留下半成品）。
`register` 默认使用原生 C 估计器；可用 `--engine python` 选择保留的 Python 参考实现。构建 C 程序：`make -C gmtsar-master/gmtsar ms_register`。也可通过 `MS_REGISTER_BIN` 指定可执行文件。
- `run_dinsar.sh --coreg multiscale`：在流水线里启用多尺度配准，默认 `legacy`（`xcorr + fitoffset`，行为不变）。复核未通过会直接中止，不会悄悄退回 legacy。

`run_dinsar.sh` 现在可以直接处理新版 XML + HDF5：`raw/` 里放 XML 和同名 `.h5`（同名 `.tiff` 只是 8 bit 预览图，有 `.h5` 时不用它），其余步骤自动完成——转换 SLC、用 `calc_dop_orb` 补轨道参数、默认 multiscale 配准、DEM 转雷达坐标、去地形相位、滤波、解缠、去趋势、地理编码。DEM 放 `topo/dem.grd` 或 `topo/dem.tif`（高程须是相对 WGS84 椭球的高，SRTM/Copernicus 的正高要先加大地水准面起伏）。
- DEM 转雷达坐标由 `topo_ra_grid.py`（SciPy Delaunay 线性内插，结果与 `gmt triangulate` 一致）完成，输出 SLC 的 1/4 分辨率网格（`phasediff` 支持比 SLC 粗的 topo 网格）；整景分辨率的网格与 `real/imag` 同样大（SM 整景 1.45 GB），在 4 GB 内存的机器上放不下。
- 实测：SM 像对 20260904→20260912（24464×14836 像元，一对 2.9 GB SLC）从 XML+HDF5 到 GeoTIFF 约 20～21 分钟（两次整景运行分别为 20 分 56 秒和 20 分 3 秒，同一输入的解缠结果逐字节一致）。整个进程放在 4.4 GB 的 cgroup 上限里运行，未触发 OOM：匿名内存峰值约 2.9 GB（出现在 `phasediff` 和随后的整景相位/振幅计算阶段，约 8 字节 × 像元数，即两个整景网格 `real`/`imag` 的大小），cgroup 总量最高 4.2 GB，其余是可回收的文件缓存。其他产品的内存需求可按这个比例估计（例如像元数约 2 倍的 SP 像对约需 5.5 GB 匿名内存，尚未实测）。
- 新版数据的默认值（环境变量可覆盖，旧版数据的默认值不变）：`DETREND_ORDER=3`（去平面；`1` 只减均值，不去倾斜）、`DISP_STRETCH_PCT=98`（位移图色标取 |位移| 的 98 分位；`100` 为最大绝对值）、`TOPO_RA_DEC=auto`（粗网格倍数；`1` 为整景分辨率）。新版数据的轨道常是实时轨道，X 波段干涉图里会留下十几周的残余倾斜（实测 SM 像对沿方位向 11.7 周），只减均值时位移图会被它主导。
- 判读注意：`corr.tif` 是 `|<I>|/<|I|>`（干涉图的相位相干度），纯噪声区的底约 0.3，不是标准相干系数（标准相干系数的整景均值约 0.46）。去平面后位移的稳健标准差约 1 cm，主要是大气，只有一对影像时无法把形变和大气分开。
- 两景的初始偏移（SM 约 406 个距离向像元）会让干涉图一侧约 1.7% 宽的条带没有重叠（SM 在最右侧）。这条带里的相位和位移只是噪声，`corr` 也不是 0（约 0.3），请裁掉或忽略。
- 绝对定位：用 DEM 模拟亮度与 SAR 幅度做互相关（山区块的相关系数 0.7～0.9，各块给出的偏移一致），DEM 与影像整体相差约 20 m（方位向）和 5 m（距离向），主要来自实时轨道和对流层延迟。SM 像对的垂直基线只有约 26 m，地形相位总共约 3 周，这个偏移对干涉相位几乎没有影响，但地理编码结果的绝对位置约有 20 m 的不确定度。

测试：`python3 -m unittest discover -s tests`（共 90 项，约 6 分钟，大头是配准的合成数据测试）。缺少 h5py、scipy，或 GMTSAR 的 `resamp`、`calc_dop_orb`、`SAT_llt2rat` 时，对应的测试会被跳过；跳过不等于通过，请看输出里的 `skipped`。

## Ray-traced atmospheric correction

For `run_dinsar.sh --atm raytrace`, provide normalized `master_weather.nc` and
`slave_weather.nc` cubes in the configured weather directory, `stage/weather`,
`stage`, or `project/weather`. Each cube's `epoch_utc` must be within 30 minutes
of its PRM acquisition time. Alternatively, place an `*era5*raw*.npz` archive
in one of those locations; the workflow builds the two cubes from the PRM epochs
and the geographic DEM bounds expanded by 0.75 degrees.

Fetch the hourly ARCO-ERA5 inputs with:

```sh
python3 gmtsar-master/gmtsar/csh/era5_weather.py fetch --bbox 107.25/112/39.25/42.25 \
  --epoch 2026-09-04T04:23:44Z --epoch 2026-09-12T04:23:44Z --out era5_raw.npz
```

The fetch bbox must cover the DEM bounds plus the 0.75 degree margin (the example is the DJ1 scene). The fetch needs `pip install xarray "zarr<3" gcsfs` and internet access; it can run on another machine and only the npz needs to be copied. ERA5 dates in 2026 are ERA5T preliminary data. The workflow's physical phase
correction subtracts the plane-removed atmospheric phase; its `--scale` option
defaults to `1.0`.

## 目录结构说明
- `gmtsar-master/`: 定制版 GMTSAR C 源码及 CSH 核心处理脚本。（需编译使用）
- `deploy/`: 系统离线与在线安装所需的打包方案及一键部署脚本。
- `scripts/`: 系统顶层调度脚本库，包括流水线执行及自检脚本。
- `docs/`: 系统的内部使用及部署参考手册。

## 安装指南
详情请参考 `docs/电建一号 DInSAR 处理软件部署与操作手册.md` 以及 `docs/DEPLOYMENT_GUIDE.md`。

## 版权及使用说明
（内置各子工具版权归原作者，本项目框架代码用于特定的处理系统配置及运行）
