# PowerChina-1 (Spacety) DInSAR Processing System

本项目是电建一号（涪城一号 / 天仪 BC3）X波段雷达卫星的高精度 DInSAR 数据自动化处理系统（V2.3）。基于深度定制的 GMTSAR，旨在实现从原始数据到高精度形变监测结果的一键式自动化处理。

## 核心特性
- **支持硬件及数据:** 针对电建一号 (DJ1) 与天仪系列 (BC3/LT1) X 波段数据定制优化。
- **自动化流水线 (`run_dinsar.sh`):** 全流程自动化处理，包含影像对齐、干涉图生成、自适应 Goldstein 滤波与解缠控制。
- **APS 大气相位校正系统:** 支持三种修正模式 (`--atm <none|empirical|gacos|st_filter>`) 以减少大气引起的相位延迟。
- **高稳健性:** 支持 snaphu OOM 自动分块重试和意外中断解缠恢复，具备生成 ArcGIS 兼容的 GeoTIFF 及自动地理编码能力。

## 目录结构说明
- `gmtsar-master/`: 定制版 GMTSAR C 源码及 CSH 核心处理脚本。（需编译使用）
- `deploy/`: 系统离线与在线安装所需的打包方案及一键部署脚本。
- `scripts/`: 系统顶层调度脚本库，包括流水线执行及自检脚本。
- `docs/`: 系统的内部使用及部署参考手册。

## 安装指南
详情请参考 `docs/电建一号 DInSAR 处理软件部署与操作手册.md` 以及 `docs/DEPLOYMENT_GUIDE.md`。

## 版权及使用说明
（内置各子工具版权归原作者，本项目框架代码用于特定的处理系统配置及运行）
