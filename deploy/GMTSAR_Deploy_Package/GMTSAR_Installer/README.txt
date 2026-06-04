PowerChina-1 DInSAR 离线部署包 v2.1
====================================

版本说明:
  v2.1 是合并更新版，包含：
  - GCC 11+ 兼容性修复 (-fcommon)
  - ALOS_preproc extern 全局变量修复
  - 完整的自动化数据处理脚本 run_dinsar.sh

部署步骤:
---------

1. 解压部署包
   tar xzf dinsar_deploy_v2.1_offline.tar.gz
   cd GMTSAR_Installer

2. 安装系统依赖 (首次需要联网或使用离线包)
   sudo apt update
   sudo apt install gmt gmt-gshhg libgmt-dev libnetcdf-dev \
       libfftw3-dev libgdal-dev build-essential autoconf gfortran

3. 运行部署脚本
   ./install.sh

4. 验证安装
   source ~/.bashrc
   make_slc_dj1 --help
   xcorr --help

数据处理:
---------

1. 准备数据目录结构:
   project/
   ├── raw/              # 放置 XML + TIFF 文件
   │   ├── image1.xml
   │   ├── image1.tiff
   │   ├── image2.xml
   │   └── image2.tiff
   └── topo/             # 放置 DEM (可选)
       └── dem.grd

2. 运行处理脚本:
   ./run_dinsar.sh /path/to/project              # 交互模式
   ./run_dinsar.sh --auto --mode A /path/to/project  # 全自动

3. 查看结果:
   project/results/
   ├── interferogram.png      # 干涉条纹图
   ├── displacement.png       # 位移图
   ├── displacement.tif       # GeoTIFF (可导入 GIS)
   └── los_displacement.grd   # 位移栅格

解缠模式说明:
-------------
Mode A (稳定): 2x 降采样 + 分块解缠，推荐 16GB 内存
Mode B (高精度): 全分辨率整图解缠，需要 32GB+ 内存

支持的卫星:
-----------
- DJ1 (电建一号) - PowerChina-1
- BC3 (涪城一号) - Spacety/TianYi
- 波长: 0.0311m (X波段, 9.65 GHz)

故障排除:
---------
- 编译失败 "multiple definition": 脚本会自动修复
- snaphu 内存不足: 使用 Mode A 或增加降采样
- 日志位置: project/logs/processing_*.log

版本: v2.1 (2026-01-28)

