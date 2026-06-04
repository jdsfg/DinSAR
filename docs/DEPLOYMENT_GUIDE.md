# 涪城一号 DInSAR 算法模型部署文档

## 文档信息

| 项目 | 说明 |
|------|------|
| **项目名称** | PowerChina-1 (Spacety) DInSAR 系统 |
| **版本** | v3.1 |
| **文档日期** | 2026-01-19 |
| **适用卫星** | 天仪BC3 (百川三号) / 涪城一号 / DJ1 (电建一号) |
| **处理模式** | C波段 StripMap DInSAR |
| **网络要求** | 支持国内网络环境（自动切换阿里云镜像源） |

---

## 1. 开发环境要求

### 1.1 操作系统

| 组件 | 版本 | 说明 |
|------|------|------|
| **操作系统** | Ubuntu 22.04.5 LTS (Jammy Jellyfish) | 64位 |
| **内核版本** | 6.8.0-87-generic | x86_64 |
| **架构** | x86_64 | AMD64兼容 |

### 1.2 编译工具链

| 工具 | 版本 | 安装命令 |
|------|------|----------|
| **GCC** | 11.4.0 | `sudo apt install build-essential` |
| **GNU Make** | 4.3 | `sudo apt install make` |
| **Autoconf** | 2.71 | `sudo apt install autoconf` |
| **CMake** | 3.22+ | `sudo apt install cmake` |

---

## 2. 依赖库清单

### 2.1 核心依赖库

| 库名称 | 版本 | 用途 | 安装命令 |
|--------|------|------|----------|
| **GMT** | 6.3.0 | 地理数据处理与可视化 | `sudo apt install gmt gmt-dev libgmt-dev` |
| **NetCDF** | 4.8.1 | 科学数据格式支持 (.grd) | `sudo apt install libnetcdf-dev` |
| **HDF5** | 1.10.7 | 高性能数据存储 | `sudo apt install libhdf5-dev` |
| **FFTW3** | 3.3.8 | 快速傅里叶变换 | `sudo apt install libfftw3-dev` |
| **LAPACK** | 3.10.0 | 线性代数计算 | `sudo apt install liblapack-dev` |
| **BLAS** | 3.10.0 | 基础线性代数运算 | `sudo apt install libblas-dev` |

### 2.2 脚本语言环境

| 语言 | 版本 | 用途 | 安装命令 |
|------|------|------|----------|
| **TCSH/CSH** | 6.21.00 | 主处理脚本执行 | `sudo apt install tcsh csh` |
| **Python3** | 3.10+ | SNAPHU 数据转换、数组处理 | `sudo apt install python3 python3-numpy` |
| **AWK** | GNU AWK | 文本处理 | `sudo apt install gawk` |

### 2.3 一键安装所有依赖

```bash
# Ubuntu 22.04 一键安装命令（国内网络可用）
sudo apt update
sudo apt install -y \
    build-essential autoconf cmake gfortran \
    libgmt-dev gmt gmt-dcw gmt-gshhg \
    libnetcdf-dev libhdf5-dev \
    libfftw3-dev liblapack-dev libblas-dev \
    libglib2.0-dev libtiff-dev \
    tcsh csh gawk \
    python3 python3-numpy \
    wget curl
```

**注意：** 如果无法访问官方源，部署脚本会自动切换到阿里云镜像源。

---

## 3. GMTSAR 安装配置

### 3.1 安装路径

```
/usr/local/gmtsar/
├── bin/                    # 可执行文件 (38个)
├── gmtsar/
│   ├── csh/               # CSH处理脚本 (87个)
│   └── *.c                # C源代码
├── preproc/               # 预处理模块
├── snaphu/                # 相位解缠模块
└── share/                 # 共享资源
```

### 3.2 一键部署（推荐）

**使用自动部署脚本（v3.1）：**

```bash
# 1. 解压源码包
cd ~
tar -xzf PowerChina1_GMTSAR_Source_v3.tar.gz

# 2. 运行一键部署
chmod +x deploy_dinsar.sh
./deploy_dinsar.sh

# 3. 加载环境变量
source ~/.bashrc

# 4. 验证安装
p2p_processing.csh
```

**脚本功能：**
- ✅ 自动检测网络环境
- ✅ 无法访问官方源时自动切换阿里云镜像
- ✅ 安装所有依赖（包括 GMT、NetCDF、FFTW3 等）
- ✅ 编译 GMTSAR
- ✅ 配置环境变量
- ✅ 验证安装

---

### 3.3 手动编译步骤

如果需要手动编译：

```bash
cd /usr/local/gmtsar

# 1. 清理旧编译
make clean

# 2. 生成配置脚本
autoconf

# 3. 配置编译选项
./configure --prefix=/usr/local/GMTSAR

# 4. 编译
make

# 5. 安装 (需要sudo权限)
sudo make install
```

### 3.3 环境变量配置

添加以下内容到 `~/.bashrc`:

```bash
# GMTSAR 环境配置
export GMTSAR_HOME="/usr/local/gmtsar"
export PATH="$GMTSAR_HOME/bin:$GMTSAR_HOME/gmtsar/csh:$PATH"

# GMT 配置
export GMT_SHAREDIR="/usr/share/gmt"
```

使配置生效:
```bash
source ~/.bashrc
```

### 3.4 验证安装

```bash
# 检查核心工具
which make_slc_s1a        # 应返回: /usr/local/gmtsar/bin/make_slc_s1a
which p2p_processing.csh  # 应返回脚本路径
gmt --version             # 应返回: 6.3.0
snaphu                    # 应返回: snaphu v2.0.7
```

---

## 4. 算法模块说明

### 4.1 核心可执行文件

| 可执行文件 | 功能 | 所属模块 |
|------------|------|----------|
| `make_slc_s1a` | 生成SLC影像（支持天仪/涪城一号） | 预处理 |
| `make_slc_dj1` | 电建一号数据预处理 | 预处理 |
| `xcorr` | 影像配准 | 配准 |
| `xcorr_fast` | **快速配准（26×加速）** | 配准 |
| `phasediff` | 干涉相位计算 | 干涉 |
| `phasefilt` | 相位滤波 | 滤波 |
| `SAT_llt2rat` | 坐标转换 | 地理编码 |
| `SAT_baseline` | 基线计算 | 几何 |
| `conv` | 卷积滤波 | 滤波 |
| `snaphu` | 相位解缠 | 解缠 |
| `proj_ra2ll.csh` | 雷达到地理坐标投影 | 地理编码 |

### 4.2 处理脚本 (CSH)

| 脚本 | 功能 | 输入 | 输出 |
|------|------|------|------|
| `p2p_processing.csh` | **端到端处理（支持DJ1/S1_STRIP）** | XML+TIFF+DEM | 干涉图+KML |
| `pre_proc.csh` | 预处理 | 原始数据 | SLC+PRM+LED |
| `dem2topo_ra.csh` | DEM处理 | dem.grd | topo_ra.grd |
| `intf.csh` | 干涉处理 | SLC对 | phase.grd |
| `filter.csh` | 滤波处理 | phase.grd | phasefilt.grd |
| `geocode.csh` | 地理编码 | *_ra.grd | *_ll.grd |

### 4.3 自定义扩展脚本

| 脚本 | 功能 | 状态 | 位置 |
|------|------|------|------|
| **snaphu.csh** | 相位解缠包装（使用Python处理NaN） | ✅ 已实现 | `gmtsar/csh/` |
| **Ph2los.csh** | 相位转LOS形变 ($d = \frac{\lambda}{4\pi}\phi$) | ✅ 已实现 | `gmtsar/csh/` |
| **Geo_all.csh** | 批量地理编码 | ✅ 已实现 | `gmtsar/csh/` |
| **Re_atm.csh** | 大气校正（ERA5/电离层） | ⏳ 占位符 | `gmtsar/csh/` |

### 4.4 关键代码修复说明

| 文件 | 修复内容 | 原因 |
|------|----------|------|
| `gmtsar/gmtsar.h` | 全局变量添加 `extern` 声明 | 修复 GCC 10+ 链接错误 |
| `gmtsar/sarglobal.h` | `verbose/debug` 改为 `extern` | 防止多重定义 |
| `config.mk` | 添加 `-fcommon` 编译标志 | GCC 10+ 兼容性 |
| `make_slc_s1a.c` | `sv[200]` → `sv[1000]` + NULL检查 | 防止天仪数据崩溃 |
| `make_s1a_tops_6par.c` | `sv[400]` → `sv[1000]` | 支持长时间跨度数据 |
| `p2p_processing.csh` | 添加 DJ1/S1_STRIP + `xcorr_fast` | 26×加速配准 |
| `snaphu.csh` | 使用Python处理NaN值 | 修复SNAPHU输入格式问题 |

---

## 5. 模块运行指南

### 5.1 完整处理流程

```
原始数据 (XML+TIFF) → 预处理 → 配准 → 干涉 → 滤波 → 解缠 → 形变反演 → 地理编码
```

### 5.2 方式一: 端到端自动处理

```bash
cd /your/project/directory

# 准备目录结构










mkdir -p raw topo intf

# 将原始数据放入 raw/ 目录
cp FC1_20231110.xml FC1_20231110.tiff raw/
cp FC1_20231121.xml FC1_20231121.tiff raw/
cp dem.grd topo/

# 创建配置文件 (复制模板)
cp $GMTSAR_HOME/config.S1_STRIP.txt config.s1a.txt

# 运行端到端处理
p2p_processing.csh S1_STRIP FC1_20231110 FC1_20231121 config.s1a.txt
```

### 5.3 方式二: 分步处理

#### Step 1: 预处理 (生成SLC)

```bash
cd raw/

# 主影像
make_slc_s1a FC1_20231110.xml FC1_20231110.tiff FC1_20231110

# 辅影像  
make_slc_s1a FC1_20231121.xml FC1_20231121.tiff FC1_20231121
```

**输出文件**:
- `FC1_20231110.SLC` - 单视复数影像
- `FC1_20231110.PRM` - 参数文件
- `FC1_20231110.LED` - 轨道文件

#### Step 2: DEM处理

```bash
cd ../topo/

# 生成雷达坐标系地形
dem2topo_ra.csh ../raw/FC1_20231110.PRM dem.grd
```

**输出文件**:
- `topo_ra.grd` - 雷达坐标系地形

#### Step 3: 配准

```bash
cd ../raw/

# 影像配准
xcorr FC1_20231110.PRM FC1_20231121.PRM -xsearch 256 -ysearch 128 -nx 30 -ny 50
```

#### Step 4: 干涉处理

```bash
cd ../intf/
mkdir 2023313_2023324
cd 2023313_2023324

# 链接文件
ln -s ../../raw/*.SLC .
ln -s ../../raw/*.PRM .
ln -s ../../raw/*.LED .
ln -s ../../topo/topo_ra.grd .

# 干涉
phasediff FC1_20231110.PRM FC1_20231121.PRM -topo topo_ra.grd
```

**输出文件**:
- `phase.grd` - 干涉相位
- `corr.grd` - 相干性图

#### Step 5: 滤波

```bash
# Goldstein滤波
filter.csh FC1_20231110.PRM FC1_20231121.PRM 500 1
```

**输出文件**:
- `phasefilt.grd` - 滤波后相位

#### Step 6: 相位解缠

```bash
# 使用自定义脚本
snaphu.csh phasefilt.grd corr.grd 0.1
```

**输出文件**:
- `unwrap.grd` - 解缠相位

#### Step 7: 相位转形变

```bash
# 相位转LOS形变 (C波段波长0.0555m)
Ph2los.csh unwrap.grd 0.0555
```

**输出文件**:
- `los_def.grd` - LOS形变量 (米)

#### Step 8: 地理编码

```bash
# 地理编码
Geo_all.csh los_def.grd trans.dat 4
```

**输出文件**:
- `los_def_ll.grd` - WGS84坐标形变图
- `los_def_ll.kml` - Google Earth叠加

---

## 6. 关键参数配置

### 6.1 天仪BC3卫星参数

| 参数 | 值 | 说明 |
|------|-----|------|
| `SC_identity` | 10 | 卫星标识 (天仪) |
| `rng_samp_rate` | 120000000 | 距离采样率 (Hz) |
| `radar_wavelength` | 0.0555 | 波长 (m, C波段) |
| `PRF` | 4105.09 | 脉冲重复频率 (Hz) |

### 6.2 处理参数建议

| 参数 | 推荐值 | 说明 |
|------|--------|------|
| `filter_wavelength` | 500 | Goldstein滤波波长 (m) |
| `corr_threshold` | 0.1 | 相干性阈值 |
| `dec_factor` | 1 | 降采样因子 |
| `geocode_resolution` | 4 | 地理编码分辨率 (arc-sec) |

### 6.3 配置文件模板 (config.s1a.txt)

```bash
# 处理阶段控制
proc_stage = 1          # 从预处理开始
skip_stage =            # 不跳过任何阶段

# 滤波参数
filter_wavelength = 500
dec_factor = 1

# 相干性
threshold_correlation = 0.1

# 地理编码
range_dec = 1
azimuth_dec = 1
```

---

## 7. 数据格式规范

### 7.1 输入数据格式

| 文件类型 | 格式 | 说明 |
|----------|------|------|
| 元数据 | XML | 天仪卫星产品元数据 |
| 影像 | GeoTIFF | 复数影像数据 |
| DEM | NetCDF (.grd) | WGS84坐标系 |

### 7.2 输出数据格式

| 文件类型 | 格式 | 说明 |
|----------|------|------|
| SLC | Binary (short complex) | 单视复数影像 |
| PRM | Text | GMTSAR参数文件 |
| LED | Text | 轨道状态向量 |
| GRD | NetCDF | GMT栅格数据 |
| KML | XML | Google Earth叠加 |

---

## 8. 故障排查

### 8.1 常见问题

| 问题 | 原因 | 解决方案 |
|------|------|----------|
| `make_slc_s1a` 段错误 | XML空指针 | 检查XML格式完整性 |
| `SAT_llt2rat` 负值范围 | precise参数 | 使用 `precise=0` |
| 干涉条纹异常 | 配准精度不足 | 增加搜索窗口参数 |
| 解缠失败 | 相干性过低 | 降低 `corr_threshold` |

### 8.2 日志检查

```bash
# 检查处理日志
tail -100 p2p_processing.log

# 检查GMT错误
gmt grdinfo phase.grd
```

---

## 9. 性能参考

### 9.1 测试环境配置

| 硬件 | 规格 |
|------|------|
| CPU | Intel Core (x86_64) |
| 内存 | 32GB+ 推荐 |
| 存储 | SSD推荐 (处理速度提升30%) |

### 9.2 处理时间参考

| 处理阶段 | 原始时间 | 优化后时间 | 说明 |
|----------|----------|------------|------|
| 预处理 | 5-10 分钟 | 5-10 分钟 | - |
| 配准（xcorr） | 25 分钟 | **1 分钟** | 使用 `xcorr_fast` (26×加速) |
| 干涉+滤波 | 5-10 分钟 | 5-10 分钟 | - |
| 解缠（SNAPHU） | 10-30 分钟 | 10-30 分钟 | 取决于相干性 |
| 地理编码 | 5-10 分钟 | 5-10 分钟 | - |
| **总计** | **50-85 分钟** | **26-65 分钟** | **配准加速显著** |

---

## 10. 技术支持

- **项目**: PowerChina-1 (Spacety) DInSAR 系统
- **版本**: v3.1
- **文档更新**: 2026-01-19
- **适用数据**: 天仪BC3 / 涪城一号 / DJ1 卫星 C波段 StripMap 模式

---

## 附录A: 依赖版本汇总表

| 组件 | 版本 | 必需 |
|------|------|------|
| Ubuntu | 22.04 LTS | ✅ |
| GCC | 11.4.0+ | ✅ |
| GMT | 6.3.0+ | ✅ |
| Python3 | 3.10+ | ✅ |
| NumPy | 1.21+ | ✅ |
| NetCDF | 4.8.1 | ✅ |
| HDF5 | 1.10.7 | ✅ |
| FFTW3 | 3.3.8 | ✅ |
| LAPACK | 3.10.0 | ✅ |
| BLAS | 3.10.0 | ✅ |
| TCSH | 6.21.00 | ✅ |
| Perl | 5.34.0 | ✅ |
| Python | 3.10.12 | ⚪ (可选) |
| SNAPHU | 2.0.7 | ✅ |
| Autoconf | 2.71 | ✅ |

## 附录B: 快速验证命令

```bash
#!/bin/bash
# 环境验证脚本

echo "=== 涪城一号 DInSAR 环境验证 ==="

# 检查核心工具
check_tool() {
    if command -v $1 &> /dev/null; then
        echo "✅ $1 已安装"
    else
        echo "❌ $1 未找到"
    fi
}

check_tool gmt
check_tool make_slc_s1a
check_tool snaphu
check_tool csh
check_tool perl

# 检查库
echo ""
echo "=== 依赖库检查 ==="
ldconfig -p | grep -q libgmt && echo "✅ GMT库" || echo "❌ GMT库"
ldconfig -p | grep -q libnetcdf && echo "✅ NetCDF库" || echo "❌ NetCDF库"
ldconfig -p | grep -q libfftw && echo "✅ FFTW库" || echo "❌ FFTW库"

echo ""
echo "=== 版本信息 ==="
gmt --version
```

---

*文档结束*
