# 电建一号 DInSAR 处理软件部署与操作手册

## 1. 系统概述
本系统基于定制版 GMTSAR，专为 **电建一号（DJ1）** 与 **涪城一号/天仪（BC3）** X 波段 SAR 数据量身打造。
系统集成了针对高分辨率 X 波段数据的参数优化，并包含了对 GCC 11+ 编译器的兼容性修复，实现了从 SLC 数据生成到干涉形变提取的自动化全流程处理。

主要功能亮点：
*   **支持多模式处理**：提供快速响应（A模式-Fast）和高精度（B模式-Quality）两种策略。
*   **自动化部署**：通过 `install.sh` 脚本实现一键式环境搭建。
*   **一键自检**：通过 `self_check.sh` 一次完成“清理→重装→冒烟→回归”。
*   **兼容性增强**：修复了 ALOS `extern` 变量和新版 Ubuntu 的编译报错。

## 2. 运行环境要求

为了确保数据处理的稳定性和效率，建议在以下环境中运行：

*   **操作系统**：Ubuntu 20.04 LTS / 22.04 LTS （已在 22.04 上通过全流程验证）。
*   **shell 环境**：Bash (推荐) 或 C-Shell (csh)。
*   **推荐硬件**：
    *   CPU: 8 核及以上建议。
    *   内存: 32GB RAM（Snaphu 解缠在高质量模式下可能占用大量内存）。
    *   存储: 至少 200GB 可用空间（SAR 影像解压和中间文件体积较大）。
*   **核心依赖包**：
    *   编译器: `gcc`, `make`, `gfortran` (需支持 -fcommon)
    *   GMT: `gmt` (版本 6.x), `libgmt-dev`
    *   库文件: `libtiff-dev` (TIFF支持), `libhdf5-dev` (NetCDF支持)
    *   其他: `csh`, `autoconf`, `ghostscript`

## 3. 快速部署指南

本系统已将复杂的编译配置封装为交互式脚本，用户无需手动编辑 Makefile。

### 3.1 核心部署脚本：`install.sh`

在软件安装包根目录下执行：

```bash
# 赋予执行权限
chmod +x install.sh

# 标准安装 (默认安装到 /usr/local/GMTSAR，需要 sudo 权限)
sudo ./install.sh
```

### 3.2 部署参数详解

脚本支持以下参数以适应不同场景：

1.  **非 Root 用户部署 (`--user`)**
    如果当前用户没有 sudo 权限，可安装到用户主目录 (`~/.local/gmtsar`)：
    ```bash
    ./install.sh --user
    ```

2.  **离线模式部署 (`--offline`)**
    在无外网环境下（如内网保密机），配合离线依赖包使用，跳过 git 拉取步骤：
    ```bash
    ./install.sh --offline
    ```

3.  **环境检查**
    脚本会自动检测并尝试修补 GCC 11+ 带来的代码兼容性问题（如 `typedef GMT_LONG` 冲突和 `extern` 全局变量问题）。

### 3.3 环境变量设置
`install.sh` 安装完成后会自动将 GMTSAR 写入 `~/.bashrc`。通常仅需重新加载：
```bash
source ~/.bashrc
```

若你需要手动修复 PATH，可执行（按安装前缀二选一）：
```bash
echo 'export PATH=/usr/local/GMTSAR/bin:$PATH' >> ~/.bashrc
# 或
echo 'export PATH=$HOME/.local/gmtsar/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

### 3.4 一键自检（推荐）

安装包内置一键自检脚本：`self_check.sh`

功能覆盖：
1. 清理旧 GMTSAR 软链接与旧安装目录；
2. 重新执行 `install.sh`；
3. 冒烟验证（`make_slc_dj1` / `xcorr` / `snaphu` / `gmt`）；
4. 自动挑选可用测试数据并执行 `run_dinsar.sh --mode A --auto`；
5. 校验结果文件（`los_displacement.grd` / `displacement.png` / `displacement.jpg`）。

使用方式：
```bash
cd /path/to/GMTSAR_Installer
chmod +x self_check.sh

# 自动选择内置测试数据
./self_check.sh

# 或指定你的 raw 数据目录
./self_check.sh /path/to/project/raw
```

## 4. 自动化处理流程

核心处理脚本为 **`run_dinsar.sh`**，它串联了该流程的所有步骤。

### 4.1 核心处理阶段

1.  **SLC 生成**：自动识别 `DJ1` 或 `BC3` 卫星元数据，调用 `make_slc_dj1` 将原始 TIFF 数据转换为 SLC（单视复数）格式。
2.  **高精度配准 (Registration)**：
    *   使用 `xcorr` 进行粗配准。
    *   执行精细配准（`xsearch` 搜索窗口自动调整为 128/256 像素），生成偏移量多项式。
3.  **相位处理**：
    *   `phasediff`：生成复数干涉图。
    *   `phasefilt`：应用 Goldstein 滤波去除相干噪声（支持 8x 强降采样以加速）。
4.  **相位解缠 (Unwrapping)**：
    *   调用 `snaphu` (Statistical-Cost, Network-Flow Algorithm for Phase Unwrapping)。
    *   支持**分块解缠（Tile Mode）**以利用多核 CPU 并节省内存。
5.  **地理编码 (Geocoding)**：
    *   将结果从雷达坐标系（Range-Azimuth）转换到地理坐标系（Lon-Lat）。
    *   生成最终交付的 PNG 和 GeoTIFF 文件。

### 4.2 处理模式选择

脚本提供了两种预设模式，通过 `--mode` 参数调用：

*   **Mode A (Fast / 快速模式)**：
    *   **适用场景**：应急响应、大范围普查、快速查看是否存在形变。
    *   **特点**：高降采样（默认 16x）、配准窗口收敛为 2 的幂（`xsearch/ysearch=64`），启用 Snaphu 分块（Tile 2x2，多线程）以缩短时间。
    *   **命令**：`./run_dinsar.sh --mode A ...`

*   **Mode B (Quality / 质量模式)**：
    *   **适用场景**：最终成果交付、微小形变分析。
    *   **特点**：低降采样 (2x/1x)，高精度配准窗口 (256)，严格的 Coherence 掩膜阈值。
    *   **命令**：`./run_dinsar.sh --mode B ...`

## 5. 部署验证

部署完成后，请执行以下命令检查三个核心二进制程序是否可用：

```bash
# 1. 检查主程序是否在路径中
which make_slc_dj1

# 2. 检查 SLC 生成模块
make_slc_dj1 2>&1 | head -n 3

# 3. 检查配准模块
xcorr 2>&1 | head -n 1

# 4. 检查解缠模块
snaphu -h | head -n 1

# 5. 推荐：一键自检（含回归）
./self_check.sh
```

当回归通过时，`/root/verify_dinsar_run/results`（或你指定目录）应至少出现以下之一：
* `los_displacement.grd`
* `displacement.png`（或 `displacement.jpg`）

## 6. 常见问题排除

### 6.1 Snaphu 报错：`TILECOSTTHRESH`
*   **现象**：日志提示 `Tile cost threshold too high` 且 `snaphu` 进程异常退出。
*   **修复**：`run_dinsar.sh` 已内置修复逻辑，启用 Tile Overlap 或降采样处理。

### 6.2 `xyz2grd` 格式报错
*   **现象**：提示 `grid dimension mismatch`。
*   **原因**：二进制文件大小与网格定义不符。
*   **修复**：脚本已强制指定 `-ZTLf` 参数（浮点型、从上到下、从小端序）。

### 6.3 `phasediff` 报错：`baseline < -9000 not set ?`
*   **现象**：在生成干涉图阶段失败。
*   **原因**：PRM 未注入基线参数。
*   **修复**：脚本在 `phasediff` 前自动调用 `SAT_baseline` 补齐 PRM 基线字段。

### 6.4 `snaphu` 报错：`extra data in file phase.bin (bad linelength?)`
*   **现象**：snaphu 读取 `phase.bin` 失败。
*   **原因**：未使用 `snaphu.conf.brief` 时，输入数据格式解释不一致。
*   **修复**：脚本自动探测并添加 `-f snaphu.conf.brief`。

### 6.5 ALOS 编译时的 `multiple definition` 错误
*   **现象**：编译阶段报错变量重复定义。
*   **修复**：`install.sh` 已包含代码自动修复补丁（Fix GCC 11+）。
