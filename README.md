# MDemon

一个强大的分子动力学(MD)仿真数据分析工具，专注于高性能计算和辐照效应分析。

## 项目简介

MDemon 是一个现代化的分子动力学分析框架，提供全面的分析工具和物理常数库。主要特性包括：

- 🚀 **高性能计算**: 使用 Cython 扩展优化核心算法
- 📊 **全面分析**: 支持径向分布函数(RDF)、角分布分析等
- ⚛️ **辐照效应**: 集成 Waligorski-Zhang 计算器用于辐照损伤分析
- 🔬 **多格式支持**: 兼容 LAMMPS、REAXFF、Extended XYZ 等主流 MD 文件格式
- 📚 **丰富常数库**: 内置物理、化学、材料常数和单位转换
- 🔧 **灵活扩展**: 模块化设计，易于扩展和定制

## 快速开始

### 环境配置

推荐使用 [uv](https://docs.astral.sh/uv/) 管理 Python 环境：

```bash
# 安装 uv (如果尚未安装)
curl -LsSf https://astral.sh/uv/install.sh | sh

# 创建并激活虚拟环境
uv venv
source .venv/bin/activate  # Linux/macOS
# 或 .venv\Scripts\activate  # Windows

# 安装项目依赖
uv pip install -e .

# 安装开发依赖 (可选)
uv pip install -e ".[dev]"
```

### 基本使用

#### 读取分子动力学文件

```python
import MDemon as md

# 读取 LAMMPS DATA 文件
u = md.Universe("system.data")

# 读取 Extended XYZ 文件（如 GPUMD 格式）
u = md.Universe("system.xyz")

# 读取多个文件（DATA + REAXFF）
u = md.Universe("system.data", "bonds.reaxff")

# 访问原子信息
print(f"Total atoms: {len(u.atoms)}")
print(f"First atom coordinate: {u.atoms[0].coordinate}")
```

#### 使用物理常数

```python
from MDemon.constants import PHYSICS, CHEMISTRY
from MDemon.constants.utils import convert_energy

# 使用物理常数
print(f"玻尔兹曼常数: {md.k_B:.6e} J/K")
print(f"阿伏伽德罗数: {md.N_A:.6e} mol⁻¹")

# 单位转换
energy_ev = 1.0
energy_j = convert_energy(energy_ev, 'eV', 'J')
print(f"1 eV = {energy_j:.6e} J")
```

#### 辐照效应分析

```python
from MDemon.utils import WaligorskiZhangCalculator
import numpy as np

# 辐照计算示例
calc = WaligorskiZhangCalculator.from_material_preset('Ga2O3')
radius_array = np.logspace(0, 3, 100)  # 1 to 1000 nm
result = calc.calculate_radial_dose(radius_array, ion_Z=8, ion_energy=10.0)
```

## 项目结构

```
MDemon/
├── analysis/           # 分析模块
│   └── single_atom/   # 单原子分析工具
│       ├── angular.py # 角分布分析
│       ├── rdf.py     # 径向分布函数
│       └── ...
├── constants/         # 物理化学常数库
│   ├── chemistry.py   # 化学常数
│   ├── conversion.py  # 单位转换
│   ├── fundamental.py # 基本物理常数
│   ├── irradiation.py # 辐照相关常数
│   ├── material.py    # 材料属性
│   └── utils.py       # 工具函数
├── core/              # 核心数据结构
│   ├── database.py    # 数据库接口
│   ├── structure.py   # 分子结构
│   ├── universe.py    # 模拟体系
│   └── source.pyx     # Cython 优化模块
├── lib/               # 高性能计算库
│   ├── c_distances.pyx # 距离计算 (Cython)
│   ├── distance.py    # 距离分析
│   └── include/       # C/C++ 头文件
├── reader/            # 文件读取器
│   ├── LAMMPS.py      # LAMMPS 格式
│   ├── REAXFF.py      # REAXFF 格式
│   ├── XYZ.py         # Extended XYZ 格式
│   └── base.py        # 基础读取器
├── selection/         # 原子选择工具
│   └── cylindrical.py # 圆柱形选择
└── utils/             # 实用工具
    ├── irradiation.py # 辐照分析工具
    ├── parallel.py    # 并行计算
    └── spatial.py     # 空间分析
```

### 主要模块说明

- **analysis**: 提供各种分析算法，包括径向分布函数、角分布分析等
- **constants**: 完整的物理、化学常数库和单位转换系统
- **core**: 核心数据结构，管理分子结构和模拟体系
- **lib**: 使用 Cython 优化的高性能计算模块
- **reader**: 支持多种 MD 文件格式的读取器
- **utils**: 实用工具，包括辐照效应分析和并行计算支持

## 功能特性

### 物理常数与单位转换

```python
from MDemon.constants.utils import convert_temperature, boltzmann_energy

# 温度转换
temp_k = convert_temperature(25.0, 'C', 'K')  # 298.15 K

# 热能计算
thermal_energy = boltzmann_energy(298.15, 'eV')  # eV 单位
```

### 材料属性查询

```python
from MDemon.constants.utils import get_material_property

# 获取 β-Ga₂O₃ 的密度
density = get_material_property('Ga2O3_beta', 'density')
print(f"β-Ga₂O₃ 密度: {density} kg/m³")
```

### 辐照效应分析

```python
from MDemon.utils import WaligorskiZhangCalculator

# 创建辐照计算器
calc = WaligorskiZhangCalculator.from_material_preset('Ga2O3')

# 计算径向剂量分布
dose_profile = calc.calculate_radial_dose(
    radius_nm=np.logspace(0, 3, 100),
    ion_Z=8,           # 氧离子
    ion_energy=10.0    # MeV/amu
)
```

## 文档和示例

### 完整文档

- [文件读取架构指南](docs/file_reading_architecture.md) - MDemon 核心文件读取逻辑详解
- [XYZ Reader 使用指南](docs/xyz_reader_guide.md) - Extended XYZ 格式读取指南
- [常数使用指南](docs/constants_guide.md) - 物理化学常数库使用说明
- [辐照分析指南](IRRADIATION_USAGE.md) - 辐照效应分析教程
- [Jupyter 使用指南](JUPYTER_GUIDE.md) - 在 Jupyter 中使用 MDemon

### 示例代码

- [XYZ Reader 示例](examples/xyz_reader_example.py) - Extended XYZ 文件读取示例
- [RDF 分析示例](examples/rdf_analysis_irradiated_kapton.py) - 径向分布函数分析
- [配位数分析示例](examples/coordination_analysis_example.py) - 配位数分析演示
- [常数使用示例](examples/constants_usage.py) - 物理常数和单位转换

## 开发和贡献

### 开发环境设置

```bash
# 克隆仓库
git clone <repository-url>
cd MDemon

# 设置开发环境
uv venv
source .venv/bin/activate
uv pip install -e ".[dev]"

# 运行测试
uv run pytest tests/

# 代码格式化
uv run black MDemon/
uv run isort MDemon/

# 代码检查
uv run ruff check MDemon/
```

### 构建 Cython 扩展

项目包含 Cython 扩展以提升性能，安装时会自动编译：

```bash
# 手动重新编译 Cython 扩展
uv run setup.py build_ext --inplace
```

## 系统要求

- Python 3.9+
- NumPy >= 1.24.0
- SciPy >= 1.10.0
- Pandas >= 2.0.0
- Matplotlib >= 3.6.0

## 许可证

本项目采用 GPL-2.0 许可证。详见 [LICENSE](LICENSE) 文件。

## 致谢

感谢所有为 MDemon 项目做出贡献的开发者和科研人员。

---

更多详细信息请参考项目文档和示例代码。如有问题，欢迎提交 Issue 或 Pull Request。
