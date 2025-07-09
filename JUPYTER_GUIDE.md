# MDemon 项目 Jupyter Notebook 使用指南

这份指南将帮助你在 MDemon 项目中使用 Jupyter Notebook 进行交互式开发和数据分析。

## 🚀 快速开始

### 1. 确保环境已准备

```bash
# 同步所有依赖（包括 Jupyter）
uv sync --extra dev

# 验证 Jupyter 是否可用
uv run jupyter --version
```

### 2. 启动 Jupyter

你有几种选择：

**选项 A: Jupyter Lab（推荐）**
```bash
uv run jupyter lab
```

**选项 B: 经典 Jupyter Notebook**
```bash
uv run jupyter notebook
```

**选项 C: 指定端口启动**
```bash
uv run jupyter lab --port=8888
```

### 3. 选择正确的 Kernel

在 Jupyter 中创建新 notebook 时，请选择：
- **"MDemon (uv)"** - 这是我们专门为项目创建的 kernel
- 或者 **"Python 3 (ipykernel)"** - 虚拟环境中的默认 kernel

## 📊 功能验证

### 测试 MDemon 导入

在 notebook cell 中运行：

```python
# 基本导入测试
import MDemon as md
print("✅ MDemon 导入成功")

# Cython 扩展测试
import MDemon.lib.c_distances
import MDemon.core.source
print("✅ Cython 扩展正常工作")

# 版本信息
import sys
print(f"Python 版本: {sys.version}")
print(f"工作目录: {import os; os.getcwd()}")
```

### 创建和分析数据

```python
# 使用测试数据
import os
test_file = os.path.join("tests", "data", "lammps", "CNT", 
                        "SWNT-8-8-graphene_hole_h2o_random-66666-1-10-30-3.data")

if os.path.exists(test_file):
    u = md.Universe(test_file)
    print(f"原子数: {len(u.atoms)}")
    print(f"分子数: {len(u.molecules)}")
```

## 🛠️ 开发工作流

### 修改 Cython 代码后的重编译

如果你修改了 `.pyx` 文件，需要重新编译：

```bash
# 强制重新安装项目（会重新编译 Cython）
uv pip install -e . --force-reinstall --no-deps

# 或者清理后重新安装
rm -f MDemon/lib/*.so MDemon/core/*.so
uv pip install -e .
```

然后在 Jupyter 中重启 kernel：
- **Kernel** → **Restart Kernel**

### 热重载模块

在开发过程中，如果修改了 Python 代码（非 Cython），可以使用：

```python
# 在 notebook 开头添加
%load_ext autoreload
%autoreload 2

# 现在修改的 Python 模块会自动重新加载
```

## 📈 数据可视化示例

### 基本绘图设置

```python
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# 设置绘图样式
plt.style.use('default')
sns.set_palette("husl")
%matplotlib inline  # 在 notebook 中内联显示图片
```

### 原子坐标可视化

```python
# 提取坐标数据
if 'u' in locals():
    coords = np.array([atom.coordinate for atom in u.atoms[:100]])
    
    # 创建 3D 投影图
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # XY 投影
    axes[0].scatter(coords[:, 0], coords[:, 1], alpha=0.6)
    axes[0].set_xlabel('X (Å)')
    axes[0].set_ylabel('Y (Å)')
    axes[0].set_title('XY 投影')
    
    # XZ 投影
    axes[1].scatter(coords[:, 0], coords[:, 2], alpha=0.6, color='orange')
    axes[1].set_xlabel('X (Å)')
    axes[1].set_ylabel('Z (Å)')
    axes[1].set_title('XZ 投影')
    
    # YZ 投影
    axes[2].scatter(coords[:, 1], coords[:, 2], alpha=0.6, color='green')
    axes[2].set_xlabel('Y (Å)')
    axes[2].set_ylabel('Z (Å)')
    axes[2].set_title('YZ 投影')
    
    plt.tight_layout()
    plt.show()
```

## 🔬 高级功能

### 距离计算性能分析

```python
from MDemon.lib.distance import distance_array
import time

# 创建测试数据
coords1 = np.random.randn(100, 3).astype(np.float32)
coords2 = np.random.randn(100, 3).astype(np.float32)

# 性能测试
start_time = time.time()
distances = distance_array(coords1, coords2)
calc_time = time.time() - start_time

print(f"计算 {len(coords1)}x{len(coords2)} 距离矩阵用时: {calc_time:.4f}秒")
print(f"结果形状: {distances.shape}")
```

### 分子网络分析

```python
import networkx as nx

if 'u' in locals() and len(u.bonds) > 0:
    # 创建分子网络图
    G = nx.Graph()
    
    # 添加原子节点
    for i, atom in enumerate(u.atoms[:50]):  # 只用前50个原子
        G.add_node(i, element=atom.element)
    
    # 添加键边
    for bond in u.bonds[:100]:  # 只用前100个键
        atoms = list(bond.atm_top)
        if len(atoms) >= 2 and all(a < 50 for a in atoms):
            G.add_edge(atoms[0], atoms[1])
    
    # 可视化网络
    plt.figure(figsize=(12, 8))
    pos = nx.spring_layout(G, k=1, iterations=50)
    nx.draw(G, pos, node_size=20, alpha=0.8, edge_color='gray', width=0.5)
    plt.title("分子键网络图")
    plt.show()
    
    print(f"网络节点数: {G.number_of_nodes()}")
    print(f"网络边数: {G.number_of_edges()}")
```

## 💡 最佳实践

### 1. 项目结构

```
your_analysis/
├── notebooks/
│   ├── 01_data_exploration.ipynb
│   ├── 02_distance_analysis.ipynb
│   └── 03_visualization.ipynb
├── data/
│   └── your_trajectory_files.data
└── results/
    └── figures/
```

### 2. Notebook 组织

```python
# 在每个 notebook 开头添加
# 标准导入
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# MDemon 导入
import MDemon as md
from MDemon.lib.distance import distance_array, self_distance_array

# 设置
%matplotlib inline
%load_ext autoreload
%autoreload 2
sns.set_style("whitegrid")
```

### 3. 数据管理

```python
# 定义数据路径
DATA_DIR = "data"
RESULTS_DIR = "results"

# 创建结果目录
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(f"{RESULTS_DIR}/figures", exist_ok=True)

# 保存图表
def save_figure(name, dpi=300):
    plt.savefig(f"{RESULTS_DIR}/figures/{name}.png", dpi=dpi, bbox_inches='tight')
    plt.savefig(f"{RESULTS_DIR}/figures/{name}.pdf", bbox_inches='tight')
```

## 🆘 故障排除

### 问题 1: Kernel 连接失败

```bash
# 重新安装 kernel
uv run python -m ipykernel install --user --name=mdemon --display-name="MDemon (uv)" --force

# 重启 Jupyter
# Ctrl+C 停止，然后重新运行 uv run jupyter lab
```

### 问题 2: 模块导入失败

```bash
# 确认项目已正确安装
uv pip list | grep mdemon

# 重新安装项目
uv pip install -e . --force-reinstall
```

### 问题 3: Cython 扩展问题

```bash
# 检查编译产物
ls -la MDemon/lib/*.so MDemon/core/*.so

# 清理重编译
find . -name "*.so" -delete
find . -name "*.c" -delete
uv pip install -e . --force-reinstall
```

### 问题 4: 内存不足

```python
# 在 notebook 中监控内存使用
import psutil
import os

def memory_usage():
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024  # MB

print(f"当前内存使用: {memory_usage():.1f} MB")
```

## 📚 更多资源

- [Jupyter Lab 用户指南](https://jupyterlab.readthedocs.io/)
- [MDAnalysis 教程](https://userguide.mdanalysis.org/) (类似项目参考)
- [NumPy 文档](https://numpy.org/doc/)
- [Matplotlib 画廊](https://matplotlib.org/stable/gallery/)

## 🔄 更新和维护

定期更新依赖：

```bash
# 更新所有依赖到最新兼容版本
uv sync --upgrade

# 检查过时的包
uv pip list --outdated
```

记住在重大更新后测试 Cython 扩展是否仍然正常工作！

---

**愉快的科学计算！** 🧬⚛️📊 