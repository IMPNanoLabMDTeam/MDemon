# XYZ Reader 使用指南

MDemon 支持读取 **Extended XYZ** 格式文件，这是一种广泛用于分子动力学模拟的标准格式，被 GPUMD、ASE 等多个软件包采用。

## 目录

1. [Extended XYZ 格式规范](#extended-xyz-格式规范)
2. [基本使用](#基本使用)
3. [支持的属性](#支持的属性)
4. [高级特性](#高级特性)
5. [示例](#示例)
6. [与其他格式的比较](#与其他格式的比较)

---

## Extended XYZ 格式规范

Extended XYZ 格式是标准 XYZ 格式的扩展版本，格式说明详见 [GPUMD 文档](https://gpumd.org/gpumd/input_files/model_xyz.html)。

### 文件结构

```
N                                           # 第1行：原子数量
Lattice="..." Properties=...               # 第2行：元数据（晶格、属性定义等）
species x y z [additional columns...]      # 第3行起：原子数据
species x y z [additional columns...]
...
```

### 第1行：原子数量

一个整数，表示文件中包含的原子总数。

```
80
```

### 第2行：元数据

包含多个 `keyword=value` 键值对，用空格分隔。值可以用双引号包围。

#### 必需关键词

1. **`Lattice`**（必需）：晶格向量
   ```
   Lattice="ax ay az bx by bz cx cy cz"
   ```
   定义三个晶格向量：
   - **a** = (ax, ay, az)
   - **b** = (bx, by, bz)
   - **c** = (cx, cy, cz)
   
   单位：Ångström (Å)

2. **`Properties`**（必需）：数据列定义
   ```
   Properties=prop1:type1:ncols1:prop2:type2:ncols2:...
   ```
   - `prop`: 属性名称（如 `species`, `pos`, `vel`）
   - `type`: 数据类型
     - `S`: 字符串（String）
     - `R`: 实数（Real）
     - `I`: 整数（Integer）
   - `ncols`: 列数

#### 可选关键词

- **`pbc`**：周期性边界条件
  ```
  pbc="T F F"  # x 方向周期，y 和 z 方向自由
  ```
  默认值：`"T T T"`（三个方向都周期）

### 第3行起：原子数据

每行包含一个原子的数据，列数和含义由 `Properties` 定义。

---

## 基本使用

### 读取 XYZ 文件

```python
import MDemon as md

# 读取 Extended XYZ 文件
u = md.Universe("system.xyz")

# 访问原子信息
print(f"Total atoms: {len(u.atoms)}")
print(f"First atom coordinate: {u.atoms[0].coordinate}")
print(f"First atom species: {u.atoms[0].species}")

# 访问盒子信息
print(f"Box dimensions: {u.box}")
```

### 文件格式自动识别

MDemon 会根据文件扩展名自动识别格式：
- `.xyz` → XYZ Reader
- `.exyz` → XYZ Reader
- `.data` → LAMMPS Reader
- `.reaxff` → REAXFF Reader

也可以手动指定格式：
```python
u = md.Universe("system.dat", format="XYZ")
```

---

## 支持的属性

MDemon XYZ Reader 支持以下属性：

### 必需属性

| 属性 | 格式 | 说明 | MDemon 对应属性 |
|------|------|------|-----------------|
| `species` | `S:1` | 原子类型/元素符号 | `atom.species` |
| `pos` | `R:3` | 位置坐标 (x, y, z) | `atom.coordinate` |

### 可选属性

| 属性 | 格式 | 说明 | MDemon 对应属性 |
|------|------|------|-----------------|
| `mass` | `R:1` | 质量（amu） | `atom.mass` |
| `charge` | `R:1` | 电荷（e） | `atom.charge` |
| `vel` | `R:3` | 速度 (vx, vy, vz) | `atom.velocity` |
| `group` | `I:N` | 分组标签（N个分组方法） | 暂不支持 |

**注意**：
- 如果文件中不包含 `mass`，MDemon 会根据元素符号自动查询周期表获取默认质量
- 如果文件中不包含 `charge`，默认为 0
- 如果文件中不包含 `vel`，默认为零向量

---

## 高级特性

### 1. 晶格类型自动识别

MDemon 会自动识别晶格类型（正交或三斜）：

```python
# 正交晶格
# Lattice="10 0 0 0 10 0 0 0 10"
# 生成：box = [10, 10, 10, 90, 90, 90, 0, 10, 0, 10, 0, 10]

# 三斜晶格
# Lattice="10 0 0 2 10 0 1 1 10"
# 生成：box = 3x3 矩阵
```

### 2. 元素符号映射

Reader 会将元素符号映射为整数类型：
- 文件中的不同元素符号会自动编号（从 1 开始）
- 相同元素符号共享相同的类型编号

```python
# 文件中：Ga Ga O O Ga O
# 映射为：1  1  2 2 1  2
```

### 3. 质量自动查询

如果文件不包含质量信息，MDemon 会自动从 `periodictable` 库查询：

```python
# Ga → 69.723 amu
# O  → 15.999 amu
```

---

## 示例

### 示例 1：基本 Extended XYZ 文件

```xyz
4
Lattice="10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 10.0" Properties=species:S:1:pos:R:3
C 0.0 0.0 0.0
C 2.5 2.5 0.0
O 5.0 0.0 0.0
O 7.5 2.5 0.0
```

读取和分析：

```python
import MDemon as md

u = md.Universe("simple.xyz")

print(f"Atoms: {len(u.atoms)}")  # 4
print(f"Species: {set([a.species for a in u.atoms])}")  # {1, 2}

# 碳原子（species=1）
carbon_atoms = [i for i, a in enumerate(u.atoms) if a.species == 1]
print(f"Carbon atoms: {carbon_atoms}")  # [0, 1]

# 氧原子（species=2）
oxygen_atoms = [i for i, a in enumerate(u.atoms) if a.species == 2]
print(f"Oxygen atoms: {oxygen_atoms}")  # [2, 3]
```

### 示例 2：包含速度的文件

```xyz
2
Lattice="5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0" Properties=species:S:1:pos:R:3:vel:R:3
C 0.0 0.0 0.0 0.001 0.002 0.003
C 2.5 0.0 0.0 -0.001 0.001 0.000
```

读取和访问速度：

```python
import MDemon as md

u = md.Universe("with_velocity.xyz")

for i, atom in enumerate(u.atoms):
    print(f"Atom {i}:")
    print(f"  Position: {atom.coordinate}")
    print(f"  Velocity: {atom.velocity}")
```

### 示例 3：包含质量和电荷的文件

```xyz
2
Lattice="5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0" Properties=species:S:1:pos:R:3:mass:R:1:charge:R:1
C 0.0 0.0 0.0 12.011 0.0
O 2.5 0.0 0.0 15.999 -0.5
```

读取和访问质量、电荷：

```python
import MDemon as md

u = md.Universe("with_mass_charge.xyz")

for atom in u.atoms:
    print(f"Species: {atom.species}")
    print(f"Mass: {atom.mass}")
    print(f"Charge: {atom.charge}")
```

### 示例 4：三斜晶格（非正交盒子）

```xyz
2
Lattice="10.0 0.0 0.0 2.0 10.0 0.0 1.0 1.0 10.0" Properties=species:S:1:pos:R:3
C 0.0 0.0 0.0
C 5.0 5.0 5.0
```

三斜晶格会被存储为 3×3 矩阵：

```python
import MDemon as md

u = md.Universe("triclinic.xyz")

print(f"Box shape: {u.box.shape}")  # (3, 3)
print(f"Lattice vectors:\n{u.box}")
# [[10.  0.  0.]
#  [ 2. 10.  0.]
#  [ 1.  1. 10.]]
```

### 示例 5：实际应用 - Ga₂O₃ 晶体

这是一个实际的 Ga₂O₃ (β-Ga₂O₃) 晶体结构文件：

```python
import MDemon as md
import numpy as np

# 读取晶体结构
u = md.Universe("beta_221_ort.xyz")

print(f"System: β-Ga₂O₃")
print(f"Total atoms: {len(u.atoms)}")

# 统计元素组成
species_list = [atom.species for atom in u.atoms]
species_counts = {}
for s in species_list:
    species_counts[s] = species_counts.get(s, 0) + 1

print(f"\nComposition:")
for species, count in species_counts.items():
    print(f"  Species {species}: {count} atoms")

# 晶格参数
box = u.box
print(f"\nLattice parameters:")
print(f"  a = {box[0]:.6f} Å")
print(f"  b = {box[1]:.6f} Å")
print(f"  c = {box[2]:.6f} Å")

# 计算密度
volume = box[0] * box[1] * box[2]  # Å³
total_mass = sum([atom.mass for atom in u.atoms])  # amu
density = total_mass / volume * 1.66054  # g/cm³

print(f"\nDensity: {density:.3f} g/cm³")
```

---

## 与其他格式的比较

### XYZ vs LAMMPS DATA

| 特性 | Extended XYZ | LAMMPS DATA |
|------|--------------|-------------|
| 文件结构 | 简洁，易读 | 复杂，多段落 |
| 晶格定义 | 直接定义晶格向量 | 定义盒子边界 |
| 原子样式 | 统一的 Properties 定义 | 多种原子样式 |
| 拓扑信息 | 不包含 | 包含键、角等 |
| 速度 | 可选，内联 | 可选，独立段落 |
| 最适用于 | 晶体结构、初始构型 | 完整分子系统 |

### XYZ vs REAXFF

| 特性 | Extended XYZ | REAXFF |
|------|--------------|--------|
| 主要用途 | 原子位置 | 反应键信息 |
| 键级 | 不包含 | 包含 |
| 孤对电子 | 不包含 | 包含 |
| 与 DATA 配合 | 独立使用 | 需要配合 DATA |

### 推荐使用场景

**使用 Extended XYZ 格式：**
- ✅ 晶体结构文件
- ✅ 与其他软件（GPUMD、ASE）交换数据
- ✅ 简单分子系统
- ✅ 原子位置和速度数据

**使用 LAMMPS DATA 格式：**
- ✅ 复杂分子系统
- ✅ 需要拓扑信息（键、角、二面角）
- ✅ 力场参数定义
- ✅ LAMMPS 模拟输入

---

## 写入 XYZ 文件

MDemon 也支持将数据写入 Extended XYZ 格式（功能开发中）：

```python
from MDemon.reader.XYZ import XYZWriter

# 创建写入器
writer = XYZWriter("output.xyz")

# 写入 Universe
writer.write(u, properties=['species', 'pos', 'vel'], 
             include_velocities=True)
```

---

## 格式参考

完整的 Extended XYZ 格式规范请参考：
- [GPUMD Documentation](https://gpumd.org/gpumd/input_files/model_xyz.html)
- [ASE Extended XYZ Format](https://wiki.fysik.dtu.dk/ase/ase/io/formatoptions.html#extxyz)

---

## 常见问题

### Q1: 为什么元素符号变成了数字？

MDemon 内部使用整数来表示原子类型（species），这样可以更高效地处理数据。元素符号会被自动映射为整数：

```python
# 文件中：C C O N C O
# 内部表示：1 1 2 3 1 2
```

如果需要元素信息，可以使用 `atom.element` 属性（从质量推断）。

### Q2: 如何处理没有质量信息的文件？

MDemon 会自动从周期表查询默认质量：

```python
# Ga → 69.723 amu (自动查询)
# O  → 15.999 amu (自动查询)
```

如果查询失败，会使用 `species * 12.0` 作为占位值。

### Q3: 三斜晶格如何处理？

三斜晶格会被存储为 3×3 矩阵：

```python
if len(u.box.shape) == 2:  # 三斜
    lattice_vectors = u.box
else:  # 正交
    lx, ly, lz = u.box[:3]
```

### Q4: 支持周期性边界条件吗？

是的！在元数据行添加 `pbc` 关键词：

```
Lattice="..." Properties=... pbc="T T F"
```

表示 x、y 方向周期，z 方向自由。

---

## 总结

Extended XYZ 格式的优势：
- ✅ **简洁明了**：易于人工阅读和编辑
- ✅ **标准化**：被多个软件包采用
- ✅ **灵活**：通过 Properties 定义支持各种数据
- ✅ **完整**：包含晶格、坐标、速度等完整信息

现在你可以使用 MDemon 轻松读取和处理 Extended XYZ 格式的文件了！

---

**文档版本**：v1.0  
**最后更新**：2025-10-25  
**作者**：MDemon Development Team

