# MDemon 配位数分析指南

## 概述

配位数分析模块是MDemon单原子分析系统的重要组成部分，专门用于分析每个原子的配位环境。与RDF和角分布分析不同，配位数分析重点关注**个体原子的配位数分布**，而不是平均行为。

## 核心特点

### 1. 个体数据保留
- **保留每个原子的具体配位数值**
- 提供完整的个体配位数分布统计
- 支持按原子类型分组的配位数分析

### 2. 与RDF/Angular分析的区别

| 特性 | 配位数分析 | RDF/Angular分析 |
|------|------------|-----------------|
| **数据保存** | 每个原子的具体配位数 | 平均分布函数 |
| **分析重点** | 个体差异和分布统计 | 整体平均行为 |
| **可视化** | 直方图、散点图、热图 | 曲线图 |
| **应用场景** | 局部结构、缺陷分析 | 整体结构特征 |

### 3. 支持的功能
- 多种原子类型的配位数计算
- 个体配位数分布统计
- 多样化的可视化方式
- 并行处理和空间预分割优化

## 使用方法

### 1. 基本导入

```python
from MDemon.analysis.single_atom import CoordinationAnalyzer, CoordinationResult
```

### 2. 创建分析器

```python
# 定义截断半径（必需参数）
cutoff_radii = {
    '1-1': 2.5,  # 类型1-类型1的截断半径
    '1-2': 3.0,  # 类型1-类型2的截断半径
    '2-2': 2.8,  # 类型2-类型2的截断半径
}

# 创建配位数分析器
analyzer = CoordinationAnalyzer(
    universe,
    cutoff_radii=cutoff_radii,
    atom_selection=atom_mask,  # 可选：原子选择掩码
    scheduler='threads',       # 并行调度器
    enable_spatial_subdivision=True  # 启用空间预分割优化
)
```

### 3. 进行分析

```python
# 分析所有选中的原子
result = analyzer.analyze_parallel()

# 分析特定原子
result = analyzer.analyze_parallel(atom_indices=[0, 1, 2, 3, 4])
```

### 4. 获取结果

```python
# 获取特定原子的配位数
coord_numbers = result.get_coordination_number(atom_index=0)
print(f"原子0的配位数: {coord_numbers}")

# 获取特定原子对特定物种的配位数
coord_to_species1 = result.get_coordination_number(atom_index=0, species='1')
print(f"原子0对物种1的配位数: {coord_to_species1}")

# 获取配位数分布统计
distribution = result.get_coordination_distribution(species='1')
print(f"物种1的配位数分布: {distribution}")
```

### 5. 可视化

```python
import matplotlib.pyplot as plt

# 创建多种可视化图表
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 配位数分布直方图
result.plot(plot_type='histogram', ax=axes[0,0])

# 平均配位数条形图
result.plot(plot_type='bar', ax=axes[0,1])

# 原子索引vs配位数散点图
result.plot(plot_type='scatter', ax=axes[1,0])

# 配位数热图
result.plot(plot_type='heatmap', ax=axes[1,1])

plt.tight_layout()
plt.show()
```

### 6. 统计信息

```python
# 打印详细的统计信息
result.print_statistics()

# 获取配位数矩阵
atom_indices, species_list, coord_matrix = result.get_coordination_matrix()
print(f"配位数矩阵形状: {coord_matrix.shape}")
```

## API 参考

### CoordinationAnalyzer 类

#### 构造函数参数
- `universe`: MDemon.Universe对象
- `cutoff_radii`: 截断半径字典（必需）
- `atom_selection`: 原子选择掩码（可选）
- `scheduler`: Dask调度器类型
- `enable_spatial_subdivision`: 是否启用空间预分割

#### 主要方法
- `analyze_single_atom(atom_index)`: 分析单个原子
- `analyze_parallel(atom_indices=None)`: 并行分析多个原子
- `compute_coordination_statistics()`: 计算统计信息

### CoordinationResult 类

#### 主要方法
- `get_coordination_number(atom_index, species=None)`: 获取配位数
- `get_coordination_distribution(species=None)`: 获取分布统计
- `get_coordination_matrix()`: 获取配位数矩阵
- `plot(plot_type='histogram')`: 绘制图表
- `print_statistics()`: 打印统计信息

#### 支持的绘图类型
- `'histogram'`: 配位数分布直方图
- `'bar'`: 平均配位数条形图
- `'scatter'`: 原子索引vs配位数散点图
- `'heatmap'`: 配位数热图

## 应用示例

### 1. 缺陷分析
```python
# 分析特定区域的配位数分布，识别缺陷
result = analyzer.analyze_parallel(atom_indices=defect_region_atoms)
distribution = result.get_coordination_distribution()

# 找出配位数异常的原子
for atom_idx in result.atom_indices:
    coord_numbers = result.get_coordination_number(atom_idx)
    if any(cn < expected_min or cn > expected_max for cn in coord_numbers.values()):
        print(f"原子 {atom_idx} 可能存在缺陷: {coord_numbers}")
```

### 2. 结构多样性分析
```python
# 分析整个系统的配位数分布多样性
result = analyzer.analyze_parallel()
distribution = result.get_coordination_distribution()

for species, data in distribution.items():
    stats = data['statistics']
    print(f"物种 {species}:")
    print(f"  配位数标准差: {stats['std']:.2f}")
    print(f"  配位数范围: {stats['min']} - {stats['max']}")
```

### 3. 界面分析
```python
# 分析界面区域的配位环境
interface_mask = create_interface_mask(universe)
result = analyzer.analyze_parallel()

# 比较界面和体相的配位数分布
interface_atoms = np.where(interface_mask)[0]
bulk_atoms = np.where(~interface_mask)[0]

# 分别统计界面和体相的配位数
interface_coords = [result.get_coordination_number(idx) for idx in interface_atoms]
bulk_coords = [result.get_coordination_number(idx) for idx in bulk_atoms]
```

## 性能优化

### 1. 空间预分割
配位数分析自动启用空间预分割优化，显著提高大系统的计算效率：

```python
analyzer = CoordinationAnalyzer(
    universe,
    cutoff_radii=cutoff_radii,
    enable_spatial_subdivision=True,  # 默认启用
    use_subdivision_masks=False       # 大系统建议False
)
```

### 2. 并行处理
支持多种并行策略：

```python
# 线程并行（适合I/O密集型）
analyzer = CoordinationAnalyzer(universe, cutoff_radii, scheduler='threads')

# 进程并行（适合CPU密集型）
analyzer = CoordinationAnalyzer(universe, cutoff_radii, scheduler='processes')

# 分布式并行（适合超大系统）
analyzer = CoordinationAnalyzer(universe, cutoff_radii, scheduler='distributed')
```

## 注意事项

1. **截断半径设置**: 必须为每个原子类型对设置合适的截断半径
2. **内存使用**: 配位数分析保留个体数据，内存使用比RDF/Angular分析稍高
3. **结果解释**: 配位数结果是整数值，反映的是邻近原子的数量
4. **可视化**: 中文标签可能在某些系统上显示不正确，建议使用英文标签

## 与其他分析的结合

配位数分析可以与RDF和角分布分析结合使用：

```python
from MDemon.analysis.single_atom import SingleAtomAnalysis

# 创建统一的分析接口
analysis = SingleAtomAnalysis(universe)

# 进行多种分析
rdf_result = analysis.rdf(atom_selection=mask, r_range=(0, 10))
angular_result = analysis.angular(atom_selection=mask, cutoff_radii=cutoffs)
coord_result = analysis.coordination(atom_selection=mask, cutoff_radii=cutoffs)

# 结合分析结果
# 例如：高配位数原子的RDF特征
high_coord_atoms = [idx for idx in coord_result.atom_indices 
                   if sum(coord_result.get_coordination_number(idx).values()) > threshold]
```

## 总结

配位数分析模块为MDemon提供了强大的局部结构分析能力，特别适用于：
- 缺陷和界面分析
- 结构多样性研究
- 局部环境表征
- 与RDF/Angular分析的互补研究

通过保留个体原子的配位数信息，该模块填补了MDemon在局部结构统计分析方面的空白，与现有的RDF和角分布分析形成了完整的单原子分析体系。
