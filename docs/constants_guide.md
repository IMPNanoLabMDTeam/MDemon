# MDemon Physical Constants Guide

## 概述

MDemon的常量模块提供了分子动力学模拟中常用的物理常量、化学常量、单位转换因子和材料属性。这个模块的设计目标是：

- **准确性**: 基于2018年CODATA国际推荐值
- **便捷性**: 提供多种访问方式和便捷函数
- **完整性**: 涵盖MD模拟的常用常量和材料属性
- **扩展性**: 易于添加新的常量和材料

## 模块结构

```
MDemon/constants/
├── __init__.py          # 主入口，便捷导入
├── fundamental.py       # 基本物理常量
├── chemistry.py         # 化学相关常量
├── conversion.py        # 单位转换因子
├── material.py          # 材料属性
└── utils.py            # 实用工具函数
```

## 基本使用方法

### 1. 导入常量

```python
# 方法1: 导入常量类别
from MDemon.constants import PHYSICS, CHEMISTRY, CONVERSION, MATERIAL

# 方法2: 直接导入最常用的常量
from MDemon.constants import k_B, N_A, c_light, ANGSTROM_TO_METER, EV_TO_JOULE

# 方法3: 导入工具函数
from MDemon.constants import utils
```

### 2. 基本物理常量

```python
from MDemon.constants import PHYSICS

# 基本常量
print(f"Boltzmann constant: {PHYSICS.k_B} J/K")
print(f"Avogadro's number: {PHYSICS.N_A} mol⁻¹")
print(f"Speed of light: {PHYSICS.c} m/s")
print(f"Elementary charge: {PHYSICS.e} C")
print(f"Planck constant: {PHYSICS.h} J⋅s")

# 粒子质量
print(f"Electron mass: {PHYSICS.m_e} kg")
print(f"Proton mass: {PHYSICS.m_p} kg")
print(f"Atomic mass unit: {PHYSICS.u} kg")
```

### 3. 化学常量

```python
from MDemon.constants import CHEMISTRY

# 气体常数
print(f"Gas constant: {CHEMISTRY.R} J/(mol⋅K)")
print(f"Gas constant: {CHEMISTRY.R_cal} cal/(mol⋅K)")

# 标准条件
print(f"STP temperature: {CHEMISTRY.STP_temperature} K")
print(f"STP pressure: {CHEMISTRY.STP_pressure} Pa")

# 分子属性
bond_length = CHEMISTRY.get_bond_length('C-C')  # 返回以米为单位的键长
vdw_radius = CHEMISTRY.get_vdw_radius('H')      # 返回以米为单位的范德华半径
electronegativity = CHEMISTRY.get_electronegativity('C')  # 返回电负性
```

### 4. 单位转换

```python
from MDemon.constants import utils

# 能量转换
energy_j = utils.convert_energy(1.0, 'eV', 'J')           # 1 eV → 焦耳
energy_kcal = utils.convert_energy(1.0, 'Hartree', 'kcal/mol')  # 1 Hartree → kcal/mol

# 长度转换
length_m = utils.convert_length(5.0, 'Å', 'm')           # 5 Å → 米
length_nm = utils.convert_length(5.0, 'Å', 'nm')        # 5 Å → 纳米

# 温度转换
temp_k = utils.convert_temperature(25.0, 'C', 'K')      # 25°C → 开尔文
temp_f = utils.convert_temperature(298.15, 'K', 'F')    # 开尔文 → 华氏度
```

### 5. 材料属性

```python
from MDemon.constants import MATERIAL

# 获取材料属性
density = MATERIAL.get_density('Ga2O3_beta')             # β-Ga₂O₃密度
melting_point = MATERIAL.get_melting_point('Si')         # 硅的熔点
thermal_cond = MATERIAL.get_thermal_conductivity('Cu')   # 铜的热导率

# 晶格参数
params = MATERIAL.get_lattice_parameters('Ga2O3_beta')
print(f"β-Ga₂O₃ lattice: a={params['a']} Å, b={params['b']} Å")

# 列出所有可用材料
materials = MATERIAL.list_available_materials()
print(f"Available materials: {materials}")
```

## 高级功能

### 1. 热力学计算

```python
from MDemon.constants import utils, CHEMISTRY

# 计算热能 k_B * T
temp = 298.15  # K
thermal_energy_j = utils.boltzmann_energy(temp, 'J')
thermal_energy_ev = utils.boltzmann_energy(temp, 'eV')

# 计算热速度
mass_h = 1.008 * CHEMISTRY.u  # 氢原子质量
v_thermal = utils.thermal_velocity(mass_h, temp)
print(f"Hydrogen thermal velocity at {temp} K: {v_thermal:.0f} m/s")
```

### 2. 电化学计算

```python
# 德拜长度计算（电解质筛选长度）
temp = 298.15  # K
ionic_strength = 0.1  # mol/L
debye_length = utils.debye_length(temp, ionic_strength, 'nm')
print(f"Debye length: {debye_length:.1f} nm")
```

### 3. 光谱学计算

```python
# 从能量计算频率
photon_energy_ev = 2.0  # eV
photon_energy_j = utils.convert_energy(photon_energy_ev, 'eV', 'J')
frequency = utils.planck_frequency(photon_energy_j, 'THz')
print(f"Photon frequency: {frequency:.1f} THz")
```

## 针对Ga₂O₃研究的特殊支持

MDemon特别针对Ga₂O₃材料研究提供了丰富的材料常量：

```python
from MDemon.constants import MATERIAL

# β-Ga₂O₃（单斜晶系）
beta_density = MATERIAL.get_density('Ga2O3_beta')
beta_thermal_cond = MATERIAL.get_thermal_conductivity('Ga2O3_beta')
beta_params = MATERIAL.get_lattice_parameters('Ga2O3_beta')

print(f"β-Ga₂O₃ properties:")
print(f"  Density: {beta_density} kg/m³")
print(f"  Thermal conductivity: {beta_thermal_cond} W/(m⋅K)")
print(f"  Lattice: a={beta_params['a']:.3f} Å, β={beta_params['beta']:.1f}°")

# α-Ga₂O₃（刚玉结构）
alpha_density = MATERIAL.get_density('Ga2O3_alpha')
alpha_params = MATERIAL.get_lattice_parameters('Ga2O3_alpha')

# 比较不同氧化物
oxides = ['Ga2O3_beta', 'Al2O3', 'ZnO', 'TiO2']
for oxide in oxides:
    density = MATERIAL.get_density(oxide)
    if density:
        print(f"{oxide}: {density} kg/m³")
```

## 最佳实践

### 1. 单位一致性

```python
# 推荐：始终使用SI单位进行计算
from MDemon.constants import utils, PHYSICS

# 输入数据转换为SI单位
energy_input_ev = 2.5  # eV
energy_si = utils.convert_energy(energy_input_ev, 'eV', 'J')

length_input_ang = 5.0  # Å
length_si = utils.convert_length(length_input_ang, 'Å', 'm')

# 使用SI单位计算
result = energy_si / (PHYSICS.k_B * 300)  # 温度当量

# 结果转换为所需单位
result_output = utils.convert_temperature(result, 'K', 'C')
```

### 2. 错误处理

```python
from MDemon.constants import utils

try:
    # 尝试转换
    result = utils.convert_energy(1.0, 'unknown_unit', 'J')
except ValueError as e:
    print(f"转换错误: {e}")

# 检查材料是否存在
density = MATERIAL.get_density('unknown_material')
if density is None:
    print("材料未找到")
```

### 3. 性能考虑

```python
# 对于重复使用的常量，可以预先提取
from MDemon.constants import PHYSICS

k_B = PHYSICS.k_B
N_A = PHYSICS.N_A

# 在循环中使用预提取的常量
for temp in temperature_range:
    thermal_energy = k_B * temp  # 比每次访问PHYSICS.k_B更快
```

## 扩展指南

### 1. 添加新材料

在`material.py`中的相应字典添加新材料数据：

```python
# 在MaterialConstants类中添加
densities = {
    # ... 现有材料 ...
    'new_material': 1234.0,  # kg/m³
}

melting_points = {
    # ... 现有材料 ...
    'new_material': 1000.0,  # K
}
```

### 2. 添加新的转换单位

在`conversion.py`中添加新的转换因子：

```python
# 新的长度单位
FURLONG_TO_METER = 201.168  # 1 furlong = 201.168 m

# 在utils.py的convert_length函数中添加支持
unit_map = {
    # ... 现有单位 ...
    'furlong': CONVERSION.FURLONG_TO_METER,
}
```

## 常见问题

### Q: 如何查看所有可用的单位？
```python
from MDemon.constants import utils
print(utils.list_available_units())
```

### Q: 如何查看所有可用的材料？
```python
from MDemon.constants import MATERIAL
print(MATERIAL.list_available_materials())
```

### Q: 常量的精度如何？
所有基本物理常量基于2018年CODATA国际推荐值，提供科学计算所需的精度。

### Q: 如何在Jupyter notebook中使用？
```python
# 在notebook中导入并使用
from MDemon.constants import *
%matplotlib inline

# 示例计算和可视化
import matplotlib.pyplot as plt
import numpy as np

temperatures = np.linspace(200, 400, 100)
thermal_energies = [utils.boltzmann_energy(T, 'eV') for T in temperatures]

plt.plot(temperatures, thermal_energies)
plt.xlabel('Temperature (K)')
plt.ylabel('Thermal Energy (eV)')
plt.title('Thermal Energy vs Temperature')
plt.show()
```

## 参考资料

- [CODATA 2018 Physical Constants](https://physics.nist.gov/cuu/Constants/)
- [IUPAC Physical Chemistry Data](https://iupac.org/)
- [Ga₂O₃ Materials Database](https://www.ga2o3.org/)

---

有关更多信息，请参阅[MDemon文档](../README.md)或查看[examples/constants_usage.py](../examples/constants_usage.py)中的完整示例。 