# MDemon 文件读取架构指南

本文档详细介绍 MDemon 框架的核心文件读取逻辑和数据组织架构。

## 目录

1. [总体架构概览](#总体架构概览)
2. [核心组件详解](#核心组件详解)
3. [数据流转过程](#数据流转过程)
4. [读取器实现细节](#读取器实现细节)
5. [实际使用示例](#实际使用示例)

---

## 总体架构概览

MDemon 采用**分层架构**设计，主要包含以下层次：

```
User API (Universe)
      ↓
  Structure Layer (Atom, Molecule, Bond, etc.)
      ↓
  Attribute Layer (Coordinate, Velocity, Mass, etc.)
      ↓
  Database Layer (Data Management & Registration)
      ↓
  Reader Layer (File Format Parsing)
      ↓
  Raw Data Files (LAMMPS .data, REAXFF .reaxff, etc.)
```

### 核心设计理念

1. **元类驱动的自动注册**：使用元类（Metaclass）自动注册 Reader、Structure、StructureAttr
2. **基于 Source 的数据存储**：所有数据通过 Source 对象统一管理
3. **延迟实例化**：Structure 对象在需要时才创建
4. **结构 ID（SID）系统**：使用元组标识不同层次的结构（如 `("Atom_Base",)` 或 `("Bond_Base", "Atom_Base")`）

---

## 核心组件详解

### 1. Reader 层：文件解析

#### 1.1 Reader 注册机制

MDemon 使用**元类 `_ReaderMeta`** 实现自动注册：

```python
# MDemon/reader/base.py
class _ReaderMeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        try:
            fmt = util.asiterable(classdict["format"])
        except KeyError:
            pass
        else:
            for fmt_name in fmt:
                fmt_name = fmt_name.upper()
                _READERS[fmt_name] = cls  # 自动注册到全局字典
```

**关键点**：
- 每个 Reader 类定义 `format` 属性（如 `format = "DATA"`）
- 元类在类定义时自动将其注册到 `_READERS` 全局字典
- 支持多格式注册（format 可以是列表）

#### 1.2 ReaderBase 基类

所有 Reader 必须继承 `ReaderBase` 并实现 `parse()` 方法：

```python
class ReaderBase(metaclass=_ReaderMeta):
    def __init__(self, filename):
        self.filename = filename
    
    def parse(self, **kwargs):
        """子类必须实现：解析文件并返回 Database 对象"""
        raise NotImplementedError
    
    # 支持上下文管理器
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False
```

#### 1.3 获取合适的 Reader

通过 `get_reader_for()` 函数自动选择 Reader：

```python
# MDemon/core/rw.py
def get_reader_for(filename, format=None):
    if format is None:
        format = util.guess_format(filename)  # 根据文件扩展名猜测
    format = format.upper()
    return _READERS[format]  # 从注册表中获取
```

### 2. Database 层：数据管理中枢

`Database` 是整个数据管理的核心，负责：
- 管理所有 `StructureAttr`（结构属性）
- 维护 `_source_register`：存储所有属性的数据源
- 协调不同 Structure 之间的关系

#### 2.1 Database 初始化

```python
class Database:
    def __init__(self, n_atoms, **kwargs):
        self._source_register = {}      # 存储所有属性数据
        self._deep_source_register = {} # 时间依赖数据
        self.attrs = []                 # 属性列表
        self.attrnames = []             # 属性名称列表
        
        self._n_atoms = n_atoms
        
        # 初始化索引
        Index(n_atoms, sid=("Atom_Base",), database=self)
        
        # 创建 Family（家族）系统
        self._base = Family(database=self)
        self._families = {"Base": self._base}
```

**关键概念 - Family（家族）**：
- MDemon 使用 "家族" 概念组织不同层次的结构
- 每个家族包含一组相关的 Structure 类（Atom, Molecule, Bond 等）
- 默认有一个 "Base" 家族

#### 2.2 属性注册流程

```python
def add_Attr(self, attr):
    """添加新属性到 Database"""
    self.attrs.append(attr)
    self.attrnames.append(attr.name)
    attr._database = self
    self.__setattr__(attr.name, attr)  # 作为 Database 的属性

def register_source(self, attr):
    """为属性注册数据源"""
    try:
        attr._source_register = self._source_register[attr.name]
    except KeyError:
        self._source_register[attr.name] = {}
        attr._source_register = self._source_register[attr.name]
```

### 3. StructureAttr 层：属性系统

`StructureAttr` 是所有结构属性的基类，包括坐标、速度、质量、电荷等。

#### 3.1 StructureAttr 类层次

```
StructureAttr (抽象基类)
    ├── StructureAttr1D (一维属性: 单一 SID)
    │   ├── Coordinate (坐标)
    │   ├── Velocity (速度)
    │   ├── ID (标识符)
    │   ├── Species (物种类型)
    │   └── ParticleAttr (粒子属性)
    │       ├── Mass (质量)
    │       ├── Charge (电荷)
    │       └── Temperature (温度)
    │
    └── StructureAttr2D (二维属性: 两个 SID)
        └── Composition (组成关系，如分子-原子)
```

#### 3.2 SID（Structure ID）系统

SID 是一个**元组**，用于标识数据的结构层次：

- **一维 SID**：`("Atom_Base",)` - 原子级别属性
- **二维 SID**：`("Molecule_Base", "Atom_Base")` - 分子-原子关系

示例：
```python
# 原子坐标：每个原子一个坐标向量
Coordinate(coords, sid=("Atom_Base",), database=dbase)

# 分子组成：哪些原子属于哪个分子
Composition(matrix, sid=("Molecule_Base", "Atom_Base"), 
            database=dbase, N=n_molecules, M=n_atoms)
```

#### 3.3 Source 对象：数据存储

属性数据通过 `Source` 对象存储（在 `source.pyx` 中实现）：

```python
class Source1D:
    """一维数据源（如坐标、速度）"""
    def __init__(self, dtype, values):
        self._dtype = dtype  # 'int', 'float', 'bool'
        self._values = np.asarray(values)
    
    @property
    def values(self):
        return self._values
    
    @values.setter
    def values(self, valix):
        """valix = [indices, new_values]"""
        ix, val = valix
        if ix is None:
            self._values = val
        else:
            self._values[ix] = val

class Source2D:
    """二维数据源（如组成关系矩阵）"""
    # 使用稀疏矩阵存储大规模关系数据
```

### 4. Structure 层：结构对象

Structure 是用户直接交互的对象，包括：

#### 4.1 Structure 类层次

```
Structure (基类)
    ├── Particle (粒子类)
    │   ├── Atom (原子)
    │   ├── Molecule (分子)
    │   └── Group (原子组)
    │
    └── Topology (拓扑类)
        ├── Chain (链)
        │   ├── Bond (键)
        │   ├── Angle (角)
        │   └── Dihedral (二面角)
        ├── Tree (树)
        │   └── Improper (不当二面角)
        └── Ring (环)
            └── Multiring (多环)
```

#### 4.2 Structure 动态类生成

MDemon 使用**混入（Mixin）机制**动态生成 Structure 类：

```python
# MDemon/core/universe.py
def _make_bases():
    """创建 Structure 基类"""
    bases = {}
    SBase = bases[Structure] = _StructureAttrContainer._subclass()
    
    # 为 Particle 创建子类
    PBase = SBase._subclass()
    for cls in _PARTICLES:
        bases[cls] = PBase._subclass()
    
    # 为 Topology 创建子类
    TBase = SBase._subclass()
    for cls in _TOPOLOGIES:
        bases[cls] = TBase._subclass()
    
    return bases
```

**为什么要动态生成？**
- 允许在运行时添加属性（如 `atom.coordinate`、`atom.velocity`）
- 每个 Universe 有自己的类实例，避免冲突
- 支持属性的动态绑定

#### 4.3 属性绑定机制

StructureAttr 通过 **property** 动态绑定到 Structure：

```python
@classmethod
def _add_prop(cls, attr):
    """将属性添加到 Structure 类"""
    sids = attr.import_sids(cls)
    
    for sid in sids:
        attrname = attr.import_attrname(sid)
        
        def make_getter(sid):
            def getter(self):
                return attr.__getitem__(self, sid)
            return getter
        
        def make_setter(sid):
            def setter(self, values):
                return attr.__setitem__(self, sid, values)
            return setter
        
        # 动态添加 property
        setattr(cls, attrname, property(make_getter(sid), make_setter(sid)))
```

这样，用户可以直接访问：
```python
atom = universe.atoms[0]
coords = atom.coordinate  # 调用 Coordinate.__getitem__(atom, ("Atom_Base",))
atom.coordinate = new_coords  # 调用 Coordinate.__setitem__()
```

### 5. Universe 层：顶层接口

`Universe` 是用户的主要入口：

```python
class Universe:
    def __init__(self, *inputfiles, **kwargs):
        self.timestep = 0
        self._class_bases = _make_bases()  # 创建动态类
        
        # 调用 Reader 解析文件
        self._database = _database_from_file_like(*inputfiles, **kwargs)
        self._database._u = self  # 双向引用
        
        # 生成 Structure 对象
        _generate_from_database(self._database)
```

---

## 数据流转过程

### 完整的文件读取流程

以读取 LAMMPS DATA 文件为例：

```python
# 用户代码
u = md.Universe("system.data", "bonds.reaxff")
```

**流程详解**：

#### 阶段 1：Universe 初始化

```python
# MDemon/core/universe.py
def __init__(self, *inputfiles, **kwargs):
    # 1. 创建动态类基础
    self._class_bases = _make_bases()
    
    # 2. 从文件创建 Database
    self._database = _database_from_file_like(*inputfiles, **kwargs)
    
    # 3. 建立双向引用
    self._database._u = self
    
    # 4. 生成所有 Structure 实例
    _generate_from_database(self._database)
```

#### 阶段 2：文件解析

```python
def _database_from_file_like(*inputfiles, **kwargs):
    for file_ in inputfiles:
        # 1. 获取合适的 Reader
        reader = get_reader_for(file_)  # 返回 DATAReader 或 REAXFFReader
        
        # 2. 使用上下文管理器打开文件
        with reader(file_) as r:
            # 3. 解析文件，返回 Database
            database = r.parse(**kwargs)
            
            # 4. 将 database 传递给下一个文件的 Reader
            kwargs["database"] = database
    
    return database
```

#### 阶段 3：LAMMPS DATA 文件解析

```python
# MDemon/reader/LAMMPS.py
class DATAReader(ReaderBase):
    format = "DATA"
    
    def parse(self, **kwargs):
        # 1. 解析文件结构
        head, sects = self.grab_datafile()
        
        # 2. 解析质量表
        masses = self._parse_masses(sects["Masses"])
        
        # 3. 解析原子数据 - 核心步骤！
        ids, species, compo, dbase, order, masses = self._parse_atoms(
            sects["Atoms"], masses
        )
        
        # 4. 解析速度
        if "Velocities" in sects:
            velocities = self._parse_vel(sects["Velocities"], order)
        else:
            velocities = np.zeros((n_atoms, 3))
        
        # 5. 创建速度和温度属性
        Velocity(velocities, sid=("Atom_Base",), database=dbase)
        Temperature.start_from_velocities(velocities, masses, 
                                         sid=("Atom_Base",), database=dbase)
        
        # 6. 解析拓扑信息（键、角、二面角）
        for sname, L, nentries in [
            ("Bond_Base", "Bonds", 2),
            ("Angle_Base", "Angles", 3),
            ("Dihedral_Base", "Dihedrals", 4),
        ]:
            if L in sects:
                # 解析并创建组成关系
                id_, type, sect = self._parse_bond_section(sects[L], nentries, mapping)
                # ...更新 database
        
        # 7. 创建盒子信息
        Box(self._parse_box(head), database=dbase)
        
        return dbase
```

#### 阶段 4：原子数据解析（最重要）

```python
def _parse_atoms(self, datalines, massdict=None):
    n_atoms = len(datalines)
    
    # 1. 识别原子样式（atomic, molecular, full 等）
    if self.style_dict is None:
        n = len(datalines[0].split())
        if n in (5, 8):      # atomic: id type x y z
            sd = {"id": 0, "type": 1, "coord": 2}
        elif n in (6, 9):    # molecular: id resid type x y z
            sd = {"id": 0, "resid": 1, "type": 2, "coord": 3}
        elif n in (7, 10):   # full: id resid type charge x y z
            sd = {"id": 0, "resid": 1, "type": 2, "charge": 3, "coord": 4}
    
    # 2. 读取原始数据
    atom_ids = np.zeros(n_atoms, dtype=np.int32)
    types = np.zeros(n_atoms, dtype=object)
    coords = np.zeros((n_atoms, 3), dtype=np.float32)
    
    for i, line in enumerate(datalines):
        line = line.split()
        atom_ids[i] = line[sd["id"]]
        types[i] = line[sd["type"]]
        coords[i] = np.array(line[sd["coord"]:sd["coord"]+3])
        if has_charge:
            charges[i] = line[sd["charge"]]
    
    # 3. 排序（LAMMPS 原子可能无序）
    order = np.argsort(atom_ids)
    atom_ids = atom_ids[order]
    types = types[order]
    coords = coords[order]
    
    # 4. 创建 Database
    dbase = Database(n_atoms)
    
    # 5. 创建各种属性（关键！）
    atm = "Atom_Base"
    mle = "Molecule_Base"
    
    # 坐标
    Coordinate(coords, sid=atm, database=dbase)
    
    # 物种类型
    species = Species(types, sid=atm, database=dbase)
    
    # 电荷（如果有）
    if has_charge:
        Charge(charges, sid=atm, database=dbase)
    
    # 质量
    masses = np.array([massdict[t] for t in types])
    Mass(masses, sid=atm, database=dbase)
    
    # 元素（从质量推断）
    Element.start_from_masses(masses, sid=atm, database=dbase)
    
    # ID
    ids = ID(atom_ids, sid=(atm,), database=dbase)
    
    # 6. 创建分子层次
    residx, resids = squash_by(resids)[:2]
    ids._update_source(resids, (mle,))
    dbase.ix._update_source(np.arange(len(resids)), (mle,))
    
    # 7. 创建分子-原子组成关系
    n_residues = len(resids)
    compomtrx = self.residx2mtrx(residx)
    compo = Composition(compomtrx, sid=(mle, atm), 
                       database=dbase, N=n_residues, M=n_atoms)
    
    return ids, species, compo, dbase, order, masses
```

**关键技术点**：
1. **自动样式识别**：根据字段数量识别 LAMMPS 原子样式
2. **数据排序**：LAMMPS 文件中原子可能无序，需要排序
3. **层次结构**：同时创建原子层（Atom_Base）和分子层（Molecule_Base）
4. **属性创建**：每个属性都是独立的对象，通过 SID 关联
5. **组成矩阵**：使用稀疏矩阵存储原子-分子关系

#### 阶段 5：REAXFF 键信息补充

```python
# MDemon/reader/REAXFF.py
class REAXFFReader(ReaderBase):
    format = "REAXFF"
    
    def parse(self, **kwargs):
        # 1. 获取已存在的 database（从 DATA 文件创建）
        dbase = kwargs["database"]
        n_atoms = dbase.n_atoms
        
        # 2. 读取 REAXFF 数据
        datalines = list(self.iterdata())
        
        charges = np.zeros(n_atoms, dtype=np.float32)
        lonepairs = np.zeros(n_atoms, dtype=np.float32)
        
        for i, line in enumerate(datalines):
            line = line.split()
            atom_ids[i] = line[0]
            nb = np.int32(line[2])  # 键数量
            charges[i] = line[6 + nb * 2]
            lonepairs[i] = line[5 + nb * 2]
        
        # 3. 更新电荷（覆盖 DATA 文件中的值）
        dbase.charge._update_source(charges, ("Atom_Base",))
        
        # 4. 添加孤对电子属性
        LonePair(lonepairs, sid=("Atom_Base",), database=dbase)
        
        # 5. 创建键和键级
        bond_ids = np.arange(1, n_bonds + 1)
        bos = np.zeros(n_bonds, dtype=np.float32)  # 键级
        
        # 解析键连接信息...
        
        # 6. 创建键级属性
        BondOrder(bos, sid=("Bond_Base",), database=dbase)
        
        # 7. 创建连接矩阵
        Connection(np.array([row, col, data]), 
                  sid=("Atom_Base", "Atom_Base"),
                  database=dbase)
        
        # 8. 计算价态
        Valence.start_from_bondorders(bos, database=dbase)
        
        return dbase
```

#### 阶段 6：生成 Structure 实例

```python
def _generate_from_database(dbase):
    """从 Database 生成 Structure 对象"""
    family = dbase.base
    
    # 1. 创建动态类（混入 Universe 特定的类）
    _make_classes(family)
    
    # 2. 将属性绑定到 Structure 类
    for attr in dbase.attrs:
        if attr.__class__ not in _UNIVERSE_ATTRS:
            family._process_attr(attr)  # 调用 _add_prop
        else:
            dbase._u._add_prop(attr)    # Universe 级别属性
    
    # 3. 实例化所有 Structure 对象
    family.instancing()
```

```python
class Family:
    def instancing(self):
        """实例化所有注册的 Structure"""
        for cls in _STRUCTURES:  # [Atom, Molecule, Bond, ...]
            cls1 = self._u.families[self.name]._classes[cls]
            sid = (cls1.import_sname(),)  # 如 ("Atom_Base",)
            attrname = cls1.__name__.lower() + "s"  # "atoms"
            
            index = self._database.ix
            if sid in index._source_register:
                n = len(index._source_register[sid].values)
                s = cls1(np.arange(n), self._u)  # 创建 Structure 实例
                
                # 添加到 Family 和 Universe
                self.__setattr__(attrname, s)
                if self.name == "Base":
                    self._u.__setattr__(attrname, s)
```

**完成后的状态**：

```python
u.atoms       # Atom([0, 1, 2, ..., n-1], u)
u.molecules   # Molecule([0, 1, 2, ..., m-1], u)
u.bonds       # Bond([0, 1, 2, ..., b-1], u)
```

每个对象都可以访问属性：

```python
u.atoms[0].coordinate  # np.array([x, y, z])
u.atoms[0].velocity    # np.array([vx, vy, vz])
u.atoms[0].mass        # float
u.molecules[0].atoms   # Atom([...])  分子中的原子
```

---

## 读取器实现细节

### LAMMPS DATA Reader

#### 文件结构

LAMMPS DATA 文件格式：
```
# Comment line

5000 atoms
100 bonds
...

10 atom types
...

0.0 50.0 xlo xhi
0.0 50.0 ylo yhi
0.0 50.0 zlo zhi

Masses

1 12.011    # C
2 1.008     # H
...

Atoms  # full

1 1 1 0.0 1.234 2.345 3.456
2 1 1 0.0 1.345 2.456 3.567
...

Velocities

1 0.001 0.002 0.003
2 0.004 0.005 0.006
...

Bonds

1 1 1 2
2 1 2 3
...
```

#### 解析流程

```python
def grab_datafile(self):
    """分割 DATA 文件为头部和各个段落"""
    f = list(self.iterdata())  # 去除注释和空行
    
    # 找到所有段落起始位置
    starts = [i for i, line in enumerate(f) 
              if line.split()[0] in SECTIONS]
    
    # 解析头部信息
    header = {}
    for line in f[:starts[0]]:
        for token in HEADERS:
            if line.endswith(token):
                header[token] = line.split(token)[0]
    
    # 分割各个段落
    sects = {f[l]: f[l+1:starts[i+1]] 
             for i, l in enumerate(starts[:-1])}
    
    return header, sects
```

#### 原子样式自动识别

```python
# 支持的原子样式：
# atomic: 5 fields    - id type x y z
# atomic: 8 fields    - id type x y z flag1 flag2 flag3
# molecular: 6 fields - id resid type x y z
# molecular: 9 fields - id resid type x y z flag1 flag2 flag3
# full: 7 fields      - id resid type charge x y z
# full: 10 fields     - id resid type charge x y z flag1 flag2 flag3

def _parse_atoms(self, datalines, massdict=None):
    n = len(datalines[0].split())
    if n in (5, 8):
        sd = {"id": 0, "type": 1, "coord": 2}
    elif n in (6, 9):
        sd = {"id": 0, "resid": 1, "type": 2, "coord": 3}
    elif n in (7, 10):
        sd = {"id": 0, "resid": 1, "type": 2, "charge": 3, "coord": 4}
```

**Flag 字段处理**：
- LAMMPS 可以在坐标后添加 3 个 flag 字段（如冻结状态）
- Reader 自动识别并跳过这些字段
- 原子样式可以通过 `atom_style` 参数手动指定

### REAXFF Bond Reader

#### 文件格式

REAXFF 键文件格式：
```
1 1 3 2 0.95 3 1.02 4 0.88 0.123 -0.456
2 1 2 1 1.00 3 0.95 0.234 -0.567
...
```

字段说明：
```
atom_id type n_bonds [neighbor1 bo1 neighbor2 bo2 ...] nlp charge
```
- `atom_id`: 原子 ID
- `type`: 原子类型
- `n_bonds`: 键数量
- `neighbor_i`: 邻居原子 ID
- `bo_i`: 键级（bond order）
- `nlp`: 孤对电子数
- `charge`: 电荷

#### 解析流程

```python
def parse(self, **kwargs):
    dbase = kwargs["database"]  # 获取已存在的 database
    n_atoms = dbase.n_atoms
    
    # 1. 第一遍扫描：读取电荷、孤对电子
    for i, line in enumerate(datalines):
        line = line.split()
        nb = np.int32(line[2])
        charges[i] = line[6 + nb * 2]      # 电荷位于键信息之后
        lonepairs[i] = line[5 + nb * 2]    # 孤对电子
    
    # 2. 第二遍扫描：构建键连接
    combi_sets = {}  # 避免重复键
    for i, line in enumerate(datalines):
        line = line.split()
        nb = np.int32(line[2])
        a1_id = atom_ids[i]
        
        for j in range(nb):
            a2_id = np.int32(line[3 + j])
            combi = frozenset((a1_id, a2_id))
            
            if combi not in combi_sets:
                # 新键
                combi_sets[combi] = n
                bos[n] = line[4 + nb + j]  # 键级
                n += 1
```

**去重机制**：
- REAXFF 文件中每个键出现两次（A-B 和 B-A）
- 使用 `frozenset` 确保每个键只记录一次
- 键级从第一次出现中获取

---

## 实际使用示例

### 示例 1：基本读取

```python
import MDemon as md

# 读取 LAMMPS DATA 文件
u = md.Universe("system.data")

# 访问原子
print(f"Total atoms: {len(u.atoms)}")
atom = u.atoms[0]
print(f"Atom 0 coordinate: {atom.coordinate}")
print(f"Atom 0 velocity: {atom.velocity}")

# 访问分子
print(f"Total molecules: {len(u.molecules)}")
mol = u.molecules[0]
print(f"Molecule 0 has {len(mol.atoms)} atoms")

# 访问键
if hasattr(u, 'bonds'):
    print(f"Total bonds: {len(u.bonds)}")
```

### 示例 2：读取多个文件

```python
# 先读取 DATA，再读取 REAXFF 补充键信息
u = md.Universe("system.data", "bonds.reaxff")

# REAXFF 文件提供了额外信息
atom = u.atoms[0]
print(f"Charge: {atom.charge}")           # 来自 REAXFF
print(f"Lone pairs: {atom.lonepair}")     # REAXFF 特有

# 键级信息
bond = u.bonds[0]
print(f"Bond order: {bond.bondorder}")    # REAXFF 特有
```

### 示例 3：手动指定原子样式

```python
# 对于非标准格式，手动指定原子样式
u = md.Universe("custom.data", 
                atom_style="id type x y z vx vy vz")
```

### 示例 4：Structure 切片和过滤

```python
# 选择特定原子
carbon_atoms = u.atoms[[i for i, a in enumerate(u.atoms) 
                       if a.element == 6]]  # 碳原子

# 获取第一个分子的坐标
mol_coords = u.molecules[0].coordinate

# 选择特定键
single_bonds = u.bonds[[i for i, b in enumerate(u.bonds)
                       if b.bondorder < 1.5]]
```

### 示例 5：数据分析应用

```python
from MDemon.analysis.single_atom import SingleAtomAnalysis

# 创建分析器
analysis = SingleAtomAnalysis(u)

# 计算 RDF
rdf_result = analysis.rdf(
    atom_indices=range(100),  # 采样前 100 个原子
    r_range=(0.5, 15.0),
    n_bins=150
)

# 绘图
rdf_result.plot()
```

---

## 架构优势

### 1. 灵活性

- **多格式支持**：通过 Reader 注册机制轻松添加新格式
- **动态属性**：Structure 属性可在运行时添加
- **可扩展性**：元类驱动的设计便于扩展

### 2. 性能

- **延迟计算**：Structure 只在访问时创建
- **稀疏存储**：大规模关系（如键）使用稀疏矩阵
- **Cython 加速**：核心 Source 对象用 Cython 实现

### 3. 可维护性

- **分层设计**：各层职责清晰
- **统一接口**：所有 Reader 遵循相同模式
- **类型安全**：使用 dtype 确保数据类型一致

### 4. 用户友好

- **直观 API**：`atom.coordinate` 而非 `get_coordinate(atom)`
- **自动推断**：文件格式和原子样式自动识别
- **上下文管理**：自动文件关闭

---

## 扩展开发指南

### 添加新的 Reader

```python
from MDemon.reader.base import ReaderBase
from MDemon.core.database import Database
from MDemon.core.structureattr import Coordinate, Species

class MyFormatReader(ReaderBase):
    format = "MYFORMAT"  # 自动注册
    
    def parse(self, **kwargs):
        # 1. 读取文件
        with open(self.filename) as f:
            data = self._parse_file(f)
        
        # 2. 创建 Database
        n_atoms = len(data['atoms'])
        dbase = Database(n_atoms)
        
        # 3. 添加属性
        Coordinate(data['coords'], sid=("Atom_Base",), database=dbase)
        Species(data['types'], sid=("Atom_Base",), database=dbase)
        
        # 4. 返回 Database
        return dbase
    
    def _parse_file(self, f):
        # 实现具体解析逻辑
        pass
```

### 添加新的 StructureAttr

```python
from MDemon.core.structureattr import StructureAttr1D

class MyAttr(StructureAttr1D):
    name = "myattr"
    _dtype = "float"
    _sid0 = ("Atom_Base",)
    
    # 继承默认的 _update_source 和 _parse_source
```

使用：
```python
# 在 Reader 中
MyAttr(values, sid=("Atom_Base",), database=dbase)

# 访问
atom.myattr  # 自动绑定！
```

---

## 总结

MDemon 的文件读取架构体现了以下设计原则：

1. **元编程**：利用元类实现自动注册和动态类生成
2. **分层抽象**：Reader → Database → StructureAttr → Structure → Universe
3. **数据驱动**：通过 SID 系统统一管理不同层次的数据
4. **延迟求值**：按需创建对象，提高性能
5. **可扩展性**：新格式和新属性的添加遵循相同模式

这种架构使得 MDemon 既强大又易于扩展，适合处理复杂的分子动力学数据。

---

**文档版本**：v1.0  
**最后更新**：2025-10-25  
**作者**：MDemon Development Team

