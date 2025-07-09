# MDemon项目迁移到uv指南

本指南将帮助你将MDemon项目从conda环境管理迁移到uv。

## 1. 安装uv

首先安装uv（如果还没有安装）：

```bash
# macOS/Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Windows (PowerShell)
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"

# 或使用pip
pip install uv
```

## 2. 迁移步骤

### 步骤1：创建uv项目

```bash
# 在项目根目录下
cd /path/to/MDemon
uv init --no-readme  # 不覆盖现有README
```

### 步骤2：同步依赖

由于我们已经创建了`pyproject.toml`，现在安装所有依赖：

```bash
# 安装基础依赖 + 开发依赖
uv sync --extra dev

# 或者只安装基础依赖
uv sync
```

### 步骤3：构建Cython扩展

```bash
# 安装构建依赖并编译Cython扩展
uv run pip install -e . --no-build-isolation

# 或者使用更现代的方式
uv pip install -e .
```

### 步骤4：验证安装

```bash
# 运行测试验证一切正常
uv run pytest

# 检查Cython扩展是否正确编译
uv run python -c "import MDemon.lib.c_distances; print('Cython扩展导入成功！')"
```

### 步骤5：更新pre-commit hooks

```bash
# 重新安装pre-commit hooks
uv run pre-commit install

# 运行一次检查确保配置正确
uv run pre-commit run --all-files
```

## 3. 常用命令对照表

| 操作 | Conda命令 | uv命令 |
|------|-----------|--------|
| 创建环境 | `conda env create -f environment.yml` | `uv sync` |
| 激活环境 | `conda activate mdemon` | 无需激活，使用`uv run` |
| 安装包 | `conda install package` | `uv add package` |
| 运行脚本 | `conda run -n mdemon python script.py` | `uv run python script.py` |
| 运行测试 | `conda run -n mdemon pytest` | `uv run pytest` |
| 安装开发依赖 | `conda env update -f environment.yml` | `uv sync --extra dev` |

## 4. 开发工作流

### 日常开发

```bash
# 添加新的依赖
uv add numpy scipy

# 添加开发依赖
uv add --dev pytest black

# 运行脚本
uv run python your_script.py

# 运行测试
uv run pytest

# 格式化代码
uv run black .

# 代码检查
uv run ruff check .
```

### 构建和安装

```bash
# 开发安装（可编辑模式）
uv pip install -e .

# 重新编译Cython扩展（如果修改了.pyx文件）
uv run python setup.py build_ext --inplace

# 或者完全重新安装
uv pip install -e . --force-reinstall --no-deps
```

### 分发和打包

```bash
# 构建wheel包
uv run python -m build

# 构建源码包
uv run python -m build --sdist
```

## 5. 环境管理

### 查看环境信息

```bash
# 查看当前环境路径
uv python dir

# 查看已安装的包
uv pip list

# 查看项目配置
uv tree
```

### 清理和重置

```bash
# 删除虚拟环境重新创建
rm -rf .venv
uv sync

# 清理构建文件
rm -rf build/ dist/ *.egg-info/
find . -name "*.c" -delete  # 删除Cython生成的C文件
find . -name "*.so" -delete  # 删除编译的扩展
```

## 6. 故障排除

### Cython编译问题

如果Cython扩展编译失败：

```bash
# 确保有C编译器
# macOS: xcode-select --install
# Ubuntu: apt install build-essential
# Windows: 安装Microsoft C++ Build Tools

# 手动编译Cython
uv run cython MDemon/lib/c_distances.pyx
uv run cython MDemon/core/source.pyx

# 重新安装
uv pip install -e . --force-reinstall
```

### NumPy相关问题

```bash
# 如果遇到NumPy版本问题
uv add "numpy>=1.24.0,<2.0.0"
uv sync
```

### 依赖冲突

```bash
# 查看依赖树
uv tree

# 解决冲突
uv lock --upgrade
uv sync
```

## 7. 性能对比

使用uv的优势：
- 🚀 **更快的依赖解析**：比conda快10-100倍
- 💾 **更小的环境**：不包含conda的开销
- 🔄 **更好的锁定文件**：uv.lock确保可重现构建
- 🛠️ **现代工具链**：与最新Python生态系统集成更好

## 8. 迁移后清理

完成迁移并确认一切正常后，可以清理旧文件：

```bash
# 备份conda环境信息（可选）
conda env export -n mdemon > conda_backup.yml

# 删除conda环境
conda env remove -n mdemon

# 可以考虑移除但不建议立即删除
# rm environment.yml
# rm setup.cfg  # 如果所有配置都移到了pyproject.toml
```

## 9. 团队协作

确保团队成员都了解新的工作流：

1. 分享这个迁移指南
2. 更新CI/CD脚本使用uv
3. 在项目README中更新安装说明
4. 考虑在`.gitignore`中添加`.venv/`（如果还没有）

这样就完成了从conda到uv的完整迁移！ 