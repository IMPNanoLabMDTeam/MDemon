#!/usr/bin/env python3
"""
XYZ Reader 使用示例

演示如何使用 MDemon 读取和分析 Extended XYZ 格式文件。

Usage:
    uv run python examples/xyz_reader_example.py
"""

from pathlib import Path

import numpy as np

import MDemon as md


def main():
    """主函数：演示 XYZ reader 功能"""

    print("=" * 70)
    print("MDemon Extended XYZ Reader 示例")
    print("=" * 70)

    # 1. 读取 XYZ 文件
    print("\n📂 步骤 1: 读取 Extended XYZ 文件")
    print("-" * 70)

    data_dir = Path(__file__).parent.parent / "tests" / "data" / "xyz" / "basic"
    xyz_file = data_dir / "beta_221_ort.xyz"

    if not xyz_file.exists():
        print(f"❌ 文件不存在: {xyz_file}")
        return

    print(f"文件路径: {xyz_file}")
    u = md.Universe(str(xyz_file))
    print(f"✅ 成功加载 Universe")

    # 2. 基本信息
    print("\n📊 步骤 2: 系统基本信息")
    print("-" * 70)

    print(f"原子总数: {len(u.atoms)}")
    print(f"分子数量: {len(u.molecules)}")

    # 3. 晶格信息
    print("\n📦 步骤 3: 晶格信息")
    print("-" * 70)

    box = u.box
    print(f"Box 类型: {'正交晶格' if len(box.shape) == 1 else '三斜晶格'}")

    if len(box.shape) == 1 and len(box) == 12:
        # 正交盒子
        lx, ly, lz = box[0], box[1], box[2]
        print(f"晶格参数:")
        print(f"  a = {lx:.6f} Å")
        print(f"  b = {ly:.6f} Å")
        print(f"  c = {lz:.6f} Å")
        print(f"  α = β = γ = 90°")

        volume = lx * ly * lz
        print(f"\n体积: {volume:.3f} Å³")
    else:
        # 三斜盒子
        print(f"晶格向量矩阵:\n{box}")

    # 4. 元素组成分析
    print("\n🔬 步骤 4: 元素组成分析")
    print("-" * 70)

    species_list = [atom.species for atom in u.atoms]
    unique_species = sorted(set(species_list))

    print(f"唯一物种数: {len(unique_species)}")
    print(f"\n物种分布:")

    species_counts = {}
    for s in species_list:
        species_counts[s] = species_counts.get(s, 0) + 1

    for species in sorted(species_counts.keys()):
        count = species_counts[species]
        percentage = (count / len(u.atoms)) * 100
        print(f"  Species {species}: {count:3d} 个原子 ({percentage:5.2f}%)")

    # 5. 质量和密度计算
    print("\n⚖️  步骤 5: 质量和密度计算")
    print("-" * 70)

    masses = [atom.mass for atom in u.atoms]
    total_mass = sum(masses)

    print(f"总质量: {total_mass:.3f} amu")
    print(f"平均原子质量: {np.mean(masses):.3f} amu")

    if len(box.shape) == 1:
        volume = box[0] * box[1] * box[2]  # Å³
        # 1 amu/Å³ = 1.66054 g/cm³
        density = total_mass / volume * 1.66054
        print(f"密度: {density:.3f} g/cm³")

    # 6. 坐标信息
    print("\n📍 步骤 6: 坐标信息")
    print("-" * 70)

    coords = np.array([atom.coordinate for atom in u.atoms])

    print(f"坐标范围:")
    print(f"  X: [{coords[:, 0].min():.3f}, {coords[:, 0].max():.3f}] Å")
    print(f"  Y: [{coords[:, 1].min():.3f}, {coords[:, 1].max():.3f}] Å")
    print(f"  Z: [{coords[:, 2].min():.3f}, {coords[:, 2].max():.3f}] Å")

    print(f"\n坐标中心:")
    center = coords.mean(axis=0)
    print(f"  ({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}) Å")

    # 7. 示例原子详细信息
    print("\n🔍 步骤 7: 示例原子详细信息")
    print("-" * 70)

    for i in [0, 10, 20]:
        if i < len(u.atoms):
            atom = u.atoms[i]
            print(f"\n原子 {i}:")
            print(f"  Species: {atom.species}")
            print(
                f"  Coordinate: ({atom.coordinate[0]:.6f}, {atom.coordinate[1]:.6f}, {atom.coordinate[2]:.6f}) Å"
            )
            print(f"  Mass: {atom.mass:.3f} amu")

            # 如果有速度信息
            if hasattr(atom, "velocity"):
                vel = atom.velocity
                vel_mag = np.linalg.norm(vel)
                print(f"  Velocity: ({vel[0]:.6f}, {vel[1]:.6f}, {vel[2]:.6f}) Å/fs")
                print(f"  |v|: {vel_mag:.6f} Å/fs")

    # 8. 按物种分类
    print("\n📑 步骤 8: 按物种分组")
    print("-" * 70)

    species_groups = {}
    for i, atom in enumerate(u.atoms):
        species = atom.species
        if species not in species_groups:
            species_groups[species] = []
        species_groups[species].append(i)

    for species in sorted(species_groups.keys()):
        indices = species_groups[species]
        print(f"\nSpecies {species}:")
        print(f"  原子数量: {len(indices)}")
        print(f"  原子索引: {indices[:5]}{'...' if len(indices) > 5 else ''}")

        # 计算该物种的平均坐标
        species_coords = coords[indices]
        species_center = species_coords.mean(axis=0)
        print(
            f"  几何中心: ({species_center[0]:.3f}, {species_center[1]:.3f}, {species_center[2]:.3f}) Å"
        )

    # 9. 统计摘要
    print("\n📈 步骤 9: 统计摘要")
    print("-" * 70)

    print(f"✓ 文件格式: Extended XYZ")
    print(f"✓ 原子总数: {len(u.atoms)}")
    print(f"✓ 物种类型: {len(unique_species)}")
    print(
        f"✓ 包含速度: {'是' if hasattr(u.atoms[0], 'velocity') and not np.allclose(u.atoms[0].velocity, 0) else '否'}"
    )
    print(f"✓ 包含电荷: {'是' if hasattr(u.atoms[0], 'charge') else '否'}")

    print("\n" + "=" * 70)
    print("✅ XYZ Reader 示例完成！")
    print("=" * 70)


def demonstrate_writer():
    """演示 XYZ Writer 功能（如果可用）"""
    print("\n" + "=" * 70)
    print("XYZ Writer 功能演示")
    print("=" * 70)

    try:
        from MDemon.reader.XYZ import XYZWriter

        # 读取文件
        data_dir = Path(__file__).parent.parent / "tests" / "data" / "xyz" / "basic"
        xyz_file = data_dir / "beta_221_ort.xyz"
        u = md.Universe(str(xyz_file))

        # 写入新文件
        output_file = Path(__file__).parent / "output.xyz"
        writer = XYZWriter(str(output_file))
        writer.write(u)

        print(f"✅ 成功写入文件: {output_file}")

        # 验证写入的文件
        u2 = md.Universe(str(output_file))
        print(f"✓ 验证: 原始原子数 = {len(u.atoms)}, 写入后原子数 = {len(u2.atoms)}")

        # 清理
        if output_file.exists():
            output_file.unlink()
            print(f"✓ 清理临时文件")

    except ImportError:
        print("XYZWriter 功能尚未实现")
    except Exception as e:
        print(f"❌ 错误: {e}")


if __name__ == "__main__":
    main()
    # demonstrate_writer()  # 取消注释以测试 Writer 功能
