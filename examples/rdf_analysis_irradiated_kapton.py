#!/usr/bin/env python3
"""
IrradiatedKapton RDF Analysis Example

使用MDemon的SingleAtomAnalysis.rdf方法计算IrradiatedKapton数据中各种元素之间的RDF分布并作图

Usage:
    uv run rdf_analysis_irradiated_kapton.py
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import MDemon as md
from MDemon.analysis.single_atom import SingleAtomAnalysis


def main():
    """主分析函数"""
    print("🧬 IrradiatedKapton RDF分析")
    print("=" * 50)

    # 1. 加载数据
    script_dir = Path(__file__).parent
    data_dir = script_dir.parent / "tests" / "data" / "lammps" / "IrradiatedKapton"
    base_name = "KAPTON5_504-33r_irradiated"
    data_file = data_dir / f"{base_name}.data"
    bond_file = data_dir / f"{base_name}.reaxff"

    print(f"📖 加载数据文件...")
    universe = md.Universe(str(data_file), str(bond_file))
    print(f"   原子数量: {len(universe.atoms)}")

    # 2. 分析系统组成
    print(f"\n🔍 分析系统组成...")
    species_count = {}
    for atom in universe.atoms:
        species = atom.species
        species_count[species] = species_count.get(species, 0) + 1

    for species, count in sorted(species_count.items()):
        percentage = (count / len(universe.atoms)) * 100
        print(f"   Species {species}: {count:,} atoms ({percentage:.1f}%)")

    # 3. 创建分析器并计算RDF
    print(f"\n🧮 计算RDF分布...")
    analysis = SingleAtomAnalysis(universe)

    # 对所有原子计算RDF (可以选择采样以提高计算效率)
    sample_size = min(1000, len(universe.atoms))  # 采样200个原子
    sample_indices = np.random.choice(len(universe.atoms), sample_size, replace=False)

    rdf_result = analysis.rdf(
        atom_indices=sample_indices.tolist(), r_range=(0.5, 15.0), n_bins=150
    )

    print(f"✅ RDF计算完成，分析了 {len(sample_indices)} 个原子")

    # 4. 获取所有可用的species组合并绘图
    rdf_result.plot()

    print(f"\n🎉 RDF分析完成!")


if __name__ == "__main__":
    main()
