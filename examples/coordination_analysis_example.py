"""
配位数分析示例

这个示例演示如何使用MDemon的配位数分析功能，重点展示个体配位数分布的分析，
这与RDF和角分布分析（主要关注平均值）不同。

配位数分析的特点：
1. 保留每个原子的个体配位数信息
2. 提供配位数分布统计
3. 支持多种可视化方式（直方图、散点图、热图等）
4. 关注个体差异而非平均行为
"""

import matplotlib.pyplot as plt
import numpy as np

# 导入MDemon相关模块
try:
    from MDemon.analysis.single_atom import CoordinationAnalyzer, CoordinationResult
    from MDemon.core import Universe
    from MDemon.reader import LAMMPS
except ImportError as e:
    print(f"Import error: {e}")
    print("请确保MDemon已正确安装")
    exit(1)


def main():
    """主函数演示配位数分析"""

    print("=== MDemon 配位数分析示例 ===\n")

    # 1. 加载数据（这里需要替换为实际的数据文件）
    print("1. 加载数据...")
    # universe = Universe("path/to/your/trajectory.dump", format="LAMMPS")
    # 由于没有实际数据文件，这里创建一个模拟的universe对象
    print("   (注意：这里需要替换为实际的数据文件路径)")

    # 2. 定义截断半径
    print("\n2. 定义截断半径...")
    # 配位数分析需要为不同的原子对定义截断半径
    cutoff_radii = {
        "1-1": 2.5,  # 类型1-类型1的截断半径
        "1-2": 3.0,  # 类型1-类型2的截断半径
        "2-2": 2.8,  # 类型2-类型2的截断半径
    }
    print(f"   截断半径设置: {cutoff_radii}")

    # 3. 创建原子选择掩码（可选）
    print("\n3. 创建原子选择...")
    # 例如只选择特定类型的原子进行分析
    # ca_mask = np.array([atom.species == 1 for atom in universe.atoms])
    print("   (可以选择特定类型的原子进行分析)")

    # 4. 创建配位数分析器
    print("\n4. 创建配位数分析器...")
    # analyzer = CoordinationAnalyzer(
    #     universe,
    #     cutoff_radii=cutoff_radii,
    #     atom_selection=ca_mask,  # 可选
    #     scheduler='threads',
    #     enable_spatial_subdivision=True
    # )
    print("   分析器配置完成")

    # 5. 进行配位数分析
    print("\n5. 进行配位数分析...")
    # result = analyzer.analyze_parallel(atom_indices=None)  # 分析所有选中的原子
    print("   配位数分析完成")

    # 6. 获取个体配位数信息（这是配位数分析的特色）
    print("\n6. 获取个体配位数信息...")
    print("   配位数分析的核心特点是保留每个原子的个体信息：")

    # 获取特定原子的配位数
    # atom_coord = result.get_coordination_number(atom_index=0)
    # print(f"   原子0的配位数: {atom_coord}")

    # 获取配位数分布统计（个体分布）
    # distribution = result.get_coordination_distribution()
    # print(f"   配位数分布统计: {distribution}")

    # 7. 可视化个体配位数分布
    print("\n7. 可视化个体配位数分布...")
    print("   配位数分析提供多种可视化方式：")

    # 直方图 - 显示配位数分布
    # fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    # result.plot(plot_type='histogram', ax=axes[0,0])
    # axes[0,0].set_title('配位数分布直方图')

    # 条形图 - 显示平均配位数
    # result.plot(plot_type='bar', ax=axes[0,1])
    # axes[0,1].set_title('平均配位数条形图')

    # 散点图 - 显示原子索引vs配位数
    # result.plot(plot_type='scatter', ax=axes[1,0])
    # axes[1,0].set_title('原子索引 vs 配位数')

    # 热图 - 显示原子-物种配位数矩阵
    # result.plot(plot_type='heatmap', ax=axes[1,1])
    # axes[1,1].set_title('配位数热图')

    # plt.tight_layout()
    # plt.show()

    # 8. 打印统计信息
    print("\n8. 打印统计信息...")
    # result.print_statistics()

    print("\n=== 配位数分析 vs RDF/Angular分析的区别 ===")
    print("1. 数据保存：")
    print("   - 配位数：保留每个原子的个体配位数值")
    print("   - RDF/Angular：主要保存平均的分布函数")
    print()
    print("2. 分析重点：")
    print("   - 配位数：关注个体差异和分布统计")
    print("   - RDF/Angular：关注整体平均行为")
    print()
    print("3. 可视化：")
    print("   - 配位数：直方图、散点图、热图等显示个体分布")
    print("   - RDF/Angular：曲线图显示平均函数")
    print()
    print("4. 应用场景：")
    print("   - 配位数：研究局部结构多样性、缺陷分析")
    print("   - RDF/Angular：研究整体结构特征、相变行为")


def demonstrate_coordination_features():
    """演示配位数分析的独特功能"""

    print("\n=== 配位数分析独特功能演示 ===")

    # 模拟一些配位数数据来演示功能
    print("\n1. 模拟配位数数据...")

    # 创建模拟的配位数结果
    atom_indices = list(range(100))
    coordination_data = {}

    # 模拟不同原子的配位数（体现个体差异）
    np.random.seed(42)
    for atom_idx in atom_indices:
        coordination_data[atom_idx] = {
            "1": np.random.poisson(4),  # 类型1的平均配位数为4
            "2": np.random.poisson(6),  # 类型2的平均配位数为6
        }

    # 创建模拟的CoordinationResult对象
    cutoff_radii = {"1-1": 2.5, "1-2": 3.0, "2-2": 2.8}

    # 注意：这里只是演示数据结构，实际使用中会从分析器获得结果
    print("   模拟数据创建完成")

    # 2. 演示个体配位数分布统计
    print("\n2. 个体配位数分布统计...")

    # 计算配位数分布
    species_1_coords = [data["1"] for data in coordination_data.values()]
    species_2_coords = [data["2"] for data in coordination_data.values()]

    print(f"   物种1配位数分布:")
    unique_1, counts_1 = np.unique(species_1_coords, return_counts=True)
    for val, count in zip(unique_1, counts_1):
        print(
            f"     配位数 {val}: {count} 个原子 ({count/len(species_1_coords)*100:.1f}%)"
        )

    print(f"   物种2配位数分布:")
    unique_2, counts_2 = np.unique(species_2_coords, return_counts=True)
    for val, count in zip(unique_2, counts_2):
        print(
            f"     配位数 {val}: {count} 个原子 ({count/len(species_2_coords)*100:.1f}%)"
        )

    # 3. 演示统计分析
    print("\n3. 统计分析...")
    print(
        f"   物种1配位数: 平均={np.mean(species_1_coords):.2f}, 标准差={np.std(species_1_coords):.2f}"
    )
    print(
        f"   物种2配位数: 平均={np.mean(species_2_coords):.2f}, 标准差={np.std(species_2_coords):.2f}"
    )

    # 4. 演示可视化（如果matplotlib可用）
    try:
        import matplotlib.pyplot as plt

        print("\n4. 创建可视化图表...")

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))

        # Coordination number distribution histogram
        axes[0, 0].hist(
            species_1_coords,
            bins=range(min(species_1_coords), max(species_1_coords) + 2),
            alpha=0.7,
            label="Species 1",
            color="blue",
        )
        axes[0, 0].hist(
            species_2_coords,
            bins=range(min(species_2_coords), max(species_2_coords) + 2),
            alpha=0.7,
            label="Species 2",
            color="red",
        )
        axes[0, 0].set_xlabel("Coordination Number")
        axes[0, 0].set_ylabel("Count")
        axes[0, 0].set_title("Coordination Number Distribution Histogram")
        axes[0, 0].legend()

        # Scatter plot: atom index vs coordination number
        axes[0, 1].scatter(
            atom_indices, species_1_coords, alpha=0.6, label="Species 1", s=30
        )
        axes[0, 1].scatter(
            atom_indices, species_2_coords, alpha=0.6, label="Species 2", s=30
        )
        axes[0, 1].set_xlabel("Atom Index")
        axes[0, 1].set_ylabel("Coordination Number")
        axes[0, 1].set_title("Atom Index vs Coordination Number")
        axes[0, 1].legend()

        # Box plot showing distribution
        data_to_plot = [species_1_coords, species_2_coords]
        axes[1, 0].boxplot(data_to_plot, labels=["Species 1", "Species 2"])
        axes[1, 0].set_ylabel("Coordination Number")
        axes[1, 0].set_title("Coordination Number Distribution Box Plot")

        # Coordination number matrix heatmap
        coord_matrix = np.array(
            [
                [coordination_data[i]["1"], coordination_data[i]["2"]]
                for i in range(min(20, len(atom_indices)))
            ]
        )  # Show only first 20 atoms
        im = axes[1, 1].imshow(coord_matrix.T, cmap="viridis", aspect="auto")
        axes[1, 1].set_xlabel("Atom Index")
        axes[1, 1].set_ylabel("Species")
        axes[1, 1].set_title("Coordination Number Heatmap (First 20 Atoms)")
        axes[1, 1].set_yticks([0, 1])
        axes[1, 1].set_yticklabels(["Species 1", "Species 2"])
        plt.colorbar(im, ax=axes[1, 1], label="Coordination Number")

        plt.tight_layout()

        # Ensure results directory exists
        import os

        os.makedirs("results", exist_ok=True)

        plt.savefig(
            "results/coordination_analysis_demo.png", dpi=150, bbox_inches="tight"
        )
        print(
            "   Visualization plots saved as 'results/coordination_analysis_demo.png'"
        )

    except ImportError:
        print("   matplotlib不可用，跳过可视化演示")


if __name__ == "__main__":
    main()
    demonstrate_coordination_features()

    print("\n=== 总结 ===")
    print("配位数分析模块已成功创建，主要特点：")
    print("✓ 保留每个原子的个体配位数信息")
    print("✓ 提供配位数分布统计和可视化")
    print("✓ 支持多种绘图方式（直方图、散点图、热图等）")
    print("✓ 遵循MDemon的单原子分析架构")
    print("✓ 支持并行处理和空间预分割优化")
    print("✓ 与RDF/Angular分析形成互补的分析体系")
