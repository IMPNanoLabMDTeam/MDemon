"""Coordination Number analysis for single atoms in MDemon

This module provides coordination number analysis capabilities for single atoms, including:
- Coordination number calculation for individual atoms based on center-neighbor atom pairs
- Support for species-specific cutoff radii for different center-neighbor combinations
- Individual coordination number distribution analysis (not averaged like RDF/Angular)
- Clear distinction between center atom species and neighbor atom species
- Parallel processing support via Dask

Classes:
    CoordinationAnalyzer: Coordination number analysis for single atoms
    CoordinationResult: Result container for coordination number analysis

Note:
    The coordination analysis requires explicit specification of center-neighbor atom pairs
    through cutoff_radii dictionary with keys like '1-1', '1-2', etc., where the first
    number represents the center atom species and the second represents the neighbor species.
"""

import time
import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from dask import delayed

from ...lib.distance import distance_array
from .base import AnalysisResult, SingleAtomAnalyzer


class CoordinationResult(AnalysisResult):
    """
    配位数分析结果容器

    这个类扩展了基础的AnalysisResult，专门用于存储和处理个体配位数分析结果。
    与RDF和Angular不同，配位数分析重点关注每个原子的个体配位数分布，而不是平均值。

    Parameters
    ----------
    atom_indices : List[int]
        分析的原子索引列表
    coordination_data : dict
        配位数数据，格式为 {atom_index: {species: coordination_number}}
    cutoff_radii : dict
        使用的截断半径字典
    metadata : dict, optional
        额外的元数据信息

    Attributes
    ----------
    distributions : dict
        预计算的配位数分布数据，格式为 {f'{center_species}-{neighbor_species}': distribution_data}
        在初始化时自动计算，避免重复计算提高性能
    """

    def __init__(
        self,
        atom_indices,
        coordination_data,
        cutoff_radii,
        metadata=None,
        universe=None,
    ):
        # 将配位数特定数据添加到基础数据中
        data = coordination_data

        super().__init__(
            analysis_type="coordination",
            atom_indices=atom_indices,
            data=data,
            metadata=metadata,
            universe=universe,
        )
        self.coordination_data = coordination_data
        self.cutoff_radii = cutoff_radii

        # 预计算配位数分布数据
        self.distributions = self._compute_distributions()

    def get_coordination_number(self, atom_index, neighbor_species=None):
        """
        获取特定原子的配位数

        Parameters
        ----------
        atom_index : int
            中心原子索引
        neighbor_species : str, optional
            邻近原子类型，如果为None则返回该原子周围所有类型的配位数

        Returns
        -------
        dict or int
            如果neighbor_species为None：返回 {neighbor_species: coordination_number} 字典
            否则返回该邻近原子类型的配位数（整数）
        """
        if atom_index not in self.coordination_data:
            raise ValueError(f"Coordination data for atom {atom_index} not available")

        atom_coordination_data = self.coordination_data[atom_index]

        if neighbor_species is None:
            return atom_coordination_data.copy()
        else:
            if neighbor_species not in atom_coordination_data:
                raise ValueError(
                    f"Neighbor species '{neighbor_species}' not found for atom {atom_index}"
                )
            return atom_coordination_data[neighbor_species]

    def get_coordination_number_by_species(self, atom_index, species=None):
        """
        [DEPRECATED] 获取特定原子的配位数 - 旧版本接口

        此方法已弃用，请使用 get_coordination_number(atom_index, neighbor_species) 代替。

        Parameters
        ----------
        atom_index : int
            中心原子索引
        species : str, optional
            邻近原子类型

        Returns
        -------
        dict or int
            配位数数据
        """
        warnings.warn(
            "get_coordination_number_by_species() is deprecated. "
            "Use get_coordination_number(atom_index, neighbor_species) instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.get_coordination_number(atom_index, neighbor_species=species)

    def get_available_species(self, atom_index=None):
        """
        获取可用的原子类型列表

        Parameters
        ----------
        atom_index : int, optional
            特定原子的索引，如果为None则返回所有原子的所有类型

        Returns
        -------
        list
            原子类型列表
        """
        if atom_index is not None:
            if atom_index not in self.coordination_data:
                return []
            return list(self.coordination_data[atom_index].keys())
        else:
            # 返回所有原子的所有类型
            all_species = set()
            for atom_data in self.coordination_data.values():
                all_species.update(atom_data.keys())
            return list(all_species)

    def _compute_distributions(self):
        """
        预计算所有配位对的分布数据

        Returns
        -------
        dict
            配位数分布数据，格式为:
            {f'{center_species}-{neighbor_species}': {'values': [coordination_numbers], 'counts': [counts], 'statistics': {...}}}
        """
        # 获取所有可用的邻近原子类型
        available_neighbor_species = self.get_available_species()

        # 获取所有中心原子的类型
        if hasattr(self, "universe") and self.universe is not None:
            center_atom_species = [
                self.universe.atoms.species[idx]
                for idx in self.atom_indices
                if idx in self.coordination_data
            ]
            available_center_species = list(set(center_atom_species))
        else:
            # 如果没有universe信息，从cutoff_radii推断中心原子类型
            available_center_species = list(
                set([key.split("-")[0] for key in self.cutoff_radii.keys()])
            )

        distribution_data = {}

        for center_sp in available_center_species:
            for neighbor_sp in available_neighbor_species:
                pair_key = f"{center_sp}-{neighbor_sp}"

                # 收集该中心-邻近原子对的所有配位数
                coordination_numbers = []
                for atom_idx in self.atom_indices:
                    if (
                        atom_idx in self.coordination_data
                        and neighbor_sp in self.coordination_data[atom_idx]
                    ):
                        # 检查该原子是否为指定的中心原子类型
                        if hasattr(self, "universe") and self.universe is not None:
                            atom_species = self.universe.atoms.species[atom_idx]
                            if atom_species == center_sp:
                                coordination_numbers.append(
                                    self.coordination_data[atom_idx][neighbor_sp]
                                )
                        else:
                            # 如果没有universe信息，假设所有原子都符合条件
                            coordination_numbers.append(
                                self.coordination_data[atom_idx][neighbor_sp]
                            )

                if coordination_numbers:
                    coordination_numbers = np.array(coordination_numbers)

                    # 计算配位数分布
                    unique_values, counts = np.unique(
                        coordination_numbers, return_counts=True
                    )

                    # 计算统计量
                    statistics = {
                        "mean": np.mean(coordination_numbers),
                        "std": np.std(coordination_numbers),
                        "min": np.min(coordination_numbers),
                        "max": np.max(coordination_numbers),
                        "median": np.median(coordination_numbers),
                        "mode": unique_values[np.argmax(counts)],
                        "total_atoms": len(coordination_numbers),
                        "center_species": center_sp,
                        "neighbor_species": neighbor_sp,
                    }

                    distribution_data[pair_key] = {
                        "values": unique_values,
                        "counts": counts,
                        "frequencies": counts / len(coordination_numbers),
                        "raw_data": coordination_numbers,
                        "statistics": statistics,
                    }

        return distribution_data

    def get_coordination_distribution(self, center_species=None, neighbor_species=None):
        """
        获取配位数分布统计

        这是配位数分析的核心功能 - 统计个体配位数的分布情况
        配位数分析需要明确区分中心原子类型和邻近原子类型

        Parameters
        ----------
        center_species : str, optional
            中心原子类型，如果为None则包含所有中心原子类型
        neighbor_species : str, optional
            邻近原子类型，如果为None则包含所有邻近原子类型

        Returns
        -------
        dict
            配位数分布统计，格式为:
            {f'{center_species}-{neighbor_species}': {'values': [coordination_numbers], 'counts': [counts], 'statistics': {...}}}
        """
        # 从预计算的distributions中筛选数据
        filtered_data = {}

        for pair_key, data in self.distributions.items():
            center_sp, neighbor_sp = pair_key.split("-")

            # 检查是否符合筛选条件
            center_match = center_species is None or center_sp == center_species
            neighbor_match = neighbor_species is None or neighbor_sp == neighbor_species

            if center_match and neighbor_match:
                filtered_data[pair_key] = data

        return filtered_data

    def get_coordination_distribution_by_species(self, species=None):
        """
        [DEPRECATED] 获取配位数分布统计 - 旧版本接口

        此方法已弃用，请使用 get_coordination_distribution(center_species, neighbor_species) 代替。
        该方法将species参数解释为邻近原子类型。

        Parameters
        ----------
        species : str, optional
            邻近原子类型，如果为None则分别统计所有类型

        Returns
        -------
        dict
            配位数分布统计
        """
        warnings.warn(
            "get_coordination_distribution_by_species() is deprecated. "
            "Use get_coordination_distribution(center_species, neighbor_species) instead. "
            "The 'species' parameter is now interpreted as 'neighbor_species'.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.get_coordination_distribution(
            center_species=None, neighbor_species=species
        )

    def get_coordination_matrix(self, species_list=None):
        """
        获取配位数矩阵，行为原子，列为不同species的配位数

        Parameters
        ----------
        species_list : list, optional
            要包含的species列表，如果为None则使用所有available species

        Returns
        -------
        tuple
            (atom_indices_array, species_array, coordination_matrix)
            coordination_matrix[i, j] 表示第i个原子对第j个species的配位数
        """
        if species_list is None:
            species_list = sorted(self.get_available_species())

        # 过滤出有数据的原子
        valid_atom_indices = []
        for atom_idx in self.atom_indices:
            if atom_idx in self.coordination_data:
                valid_atom_indices.append(atom_idx)

        if not valid_atom_indices:
            return np.array([]), np.array([]), np.array([])

        # 创建配位数矩阵
        coordination_matrix = np.zeros(
            (len(valid_atom_indices), len(species_list)), dtype=int
        )

        for i, atom_idx in enumerate(valid_atom_indices):
            atom_data = self.coordination_data[atom_idx]
            for j, species in enumerate(species_list):
                if species in atom_data:
                    coordination_matrix[i, j] = atom_data[species]

        return np.array(valid_atom_indices), np.array(species_list), coordination_matrix

    def plot(
        self,
        center_species=None,
        neighbor_species=None,
        plot_type="histogram",
        ax=None,
        **kwargs,
    ):
        """
        绘制配位数分布图。默认绘制直方图显示个体配位数分布。

        Parameters
        ----------
        center_species : str, optional
            中心原子类型，如果为None则绘制所有中心原子类型
        neighbor_species : str, optional
            邻近原子类型，如果为None则绘制所有邻近原子类型
        plot_type : str, optional
            绘图类型：'histogram', 'bar', 'scatter', 'heatmap'
            默认为 'histogram'
        ax : matplotlib.axes.Axes, optional
            用于绘图的matplotlib轴对象
        **kwargs : dict
            传递给绘图函数的其他关键字参数

        Returns
        -------
        matplotlib.axes.Axes
            包含配位数分布图的matplotlib轴对象
        """
        return self.plot_coordination_distribution(
            center_species=center_species,
            neighbor_species=neighbor_species,
            plot_type=plot_type,
            ax=ax,
            **kwargs,
        )

    def plot_coordination_distribution(
        self,
        center_species=None,
        neighbor_species=None,
        plot_type="histogram",
        ax=None,
        show_statistics=True,
        **kwargs,
    ):
        """
        绘制配位数分布图

        Parameters
        ----------
        center_species : str, optional
            中心原子类型
        neighbor_species : str, optional
            邻近原子类型
        plot_type : str, optional
            绘图类型：'histogram', 'bar', 'scatter', 'heatmap'
        ax : matplotlib.axes.Axes, optional
            绘图轴
        show_statistics : bool, optional
            是否在图上显示统计信息
        **kwargs : dict
            传递给绘图函数的额外参数
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib is required for plotting")

        distribution_data = self.get_coordination_distribution(
            center_species, neighbor_species
        )

        if not distribution_data:
            warnings.warn("No coordination data available for plotting.")
            return

        if plot_type == "histogram":
            return self._plot_histogram(
                distribution_data, ax, show_statistics, **kwargs
            )
        elif plot_type == "bar":
            return self._plot_bar(distribution_data, ax, show_statistics, **kwargs)
        elif plot_type == "scatter":
            return self._plot_scatter(distribution_data, ax, **kwargs)
        elif plot_type == "heatmap":
            return self._plot_heatmap(ax, **kwargs)
        else:
            raise ValueError(f"Unsupported plot_type: {plot_type}")

    def _plot_histogram(
        self, distribution_data, ax=None, show_statistics=True, **kwargs
    ):
        """绘制配位数直方图"""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib is required for plotting")

        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))

        colors = plt.cm.tab10(np.linspace(0, 1, len(distribution_data)))

        for i, (pair_key, data) in enumerate(distribution_data.items()):
            values = data["values"]
            counts = data["counts"]
            stats = data["statistics"]

            # 绘制直方图
            ax.bar(
                values + i * 0.1,
                counts,
                width=0.8,
                alpha=0.7,
                color=colors[i],
                label=f"{pair_key}",
                **kwargs,
            )

            if show_statistics:
                # 添加统计信息到图例
                mean_val = stats["mean"]
                std_val = stats["std"]
                ax.axvline(
                    mean_val,
                    color=colors[i],
                    linestyle="--",
                    alpha=0.8,
                    label=f"{pair_key}: μ={mean_val:.1f}, σ={std_val:.1f}",
                )

        ax.set_xlabel("Coordination Number")
        ax.set_ylabel("Count")
        ax.set_title("Coordination Number Distribution")
        ax.legend()
        ax.grid(True, alpha=0.3)

        return ax

    def _plot_bar(self, distribution_data, ax=None, show_statistics=True, **kwargs):
        """绘制配位数条形图"""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib is required for plotting")

        if ax is None:
            fig, ax = plt.subplots(figsize=(12, 6))

        pair_keys = list(distribution_data.keys())
        x_pos = np.arange(len(pair_keys))

        means = [
            distribution_data[pair_key]["statistics"]["mean"] for pair_key in pair_keys
        ]
        stds = [
            distribution_data[pair_key]["statistics"]["std"] for pair_key in pair_keys
        ]

        bars = ax.bar(x_pos, means, yerr=stds, capsize=5, alpha=0.7, **kwargs)

        ax.set_xlabel("Center-Neighbor Pair")
        ax.set_ylabel("Average Coordination Number")
        ax.set_title("Average Coordination Number by Center-Neighbor Pair")
        ax.set_xticks(x_pos)
        ax.set_xticklabels(pair_keys, rotation=45)
        ax.grid(True, alpha=0.3, axis="y")

        if show_statistics:
            # 在柱子上显示数值
            for i, (bar, mean, std) in enumerate(zip(bars, means, stds)):
                height = bar.get_height()
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    height + std + 0.1,
                    f"{mean:.1f}±{std:.1f}",
                    ha="center",
                    va="bottom",
                )

        return ax

    def _plot_scatter(self, distribution_data, ax=None, **kwargs):
        """绘制配位数散点图（原子索引 vs 配位数）"""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib is required for plotting")

        if ax is None:
            fig, ax = plt.subplots(figsize=(12, 8))

        colors = plt.cm.tab10(np.linspace(0, 1, len(distribution_data)))

        for i, (pair_key, data) in enumerate(distribution_data.items()):
            stats = data["statistics"]
            center_sp = stats["center_species"]
            neighbor_sp = stats["neighbor_species"]

            # 获取原子索引和对应的配位数
            atom_indices = []
            coordination_numbers = []

            for atom_idx in self.atom_indices:
                if (
                    atom_idx in self.coordination_data
                    and neighbor_sp in self.coordination_data[atom_idx]
                ):
                    # 检查该原子是否为指定的中心原子类型
                    if hasattr(self, "universe") and self.universe is not None:
                        atom_species = self.universe.atoms.species[atom_idx]
                        if atom_species == center_sp:
                            atom_indices.append(atom_idx)
                            coordination_numbers.append(
                                self.coordination_data[atom_idx][neighbor_sp]
                            )
                    else:
                        atom_indices.append(atom_idx)
                        coordination_numbers.append(
                            self.coordination_data[atom_idx][neighbor_sp]
                        )

            ax.scatter(
                atom_indices,
                coordination_numbers,
                color=colors[i],
                alpha=0.7,
                s=50,
                label=f"{pair_key}",
                **kwargs,
            )

        ax.set_xlabel("Atom Index")
        ax.set_ylabel("Coordination Number")
        ax.set_title("Coordination Number vs Atom Index")
        ax.legend()
        ax.grid(True, alpha=0.3)

        return ax

    def _plot_heatmap(self, ax=None, **kwargs):
        """绘制配位数热图（原子 vs species）"""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib is required for plotting")

        atom_indices, species_list, coord_matrix = self.get_coordination_matrix()

        if len(atom_indices) == 0:
            warnings.warn("No data available for heatmap.")
            return

        if ax is None:
            fig, ax = plt.subplots(
                figsize=(max(8, len(species_list)), max(6, len(atom_indices) // 10))
            )

        im = ax.imshow(coord_matrix, cmap="viridis", aspect="auto", **kwargs)

        # 设置标签
        ax.set_xticks(range(len(species_list)))
        ax.set_xticklabels(species_list)
        ax.set_ylabel("Atom Index")
        ax.set_title("Coordination Number Heatmap")

        # 添加颜色条
        plt.colorbar(im, ax=ax, label="Coordination Number")

        # 如果原子数量不太多，显示y轴标签
        if len(atom_indices) <= 50:
            ax.set_yticks(range(len(atom_indices)))
            ax.set_yticklabels(atom_indices)

        return ax

    def print_statistics(self, center_species=None, neighbor_species=None):
        """打印配位数统计信息"""
        distribution_data = self.get_coordination_distribution(
            center_species, neighbor_species
        )

        print("\n=== Coordination Number Statistics ===")
        for pair_key, data in distribution_data.items():
            stats = data["statistics"]
            center_sp = stats["center_species"]
            neighbor_sp = stats["neighbor_species"]
            print(f"\nCenter-Neighbor Pair {center_sp}-{neighbor_sp}:")
            print(f"  Total atoms analyzed: {stats['total_atoms']}")
            print(f"  Mean coordination: {stats['mean']:.2f}")
            print(f"  Standard deviation: {stats['std']:.2f}")
            print(f"  Range: {stats['min']} - {stats['max']}")
            print(f"  Median: {stats['median']:.1f}")
            print(f"  Mode: {stats['mode']}")

            # 显示分布
            values = data["values"]
            frequencies = data["frequencies"]
            print(f"  Distribution:")
            for val, freq in zip(values, frequencies):
                print(
                    f"    CN={val}: {freq*100:.1f}% ({int(freq*stats['total_atoms'])} atoms)"
                )
        print("=" * 40)


class CoordinationAnalyzer(SingleAtomAnalyzer):
    """
    配位数分析器

    这个类提供了单原子配位数分析的核心功能，包括：
    - 单个原子的配位数计算
    - 并行处理多个原子的配位数
    - 支持不同原子类型的特定截断半径
    - 重点关注个体配位数分布而非平均值

    Parameters
    ----------
    universe : MDemon.Universe
        分析的宇宙对象
    cutoff_radii : dict
        中心-邻近原子对的截断半径字典，格式为 {'1-1': 2.5, '1-2': 3.0, ...}
        其中键表示中心-邻近原子对：
        - '1-1' 表示中心原子类型1，邻近原子类型1
        - '1-2' 表示中心原子类型1，邻近原子类型2
        此参数是必需的，因为截断半径是高度特异的
    atom_selection : numpy.ndarray, optional
        原子选择掩码，长度为len(universe.atoms)的布尔数组
    scheduler : str, optional
        Dask调度器类型，Default: 'threads'
    enable_spatial_subdivision : bool, optional
        是否启用空间预分割优化，Default: True
    **kwargs : dict
        传递给父类的其他参数
    """

    def __init__(
        self,
        universe,
        cutoff_radii,
        atom_selection=None,
        scheduler="threads",
        enable_spatial_subdivision=True,
        **kwargs,
    ):
        # 验证截断半径参数
        if cutoff_radii is None:
            raise ValueError(
                "cutoff_radii is required and must be a dictionary with format "
                "{'1-1': 2.5, '1-2': 3.0, ...} where keys represent atom type pairs "
                "(e.g., '1-1' for same type, '1-2' for different types)"
            )

        if not isinstance(cutoff_radii, dict):
            raise TypeError("cutoff_radii must be a dictionary")

        self.cutoff_radii = cutoff_radii

        # 使用固定的10埃作为空间预分割的截断半径
        # 完全参考angular.py的策略，避免使用过小的cutoff影响空间预分割效果
        r_cutoff = 10.0 if enable_spatial_subdivision else None

        super().__init__(
            universe,
            atom_selection=atom_selection,
            scheduler=scheduler,
            enable_spatial_subdivision=enable_spatial_subdivision,
            r_cutoff=r_cutoff,
            **kwargs,
        )

        # 在初始化时获取系统中的所有唯一物种类型
        self.unique_species = np.unique(self.universe.atoms.species)

    def _get_cutoff_for_pair(self, center_species, neighbor_species):
        """
        获取特定中心-邻近原子对的截断半径

        Parameters
        ----------
        center_species : str
            中心原子类型
        neighbor_species : str
            邻近原子类型

        Returns
        -------
        float
            截断半径

        Notes
        -----
        键格式为 "center-neighbor"，其中：
        - "1-1" 表示中心原子类型1，邻近原子类型1
        - "1-2" 表示中心原子类型1，邻近原子类型2
        """
        pair_key = f"{center_species}-{neighbor_species}"

        if pair_key in self.cutoff_radii:
            return self.cutoff_radii[pair_key]
        else:
            # 抛出错误，不再使用默认值
            raise ValueError(
                f"No cutoff radius found for center-neighbor pair {pair_key}. "
                f"Available keys: {list(self.cutoff_radii.keys())}"
            )

    def analyze_single_atom(self, atom_index, **kwargs):
        """
        分析单个原子的配位数，返回按原子类型分组的结果

        Parameters
        ----------
        atom_index : int
            中心原子的索引
        **kwargs : dict
            其他参数

        Returns
        -------
        dict
            按原子类型分组的配位数字典，格式为 {species: coordination_number}
        """
        if self.show_running_time:
            start_total = time.perf_counter()
            print(f"\n[COORDINATION DEBUG] === Analyzing atom {atom_index} ===")

        # 获取原子坐标和原子类型
        if self.show_running_time:
            t1 = time.perf_counter()

        atoms = self.universe.atoms
        all_coordinates = atoms.coordinate
        all_species = atoms.species
        center_species = all_species[atom_index]

        if self.show_running_time:
            t2 = time.perf_counter()
            print(f"[COORDINATION DEBUG] 获取原子坐标和类型: {(t2-t1)*1000:.3f} ms")
            print(f"[COORDINATION DEBUG] 中心原子类型: {center_species}")

        # 使用空间预分割获取相关原子（如果启用）
        if self.show_running_time:
            t3 = time.perf_counter()

        relevant_atom_indices, combined_mask = self._get_relevant_atoms_for_analysis(
            atom_index
        )

        if self.show_running_time:
            t4 = time.perf_counter()
            print(f"[COORDINATION DEBUG] 空间预分割获取相关原子: {(t4-t3)*1000:.3f} ms")

        # 筛选坐标和物种信息
        if self.show_running_time:
            t5 = time.perf_counter()

        if combined_mask is not None:
            coordinates = all_coordinates[combined_mask]
            species = all_species[combined_mask]
            # 找到中心原子在筛选后数组中的新位置
            center_index_in_filtered = (
                np.where(combined_mask)[0].tolist().index(atom_index)
            )
        else:
            coordinates = all_coordinates[relevant_atom_indices]
            species = all_species[relevant_atom_indices]
            center_index_in_filtered = relevant_atom_indices.tolist().index(atom_index)

        if self.show_running_time:
            t6 = time.perf_counter()
            print(f"[COORDINATION DEBUG] 筛选坐标和物种信息: {(t6-t5)*1000:.3f} ms")
            print(f"[COORDINATION DEBUG] 筛选后原子数: {len(coordinates)}")

        # 计算配位数
        if self.show_running_time:
            t7 = time.perf_counter()

        coordination_results = self._compute_coordination_numbers(
            center_index_in_filtered, coordinates, species
        )

        if self.show_running_time:
            t8 = time.perf_counter()
            print(f"[COORDINATION DEBUG] 配位数计算: {(t8-t7)*1000:.3f} ms")
            print(
                f"[COORDINATION DEBUG] 计算得到的物种类型数: {len(coordination_results)}"
            )

        if self.show_running_time:
            t9 = time.perf_counter()
            print(
                f"[COORDINATION DEBUG] === 总时间: {(t9-start_total)*1000:.3f} ms ===\n"
            )

        return coordination_results

    def _compute_coordination_numbers(self, center_index, coordinates, species):
        """
        计算配位数

        Parameters
        ----------
        center_index : int
            中心原子在坐标数组中的索引
        coordinates : numpy.ndarray
            原子坐标数组，形状为(n_atoms, 3)
        species : numpy.ndarray
            原子类型数组，形状为(n_atoms,)

        Returns
        -------
        dict
            包含各种原子类型配位数的字典，格式为 {species: coordination_number}
        """
        if self.show_running_time:
            start_coord = time.perf_counter()
            print(f"  [COORDINATION DEBUG] 开始配位数计算")

        center_pos = coordinates[center_index]  # 中心原子位置，1D数组
        center_species = species[center_index]
        box_dims = self.universe.box[:6]

        # 计算所有原子到中心原子的距离
        if self.show_running_time:
            t1 = time.perf_counter()

        distances = distance_array(
            reference=center_pos.reshape(1, -1), configuration=coordinates, box=box_dims
        ).flatten()

        if self.show_running_time:
            t2 = time.perf_counter()
            print(f"  [COORDINATION DEBUG] 距离计算: {(t2-t1)*1000:.3f} ms")

        # 移除中心原子本身
        tolerance = 1e-10
        valid_mask = distances > tolerance
        distances = distances[valid_mask]
        nearby_species = species[valid_mask]

        if self.show_running_time:
            t3 = time.perf_counter()
            print(f"  [COORDINATION DEBUG] 过滤中心原子: {(t3-t2)*1000:.3f} ms")
            print(f"  [COORDINATION DEBUG] 邻近原子数: {len(distances)}")

        coordination_results = {}

        # 对于每个物种类型，计算对应的配位数
        if self.show_running_time:
            t4 = time.perf_counter()

        for sp in self.unique_species:
            try:
                # 获取对应的截断半径
                cutoff = self._get_cutoff_for_pair(center_species, sp)

                # 筛选出在截断半径内且类型为sp的原子
                sp_mask = (nearby_species == sp) & (distances <= cutoff)
                coordination_number = np.sum(sp_mask)

                coordination_results[sp] = coordination_number

                if self.show_running_time:
                    print(
                        f"  [COORDINATION DEBUG] {center_species}-{sp}: {coordination_number} 个邻近原子"
                    )

            except ValueError:
                # 如果没有找到对应的截断半径配置，跳过此物种
                continue

        if self.show_running_time:
            t5 = time.perf_counter()
            print(f"  [COORDINATION DEBUG] 配位数计算: {(t5-t4)*1000:.3f} ms")
            print(
                f"  [COORDINATION DEBUG] 配位数计算总时间: {(t5-start_coord)*1000:.3f} ms"
            )

        return coordination_results

    def _create_analysis_result(self, result_dict, **kwargs):
        """
        创建配位数特定的分析结果对象

        重写基类方法以返回CoordinationResult而不是AnalysisResult

        Parameters
        ----------
        result_dict : dict
            分析结果字典，键为原子索引，值为配位数字典
        **kwargs : dict
            额外的元数据

        Returns
        -------
        CoordinationResult
            配位数分析结果对象
        """
        atom_indices = list(result_dict.keys())

        # 结果数据已经是正确的格式：{atom_index: {species: coordination_number}}
        coordination_data = result_dict

        # 创建元数据
        metadata = {
            "cutoff_radii": self.cutoff_radii,
            "timestep": self.universe.timestep,
            "scheduler": self.parallel_manager.scheduler,
            "n_atoms": len(atom_indices),
        }
        metadata.update(kwargs)

        return CoordinationResult(
            universe=self.universe,
            atom_indices=atom_indices,
            coordination_data=coordination_data,
            cutoff_radii=self.cutoff_radii,
            metadata=metadata,
        )

    def compute_coordination_statistics(
        self, atom_indices=None, species=None, **kwargs
    ):
        """
        计算配位数统计信息

        Parameters
        ----------
        atom_indices : List[int], optional
            要分析的原子索引列表
        species : str, optional
            原子类型，如果为None则统计所有类型
        **kwargs : dict
            其他参数

        Returns
        -------
        dict
            配位数统计信息
        """
        # 并行分析所有原子
        result = self.analyze_parallel(atom_indices=atom_indices, **kwargs)

        # 获取配位数分布
        distribution_data = result.get_coordination_distribution(species)

        return distribution_data
