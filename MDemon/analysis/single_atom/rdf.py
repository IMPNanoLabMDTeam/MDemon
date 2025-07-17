"""
Radial Distribution Function (RDF) analysis for single atoms in MDemon

This module provides RDF analysis capabilities for single atoms, including:
- Basic RDF calculation for single atoms
- Parallel processing support via Dask

Classes:
    RDFAnalyzer: Basic RDF analysis for single atoms
    RDFResult: Result container for RDF analysis
"""

import time
import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from dask import delayed

from ...lib.distance import distance_array
from .base import AnalysisResult, SingleAtomAnalyzer


class RDFResult(AnalysisResult):
    """
    RDF分析结果容器

    这个类扩展了基础的AnalysisResult，专门用于存储和处理按原子类型分组的RDF分析结果。

    Parameters
    ----------
    atom_indices : List[int]
        分析的原子索引列表
    r_bins : numpy.ndarray
        径向距离分箱边界
    species_rdf_data : dict
        按原子类型分组的RDF数据，格式为 {atom_index: {species: rdf_values}}
    metadata : dict, optional
        额外的元数据信息
    """

    def __init__(
        self, atom_indices, r_bins, species_rdf_data, metadata=None, universe=None
    ):
        # 将RDF特定数据添加到基础数据中
        data = species_rdf_data

        super().__init__(
            analysis_type="rdf",
            atom_indices=atom_indices,
            data=data,
            metadata=metadata,
            universe=universe,
        )
        self.r_bins = r_bins
        self.species_rdf_data = species_rdf_data

    def get_species_rdf(self, atom_index, species=None):
        """
        获取特定原子与特定原子类型的RDF

        Parameters
        ----------
        atom_index : int
            原子索引
        species : str, optional
            原子类型，如果为None则返回所有类型

        Returns
        -------
        dict or tuple
            如果species为None：返回 {species: (r_values, rdf_values)} 字典
            否则返回 (r_values, rdf_values) 元组
        """
        if atom_index not in self.species_rdf_data:
            raise ValueError(f"RDF data for atom {atom_index} not available")

        atom_species_data = self.species_rdf_data[atom_index]
        r_values = self.r_bins[:-1]

        if species is None:
            # 返回所有类型的RDF
            result = {}
            for sp, rdf_vals in atom_species_data.items():
                result[sp] = (r_values, rdf_vals)
            return result
        else:
            # 返回特定类型的RDF
            if species not in atom_species_data:
                raise ValueError(f"Species '{species}' not found for atom {atom_index}")
            return r_values, atom_species_data[species]

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
            if atom_index not in self.species_rdf_data:
                return []
            return list(self.species_rdf_data[atom_index].keys())
        else:
            # 返回所有原子的所有类型
            all_species = set()
            for atom_data in self.species_rdf_data.values():
                all_species.update(atom_data.keys())
            return list(all_species)

    def get_first_peak(self, atom_index, species):
        """
        获取特定原子和原子类型的第一个峰位置

        Parameters
        ----------
        atom_index : int
            原子索引
        species : str
            原子类型

        Returns
        -------
        float
            第一峰位置
        """
        r_values, rdf_vals = self.get_species_rdf(atom_index, species)
        peak_idx = np.argmax(rdf_vals)
        return r_values[peak_idx]

    def plot(self, species_pairs=None, ax=None, show_std=True, **kwargs):
        """
        绘制平均RDF图。默认绘制所有物种对的组合图。

        此方法计算并绘制所选物种对之间的平均径向分布函数(RDF)。
        它提供了显示标准差区间的选项，并允许在给定的matplotlib轴上进行绘制。

        Parameters
        ----------
        species_pairs : list of tuple, optional
            要绘制的物种对列表，例如 `[('C', 'H'), ('O', 'O')]`。
            如果为 `None` (默认)，将绘制所有可用物种对的组合图。
        ax : matplotlib.axes.Axes, optional
            用于绘图的matplotlib轴对象。如果为 `None`，将创建一个新的图和轴。
        show_std : bool, optional
            是否显示标准差的阴影区域。默认为 `True`。
        **kwargs : dict
            传递给 `matplotlib.pyplot.plot` 的其他关键字参数。

        Returns
        -------
        matplotlib.axes.Axes
            包含RDF图的matplotlib轴对象。

        Raises
        ------
        ValueError
            如果没有可用于绘图的数据。
        ImportError
            如果系统中未安装matplotlib。
        """
        return self.plot_average_rdf(
            species_pairs=species_pairs,
            ax=ax,
            show_std=show_std,
            **kwargs,
        )

    def plot_average_rdf(
        self, species_pairs=None, plot_type="combined", show_std=True, **kwargs
    ):
        """
        计算并绘制物种对之间的平均RDF。

        此高级绘图功能可以为指定的物种对生成单独的或组合的RDF图。
        支持显示平均RDF的标准差。

        Parameters
        ----------
        species_pairs : list of tuple, optional
            要绘制的物种对列表, e.g., `[('C', 'H'), ('O', 'O')]`。
            如果为 `None`, 将使用所有可能的物种对。
        plot_type : {'combined', 'separate'}, optional
            绘图类型。'combined' 将所有RDF绘制在一张图上，'separate'
            为每对物种创建一个单独的图。默认为 'combined'。
        show_std : bool, optional
            是否显示标准差的阴影区域。默认为 `True`。
        **kwargs : dict
            传递给 `matplotlib.pyplot.plot` 的其他关键字参数。

        Returns
        -------
        list of matplotlib.axes.Axes or matplotlib.axes.Axes
            如果 `plot_type` 是 'separate'，返回轴对象列表。
            如果 `plot_type` 是 'combined'，返回单个轴对象。
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib is required for plotting")

        averaged_data = self._calculate_average_rdf(species_pairs)

        if not averaged_data:
            warnings.warn("No data available for plotting.")
            return

        if plot_type == "combined":
            fig, ax = plt.subplots(figsize=(12, 8))
            colors = plt.cm.tab10(np.linspace(0, 1, len(averaged_data)))

            # 过滤掉可能冲突的参数
            plot_kwargs = {k: v for k, v in kwargs.items() if k not in ["ax"]}

            for i, ((sp1, sp2), data) in enumerate(averaged_data.items()):
                r_vals = data["r_values"]
                avg_rdf = data["avg_rdf"]
                n_samples = data["n_samples"]

                ax.plot(
                    r_vals,
                    avg_rdf,
                    color=colors[i],
                    linewidth=2,
                    label=f"Species {sp1}-{sp2} (n={n_samples})",
                    **plot_kwargs,
                )

            ax.set_xlabel("Distance (Å)")
            ax.set_ylabel("g(r)")
            ax.set_title("RDF Comparison - All Species Pairs")
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            return ax

        elif plot_type == "separate":
            axes = []
            # 过滤掉可能冲突的参数
            plot_kwargs = {k: v for k, v in kwargs.items() if k not in ["ax"]}

            for (sp1, sp2), data in averaged_data.items():
                fig, ax = plt.subplots(figsize=(10, 6))

                r_vals = data["r_values"]
                avg_rdf = data["avg_rdf"]
                std_rdf = data["std_rdf"]
                n_samples = data["n_samples"]

                ax.plot(
                    r_vals,
                    avg_rdf,
                    linewidth=2,
                    label=f"Species {sp1}-{sp2} (n={n_samples})",
                    **plot_kwargs,
                )
                if show_std:
                    ax.fill_between(
                        r_vals, avg_rdf - std_rdf, avg_rdf + std_rdf, alpha=0.3
                    )

                ax.set_xlabel("Distance (Å)")
                ax.set_ylabel("g(r)")
                ax.set_title(f"RDF: Species {sp1} - Species {sp2}")
                ax.legend()
                ax.grid(True, alpha=0.3)
                axes.append(ax)
            return axes

    def _calculate_average_rdf(self, species_pairs=None):
        """计算指定物种对的平均RDF、标准差和样本数。"""
        if self.universe is None:
            raise ValueError("Universe object is required for species-based averaging.")

        available_species = self.get_available_species()

        if species_pairs is None:
            pairs = []
            for i, sp1 in enumerate(available_species):
                for j, sp2 in enumerate(available_species):
                    if i <= j:
                        pairs.append((sp1, sp2))
            species_pairs = pairs

        averaged_data = {}
        analyzed_atom_indices = set(self.atom_indices)

        for sp1, sp2 in species_pairs:
            rdf_values_list = []

            sp1_mask = self.universe.atoms.species == sp1
            sp1_indices_in_universe = set(np.where(sp1_mask)[0])
            sp1_analyzed_indices = list(
                analyzed_atom_indices.intersection(sp1_indices_in_universe)
            )

            if not sp1_analyzed_indices:
                continue

            for atom_idx in sp1_analyzed_indices:
                if (
                    atom_idx in self.species_rdf_data
                    and sp2 in self.species_rdf_data[atom_idx]
                ):
                    rdf_values_list.append(self.species_rdf_data[atom_idx][sp2])

            if rdf_values_list:
                rdf_array = np.array(rdf_values_list)
                avg_rdf = np.mean(rdf_array, axis=0)
                std_rdf = np.std(rdf_array, axis=0)

                averaged_data[(sp1, sp2)] = {
                    "r_values": self.r_bins[:-1],
                    "avg_rdf": avg_rdf,
                    "std_rdf": std_rdf,
                    "n_samples": len(rdf_values_list),
                }
        self.averaged_data = averaged_data
        return averaged_data

    def plot_species_comparison(self, atom_index, ax=None, **kwargs):
        """
        绘制单个原子与不同原子类型的RDF比较图

        Parameters
        ----------
        atom_index : int
            原子索引
        ax : matplotlib.axes.Axes, optional
            绘图轴
        **kwargs : dict
            传递给plot的额外参数
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib is required for plotting")

        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))

        if atom_index not in self.species_rdf_data:
            raise ValueError(f"RDF data for atom {atom_index} not available")

        # 绘制各个原子类型的RDF
        species_data = self.get_species_rdf(atom_index)
        for species, (r_vals, rdf_vals) in species_data.items():
            ax.plot(r_vals, rdf_vals, label=f"vs {species}", alpha=0.8)

        ax.set_xlabel("Distance (Å)")
        ax.set_ylabel("g(r)")
        ax.set_title(f"RDF for Atom {atom_index} by Species")
        ax.legend()
        ax.grid(True, alpha=0.3)

        return ax


class RDFAnalyzer(SingleAtomAnalyzer):
    """
    径向分布函数分析器

    这个类提供了单原子RDF分析的核心功能，包括：
    - 单个原子的RDF计算
    - 并行处理多个原子的RDF
    - 与不同类型的参考原子计算RDF

    Parameters
    ----------
    universe : MDemon.Universe
        分析的宇宙对象
    atom_selection : numpy.ndarray, optional
        原子选择掩码，长度为len(universe.atoms)的布尔数组
    r_range : tuple, optional
        径向距离范围 (r_min, r_max)，Default: (0.0, 10.0)
    n_bins : int, optional
        距离分箱数量，Default: 100
    scheduler : str, optional
        Dask调度器类型，Default: 'threads'
    enable_spatial_subdivision : bool, optional
        是否启用空间预分割优化，Default: True
    use_subdivision_masks : bool, optional
        是否在空间预分割中使用掩码，Default: False。对于大体系建议保持False
    **kwargs : dict
        传递给父类的其他参数
    """

    def __init__(
        self,
        universe,
        atom_selection=None,
        r_range=(0.0, 10.0),
        n_bins=100,
        scheduler="threads",
        enable_spatial_subdivision=True,
        use_subdivision_masks=False,
        **kwargs,
    ):
        # 使用r_range的最大值作为空间预分割的截断半径
        r_cutoff = r_range[1] if enable_spatial_subdivision else None

        super().__init__(
            universe,
            atom_selection=atom_selection,
            scheduler=scheduler,
            enable_spatial_subdivision=enable_spatial_subdivision,
            r_cutoff=r_cutoff if r_cutoff > 10 else 10,
            use_subdivision_masks=use_subdivision_masks,
            **kwargs,
        )
        self.r_range = r_range
        self.n_bins = n_bins
        self.r_bins = np.linspace(r_range[0], r_range[1], n_bins + 1)
        self.dr = self.r_bins[1] - self.r_bins[0]

        # 预计算并缓存球壳体积，避免重复计算
        self._shell_volumes = self._compute_shell_volumes()

    def _compute_shell_volumes(self):
        """
        计算球壳体积并缓存

        Returns
        -------
        numpy.ndarray
            每个球壳的体积数组
        """
        r_inner = self.r_bins[:-1]
        r_outer = self.r_bins[1:]
        shell_volumes = (4.0 / 3.0) * np.pi * (r_outer**3 - r_inner**3)

        # 避免除零
        shell_volumes = np.where(shell_volumes > 0, shell_volumes, 1.0)

        return shell_volumes

    def analyze_single_atom(self, atom_index, **kwargs):
        """
        分析单个原子的RDF，返回按原子类型分组的结果

        Parameters
        ----------
        atom_index : int
            中心原子的索引
        **kwargs : dict
            其他参数

        Returns
        -------
        dict
            按原子类型分组的RDF字典，格式为 {species: (r_values, rdf_values)}
        """
        if self.show_running_time:
            start_total = time.perf_counter()
            print(f"\n[DEBUG] === Analyzing atom {atom_index} ===")

        # 获取原子坐标和原子类型 - 使用MDemon规范的属性名称
        if self.show_running_time:
            t1 = time.perf_counter()

        atoms = self.universe.atoms
        all_coordinates = atoms.coordinate
        all_species = atoms.species

        if self.show_running_time:
            t2 = time.perf_counter()
            print(f"[DEBUG] 获取原子坐标和类型: {(t2-t1)*1000:.3f} ms")
            print(f"[DEBUG] 总原子数: {len(all_coordinates)}")

        # 使用空间预分割获取相关原子（如果启用）
        if self.show_running_time:
            t3 = time.perf_counter()

        relevant_atom_indices, combined_mask = self._get_relevant_atoms_for_analysis(
            atom_index
        )

        if self.show_running_time:
            t4 = time.perf_counter()
            print(f"[DEBUG] 空间预分割获取相关原子: {(t4-t3)*1000:.3f} ms")
            print(
                f"[DEBUG] 相关原子数量: {len(relevant_atom_indices) if relevant_atom_indices is not None else 'All'}"
            )

        # 筛选坐标和物种信息
        if self.show_running_time:
            t5 = time.perf_counter()

        if combined_mask is not None:
            coordinates = all_coordinates[combined_mask]
            species = all_species[combined_mask]
        else:
            coordinates = all_coordinates[relevant_atom_indices]
            species = all_species[relevant_atom_indices]

        if self.show_running_time:
            t6 = time.perf_counter()
            print(f"[DEBUG] 筛选坐标和物种信息: {(t6-t5)*1000:.3f} ms")
            print(f"[DEBUG] 筛选后原子数: {len(coordinates)}")

        # 获取中心原子在筛选后数组中的新位置
        # 中心原子应该在相关原子中，找到它的新索引
        if self.show_running_time:
            t7 = time.perf_counter()

        original_center_pos = all_coordinates[atom_index : atom_index + 1]  # 保持2D形状

        if self.show_running_time:
            t8 = time.perf_counter()
            print(f"[DEBUG] 获取中心原子位置: {(t8-t7)*1000:.3f} ms")

        # 计算RDF（传入筛选后的坐标和物种）
        if self.show_running_time:
            t9 = time.perf_counter()

        rdf_results = self._compute_rdf_optimized(
            original_center_pos, coordinates, species
        )

        if self.show_running_time:
            t10 = time.perf_counter()
            print(f"[DEBUG] RDF计算: {(t10-t9)*1000:.3f} ms")
            print(f"[DEBUG] 计算得到的物种数: {len(rdf_results)}")

        # 重新格式化结果
        if self.show_running_time:
            t11 = time.perf_counter()

        result = {}
        r_values = self.r_bins[:-1]

        for species_type, rdf_values in rdf_results.items():
            result[species_type] = (r_values, rdf_values)

        if self.show_running_time:
            t12 = time.perf_counter()
            print(f"[DEBUG] 格式化结果: {(t12-t11)*1000:.3f} ms")
            print(f"[DEBUG] === 总时间: {(t12-start_total)*1000:.3f} ms ===\n")

        return result

    def get_species_rdf(self, atom_index, species=None, **kwargs):
        """
        获取单个原子与特定原子类型的RDF

        Parameters
        ----------
        atom_index : int
            中心原子的索引
        species : str, optional
            原子类型，如果为None则返回所有类型
        **kwargs : dict
            其他参数

        Returns
        -------
        dict or tuple
            如果species为None：返回完整的RDF字典
            否则返回特定species的(r_values, rdf_values)元组
        """
        rdf_dict = self.analyze_single_atom(atom_index, **kwargs)

        if species is None:
            return rdf_dict
        else:
            if species not in rdf_dict:
                raise ValueError(f"Species '{species}' not found for atom {atom_index}")
            return rdf_dict[species]

    def _compute_rdf_optimized(
        self, center_pos, reference_coordinates, reference_species
    ):
        """
        优化的RDF计算方法，分别计算与不同原子类型的RDF

        Parameters
        ----------
        center_pos : numpy.ndarray
            中心原子位置，形状为(1, 3)
        reference_coordinates : numpy.ndarray
            参考原子位置数组，形状为(n_atoms, 3)
        reference_species : numpy.ndarray
            参考原子类型数组，形状为(n_atoms,)
        Returns
        -------
        dict
            包含各种原子类型RDF的字典，格式为 {species: rdf_values}
        """
        if self.show_running_time:
            start_rdf = time.perf_counter()
            print(
                f"  [RDF DEBUG] 开始RDF计算，参考原子数: {len(reference_coordinates)}"
            )

        # 使用MDemon的高效距离计算函数
        if self.show_running_time:
            t1 = time.perf_counter()

        box_dims = self.universe.box[:6]

        if self.show_running_time:
            t2 = time.perf_counter()
            print(f"  [RDF DEBUG] 获取盒子维度: {(t2-t1)*1000:.3f} ms")

        # 计算距离
        if self.show_running_time:
            t3 = time.perf_counter()

        distances = distance_array(
            reference=center_pos, configuration=reference_coordinates, box=box_dims
        ).flatten()  # 展平成一维数组

        if self.show_running_time:
            t4 = time.perf_counter()
            print(f"  [RDF DEBUG] 距离矩阵计算: {(t4-t3)*1000:.3f} ms")
            print(f"  [RDF DEBUG] 计算的距离数: {len(distances)}")

        # 移除距离为0的原子（即中心原子本身）
        # 使用小的容差值来处理浮点数精度问题
        if self.show_running_time:
            t5 = time.perf_counter()

        tolerance = 1e-10
        valid_mask = distances > tolerance
        distances = distances[valid_mask]
        reference_species_filtered = reference_species[valid_mask]

        if self.show_running_time:
            t6 = time.perf_counter()
            print(f"  [RDF DEBUG] 过滤距离为0的原子: {(t6-t5)*1000:.3f} ms")
            print(f"  [RDF DEBUG] 有效距离数: {len(distances)}")

        if len(distances) == 0:
            warnings.warn("No reference atoms found for RDF calculation")
            return {}

        # 获取所有唯一的原子类型
        if self.show_running_time:
            t7 = time.perf_counter()

        unique_species = np.unique(reference_species_filtered)

        if self.show_running_time:
            t8 = time.perf_counter()
            print(f"  [RDF DEBUG] 获取唯一物种: {(t8-t7)*1000:.3f} ms")
            print(f"  [RDF DEBUG] 物种类型: {unique_species}")

        rdf_results = {}

        # 分别计算每种原子类型的RDF
        for i, species in enumerate(unique_species):
            if self.show_running_time:
                t_species_start = time.perf_counter()

            species_mask = reference_species_filtered == species
            species_distances = distances[species_mask]

            if self.show_running_time:
                t_mask = time.perf_counter()
                print(
                    f"  [RDF DEBUG] 物种 {species} 掩码筛选: {(t_mask-t_species_start)*1000:.3f} ms"
                )
                print(f"  [RDF DEBUG] 物种 {species} 距离数: {len(species_distances)}")

            if len(species_distances) > 0:
                # 计算直方图
                if self.show_running_time:
                    t_hist_start = time.perf_counter()

                hist, _ = np.histogram(species_distances, bins=self.r_bins)

                if self.show_running_time:
                    t_hist_end = time.perf_counter()
                    print(
                        f"  [RDF DEBUG] 物种 {species} 直方图计算: {(t_hist_end-t_hist_start)*1000:.3f} ms"
                    )

                # 归一化
                if self.show_running_time:
                    t_norm_start = time.perf_counter()

                rdf_values = self._normalize_rdf(hist, len(species_distances))
                rdf_results[species] = rdf_values

                if self.show_running_time:
                    t_norm_end = time.perf_counter()
                    print(
                        f"  [RDF DEBUG] 物种 {species} 归一化: {(t_norm_end-t_norm_start)*1000:.3f} ms"
                    )
                    print(
                        f"  [RDF DEBUG] 物种 {species} 总处理时间: {(t_norm_end-t_species_start)*1000:.3f} ms"
                    )
            else:
                rdf_results[species] = np.zeros(len(self.r_bins) - 1)
                if self.show_running_time:
                    print(f"  [RDF DEBUG] 物种 {species} 无有效距离，填充零值")

        if self.show_running_time:
            end_rdf = time.perf_counter()
            print(f"  [RDF DEBUG] RDF计算总时间: {(end_rdf-start_rdf)*1000:.3f} ms")

        return rdf_results

    def _normalize_rdf(self, hist, n_reference):
        """
        归一化RDF

        使用预计算的球壳体积进行归一化，避免重复计算
        """
        if self.show_running_time:
            start_norm = time.perf_counter()

        # 使用缓存的球壳体积
        if self.show_running_time:
            t1 = time.perf_counter()

        shell_volumes = self._shell_volumes

        if self.show_running_time:
            t2 = time.perf_counter()
            print(f"    [NORM DEBUG] 获取缓存的球壳体积: {(t2-t1)*1000:.3f} ms")

        # 计算数密度
        if self.show_running_time:
            t3 = time.perf_counter()

        if n_reference > 0:
            number_density = hist / shell_volumes
            # 这里需要知道系统的总体积来计算正确的RDF
            # 暂时使用简化的归一化
            rdf_values = number_density / (n_reference / np.sum(shell_volumes))
        else:
            rdf_values = np.zeros_like(hist, dtype=float)

        if self.show_running_time:
            t4 = time.perf_counter()
            print(f"    [NORM DEBUG] 数密度计算和归一化: {(t4-t3)*1000:.3f} ms")
            print(f"    [NORM DEBUG] 归一化总时间: {(t4-start_norm)*1000:.3f} ms")

        return rdf_values

    def _create_analysis_result(self, result_dict, **kwargs):
        """
        创建RDF特定的分析结果对象

        重写基类方法以返回RDFResult而不是AnalysisResult

        Parameters
        ----------
        result_dict : dict
            分析结果字典，键为原子索引，值为按species分组的RDF字典
        **kwargs : dict
            额外的元数据

        Returns
        -------
        RDFResult
            RDF分析结果对象
        """
        atom_indices = list(result_dict.keys())

        # 转换结果数据格式
        species_rdf_data = {}

        for atom_idx in atom_indices:
            rdf_dict = result_dict[atom_idx]
            species_data = {}

            # 提取每个species的RDF值
            for species, (r_values, rdf_values) in rdf_dict.items():
                species_data[species] = rdf_values

            species_rdf_data[atom_idx] = species_data

        # 创建元数据
        metadata = {
            "r_range": self.r_range,
            "n_bins": self.n_bins,
            "timestep": self.universe.timestep,
            "scheduler": self.parallel_manager.scheduler,
            "n_atoms": len(atom_indices),
        }
        metadata.update(kwargs)

        return RDFResult(
            universe=self.universe,
            atom_indices=atom_indices,
            r_bins=self.r_bins,
            species_rdf_data=species_rdf_data,
            metadata=metadata,
        )

    def compute_average_rdf(self, atom_indices=None, species=None, **kwargs):
        """
        计算平均RDF

        Parameters
        ----------
        atom_indices : List[int], optional
            要分析的原子索引列表
        species : str
            原子类型（必须指定，因为不再支持总RDF）
        **kwargs : dict
            其他参数

        Returns
        -------
        tuple
            (r_values, average_rdf) - 距离值和平均RDF值
        """
        if species is None:
            raise ValueError(
                "Must specify species for RDF calculation. Total RDF is no longer supported."
            )

        # 并行分析所有原子
        result = self.analyze_parallel(atom_indices=atom_indices, **kwargs)

        # 收集特定species的RDF数据
        species_rdf_values = []
        for atom_idx in result.atom_indices:
            if atom_idx in result.species_rdf_data:
                atom_species_data = result.species_rdf_data[atom_idx]
                if species in atom_species_data:
                    species_rdf_values.append(atom_species_data[species])

        if not species_rdf_values:
            raise ValueError(f"No RDF data found for species '{species}'")

        # 计算平均值
        species_rdf_values = np.array(species_rdf_values)
        average_rdf = np.mean(species_rdf_values, axis=0)

        return result.r_bins[:-1], average_rdf
