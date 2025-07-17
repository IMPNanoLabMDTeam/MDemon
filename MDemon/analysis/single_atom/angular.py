"""
Angular Distribution Function analysis for single atoms in MDemon

This module provides angular distribution analysis capabilities for single atoms, including:
- Angular distribution calculation between central atom and two neighboring atoms
- Support for a-a-a and a-b-a triplet types only
- Species-specific cutoff radii for different atom pairs
- Parallel processing support via Dask

Classes:
    AngularAnalyzer: Angular distribution analysis for single atoms
    AngularResult: Result container for angular distribution analysis
"""

import time
import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from dask import delayed

from ...lib.distance import distance_array
from .base import AnalysisResult, SingleAtomAnalyzer


class AngularResult(AnalysisResult):
    """
    角分布分析结果容器

    这个类扩展了基础的AnalysisResult，专门用于存储和处理按原子类型分组的角分布分析结果。

    Parameters
    ----------
    atom_indices : List[int]
        分析的原子索引列表
    angle_bins : numpy.ndarray
        角度分箱边界（弧度）
    triplet_angular_data : dict
        按三原子组合类型分组的角分布数据，格式为 {atom_index: {triplet_type: angle_distribution}}
    metadata : dict, optional
        额外的元数据信息
    """

    def __init__(
        self,
        atom_indices,
        angle_bins,
        triplet_angular_data,
        metadata=None,
        universe=None,
    ):
        # 将角分布特定数据添加到基础数据中
        data = triplet_angular_data

        super().__init__(
            analysis_type="angular",
            atom_indices=atom_indices,
            data=data,
            metadata=metadata,
            universe=universe,
        )
        self.angle_bins = angle_bins
        self.triplet_angular_data = triplet_angular_data

    def get_triplet_angular(self, atom_index, triplet_type=None):
        """
        获取特定原子与特定三原子类型的角分布

        Parameters
        ----------
        atom_index : int
            原子索引
        triplet_type : str, optional
            三原子类型（如'C-C-C', 'C-N-C'），如果为None则返回所有类型

        Returns
        -------
        dict or tuple
            如果triplet_type为None：返回 {triplet_type: (angle_values, angular_distribution)} 字典
            否则返回 (angle_values, angular_distribution) 元组
        """
        if atom_index not in self.triplet_angular_data:
            raise ValueError(f"Angular data for atom {atom_index} not available")

        atom_triplet_data = self.triplet_angular_data[atom_index]
        angle_values = self.angle_bins[:-1]  # 使用分箱的左边界作为角度值

        if triplet_type is None:
            # 返回所有类型的角分布
            result = {}
            for tp, angular_vals in atom_triplet_data.items():
                result[tp] = (angle_values, angular_vals)
            return result
        else:
            # 返回特定类型的角分布
            if triplet_type not in atom_triplet_data:
                raise ValueError(
                    f"Triplet type '{triplet_type}' not found for atom {atom_index}"
                )
            return angle_values, atom_triplet_data[triplet_type]

    def get_available_triplet_types(self, atom_index=None):
        """
        获取可用的三原子类型列表

        Parameters
        ----------
        atom_index : int, optional
            特定原子的索引，如果为None则返回所有原子的所有类型

        Returns
        -------
        list
            三原子类型列表
        """
        if atom_index is not None:
            if atom_index not in self.triplet_angular_data:
                return []
            return list(self.triplet_angular_data[atom_index].keys())
        else:
            # 返回所有原子的所有类型
            all_triplet_types = set()
            for atom_data in self.triplet_angular_data.values():
                all_triplet_types.update(atom_data.keys())
            return list(all_triplet_types)

    def get_peak_angle(self, atom_index, triplet_type):
        """
        获取特定原子和三原子类型的峰值角度

        Parameters
        ----------
        atom_index : int
            原子索引
        triplet_type : str
            三原子类型

        Returns
        -------
        float
            峰值角度（弧度）
        """
        angle_values, angular_vals = self.get_triplet_angular(atom_index, triplet_type)
        peak_idx = np.argmax(angular_vals)
        return angle_values[peak_idx]

    def plot(
        self, triplet_types=None, ax=None, show_std=True, angle_unit="degree", **kwargs
    ):
        """
        绘制平均角分布图。默认绘制所有三原子类型的组合图。

        Parameters
        ----------
        triplet_types : list of str, optional
            要绘制的三原子类型列表，例如 `['C-C-C', 'C-N-C']`。
            如果为 `None` (默认)，将绘制所有可用三原子类型的组合图。
        ax : matplotlib.axes.Axes, optional
            用于绘图的matplotlib轴对象。如果为 `None`，将创建一个新的图和轴。
        show_std : bool, optional
            是否显示标准差的阴影区域。默认为 `True`。
        angle_unit : str, optional
            角度单位：'degree' 或 'radian'。默认为 'degree'。
        **kwargs : dict
            传递给 `matplotlib.pyplot.plot` 的其他关键字参数。

        Returns
        -------
        matplotlib.axes.Axes
            包含角分布图的matplotlib轴对象。
        """
        return self.plot_average_angular(
            triplet_types=triplet_types,
            ax=ax,
            show_std=show_std,
            angle_unit=angle_unit,
            **kwargs,
        )

    def plot_average_angular(
        self,
        triplet_types=None,
        plot_type="combined",
        show_std=True,
        angle_unit="degree",
        **kwargs,
    ):
        """
        计算并绘制三原子类型之间的平均角分布。

        Parameters
        ----------
        triplet_types : list of str, optional
            要绘制的三原子类型列表, e.g., `['C-C-C', 'C-N-C']`。
            如果为 `None`, 将使用所有可能的三原子类型。
        plot_type : {'combined', 'separate'}, optional
            绘图类型。'combined' 将所有角分布绘制在一张图上，'separate'
            为每个三原子类型创建一个单独的图。默认为 'combined'。
        show_std : bool, optional
            是否显示标准差的阴影区域。默认为 `True`。
        angle_unit : str, optional
            角度单位：'degree' 或 'radian'。默认为 'degree'。
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

        averaged_data = self._calculate_average_angular(triplet_types)

        if not averaged_data:
            warnings.warn("No data available for plotting.")
            return

        if plot_type == "combined":
            fig, ax = plt.subplots(figsize=(12, 8))
            colors = plt.cm.tab10(np.linspace(0, 1, len(averaged_data)))

            # 过滤掉可能冲突的参数
            plot_kwargs = {k: v for k, v in kwargs.items() if k not in ["ax"]}

            for i, (triplet_type, data) in enumerate(averaged_data.items()):
                angle_vals = data["angle_values"]
                avg_angular = data["avg_angular"]
                n_samples = data["n_samples"]

                # 转换角度单位
                if angle_unit == "degree":
                    angle_vals = np.degrees(angle_vals)
                    xlabel = "Angle (degrees)"
                else:
                    xlabel = "Angle (radians)"

                ax.plot(
                    angle_vals,
                    avg_angular,
                    color=colors[i],
                    linewidth=2,
                    label=f"Triplet {triplet_type} (n={n_samples})",
                    **plot_kwargs,
                )

            ax.set_xlabel(xlabel)
            ax.set_ylabel("P(θ)")
            ax.set_title("Angular Distribution Comparison - All Triplet Types")
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            return ax

        elif plot_type == "separate":
            axes = []
            # 过滤掉可能冲突的参数
            plot_kwargs = {k: v for k, v in kwargs.items() if k not in ["ax"]}

            for triplet_type, data in averaged_data.items():
                fig, ax = plt.subplots(figsize=(10, 6))

                angle_vals = data["angle_values"]
                avg_angular = data["avg_angular"]
                std_angular = data["std_angular"]
                n_samples = data["n_samples"]

                # 转换角度单位
                if angle_unit == "degree":
                    angle_vals = np.degrees(angle_vals)
                    xlabel = "Angle (degrees)"
                else:
                    xlabel = "Angle (radians)"

                ax.plot(
                    angle_vals,
                    avg_angular,
                    linewidth=2,
                    label=f"Triplet {triplet_type} (n={n_samples})",
                    **plot_kwargs,
                )
                if show_std:
                    ax.fill_between(
                        angle_vals,
                        avg_angular - std_angular,
                        avg_angular + std_angular,
                        alpha=0.3,
                    )

                ax.set_xlabel(xlabel)
                ax.set_ylabel("P(θ)")
                ax.set_title(f"Angular Distribution: Triplet {triplet_type}")
                ax.legend()
                ax.grid(True, alpha=0.3)
                axes.append(ax)
            return axes

    def _calculate_average_angular(self, triplet_types=None):
        """计算指定三原子类型的平均角分布、标准差和样本数。"""
        if self.universe is None:
            raise ValueError("Universe object is required for triplet-based averaging.")

        available_triplet_types = self.get_available_triplet_types()

        if triplet_types is None:
            triplet_types = available_triplet_types

        averaged_data = {}
        analyzed_atom_indices = set(self.atom_indices)

        for triplet_type in triplet_types:
            if triplet_type not in available_triplet_types:
                continue

            angular_values_list = []

            for atom_idx in analyzed_atom_indices:
                if (
                    atom_idx in self.triplet_angular_data
                    and triplet_type in self.triplet_angular_data[atom_idx]
                ):
                    angular_values_list.append(
                        self.triplet_angular_data[atom_idx][triplet_type]
                    )

            if angular_values_list:
                angular_array = np.array(angular_values_list)
                avg_angular = np.mean(angular_array, axis=0)
                avg_angular = avg_angular / np.sum(avg_angular)
                std_angular = np.std(angular_array, axis=0)

                averaged_data[triplet_type] = {
                    "angle_values": self.angle_bins[:-1],
                    "avg_angular": avg_angular,
                    "std_angular": std_angular,
                    "n_samples": len(angular_values_list),
                }
        self.averaged_data = averaged_data
        return averaged_data

    def plot_triplet_comparison(
        self, atom_index, ax=None, angle_unit="degree", **kwargs
    ):
        """
        绘制单个原子与不同三原子类型的角分布比较图

        Parameters
        ----------
        atom_index : int
            原子索引
        ax : matplotlib.axes.Axes, optional
            绘图轴
        angle_unit : str, optional
            角度单位：'degree' 或 'radian'。默认为 'degree'。
        **kwargs : dict
            传递给plot的额外参数
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib is required for plotting")

        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))

        if atom_index not in self.triplet_angular_data:
            raise ValueError(f"Angular data for atom {atom_index} not available")

        # 绘制各个三原子类型的角分布
        triplet_data = self.get_triplet_angular(atom_index)
        for triplet_type, (angle_vals, angular_vals) in triplet_data.items():
            # 转换角度单位
            if angle_unit == "degree":
                angle_vals = np.degrees(angle_vals)
                xlabel = "Angle (degrees)"
            else:
                xlabel = "Angle (radians)"

            ax.plot(
                angle_vals, angular_vals, label=f"Triplet {triplet_type}", alpha=0.8
            )

        ax.set_xlabel(xlabel)
        ax.set_ylabel("P(θ)")
        ax.set_title(f"Angular Distribution for Atom {atom_index} by Triplet Type")
        ax.legend()
        ax.grid(True, alpha=0.3)

        return ax


class AngularAnalyzer(SingleAtomAnalyzer):
    """
    角分布分析器

    这个类提供了单原子角分布分析的核心功能，包括：
    - 单个原子的角分布计算
    - 并行处理多个原子的角分布
    - 支持a-a-a和a-b-a两种三原子类型
    - 不同原子对的特定截断半径

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
    angle_range : tuple, optional
        角度范围（弧度） (angle_min, angle_max)，Default: (0.0, np.pi)
    n_bins : int, optional
        角度分箱数量，Default: 90
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
        angle_range=(0.0, np.pi),
        n_bins=90,
        scheduler="threads",
        enable_spatial_subdivision=True,
        **kwargs,
    ):
        # 验证截断半径参数
        if cutoff_radii is None:
            raise ValueError(
                "cutoff_radii is required and must be a dictionary with format "
                "{'1-1': 2.5, '1-2': 3.0, ...} where keys represent triplet types "
                "(e.g., '1-1' for 1-1-1 angles, '1-2' for 2-1-2 angles)"
            )

        if not isinstance(cutoff_radii, dict):
            raise TypeError("cutoff_radii must be a dictionary")

        self.cutoff_radii = cutoff_radii

        # 使用固定的10埃作为空间预分割的截断半径
        r_cutoff = 10.0 if enable_spatial_subdivision else None

        super().__init__(
            universe,
            atom_selection=atom_selection,
            scheduler=scheduler,
            enable_spatial_subdivision=enable_spatial_subdivision,
            r_cutoff=r_cutoff,
            **kwargs,
        )

        self.angle_range = angle_range
        self.n_bins = n_bins
        self.angle_bins = np.linspace(angle_range[0], angle_range[1], n_bins + 1)
        self.dangle = self.angle_bins[1] - self.angle_bins[0]

        # 在初始化时获取系统中的所有唯一物种类型
        self.unique_species = np.unique(self.universe.atoms.species)

    def _get_cutoff_for_triplet(self, center_species, neighbor_species):
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
        分析单个原子的角分布，返回按三原子类型分组的结果

        Parameters
        ----------
        atom_index : int
            中心原子的索引
        **kwargs : dict
            其他参数

        Returns
        -------
        dict
            按三原子类型分组的角分布字典，格式为 {triplet_type: (angle_values, angular_distribution)}
        """
        if self.show_running_time:
            start_total = time.perf_counter()
            print(f"\n[ANGULAR DEBUG] === Analyzing atom {atom_index} ===")

        # 获取原子坐标和原子类型
        if self.show_running_time:
            t1 = time.perf_counter()

        atoms = self.universe.atoms
        all_coordinates = atoms.coordinate
        all_species = atoms.species
        center_species = all_species[atom_index]

        if self.show_running_time:
            t2 = time.perf_counter()
            print(f"[ANGULAR DEBUG] 获取原子坐标和类型: {(t2-t1)*1000:.3f} ms")
            print(f"[ANGULAR DEBUG] 中心原子类型: {center_species}")

        # 使用空间预分割获取相关原子（如果启用）
        if self.show_running_time:
            t3 = time.perf_counter()

        relevant_atom_indices, combined_mask = self._get_relevant_atoms_for_analysis(
            atom_index
        )

        if self.show_running_time:
            t4 = time.perf_counter()
            print(f"[ANGULAR DEBUG] 空间预分割获取相关原子: {(t4-t3)*1000:.3f} ms")

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
            print(f"[ANGULAR DEBUG] 筛选坐标和物种信息: {(t6-t5)*1000:.3f} ms")
            print(f"[ANGULAR DEBUG] 筛选后原子数: {len(coordinates)}")

        # 计算角分布
        if self.show_running_time:
            t7 = time.perf_counter()

        angular_results = self._compute_angular_distribution(
            center_index_in_filtered, coordinates, species
        )

        if self.show_running_time:
            t8 = time.perf_counter()
            print(f"[ANGULAR DEBUG] 角分布计算: {(t8-t7)*1000:.3f} ms")
            print(f"[ANGULAR DEBUG] 计算得到的三原子类型数: {len(angular_results)}")

        # 重新格式化结果
        if self.show_running_time:
            t9 = time.perf_counter()

        result = {}
        angle_values = self.angle_bins[:-1]

        for triplet_type, angular_values in angular_results.items():
            result[triplet_type] = (angle_values, angular_values)

        if self.show_running_time:
            t10 = time.perf_counter()
            print(f"[ANGULAR DEBUG] 格式化结果: {(t10-t9)*1000:.3f} ms")
            print(f"[ANGULAR DEBUG] === 总时间: {(t10-start_total)*1000:.3f} ms ===\n")

        return result

    def _compute_angular_distribution(self, center_index, coordinates, species):
        """
        计算角分布

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
            包含各种三原子类型角分布的字典，格式为 {triplet_type: angular_distribution}
        """
        if self.show_running_time:
            start_angular = time.perf_counter()
            print(f"  [ANGULAR DEBUG] 开始角分布计算")

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
            print(f"  [ANGULAR DEBUG] 距离计算: {(t2-t1)*1000:.3f} ms")

        # 移除中心原子本身
        tolerance = 1e-10
        valid_mask = distances > tolerance
        distances = distances[valid_mask]
        nearby_coordinates = coordinates[valid_mask]
        nearby_species = species[valid_mask]

        if self.show_running_time:
            t3 = time.perf_counter()
            print(f"  [ANGULAR DEBUG] 过滤中心原子: {(t3-t2)*1000:.3f} ms")
            print(f"  [ANGULAR DEBUG] 邻近原子数: {len(nearby_coordinates)}")

        angular_results = {}

        # 对于每个物种类型，计算对应的角分布
        if self.show_running_time:
            t4 = time.perf_counter()

        for sp in self.unique_species:
            # 构建中心-邻近原子键
            triplet_key = f"{center_species}-{sp}"

            try:
                # 获取对应的截断半径
                cutoff = self._get_cutoff_for_triplet(center_species, sp)

                # 筛选出在截断半径内且类型为sp的原子
                sp_mask = (nearby_species == sp) & (distances <= cutoff)

                if np.sum(sp_mask) >= 2:  # 至少需要两个邻近原子
                    # 获取目标原子的坐标
                    target_coords = nearby_coordinates[sp_mask]

                    # 计算从邻近原子到中心原子的向量
                    vectors = target_coords - center_pos  # shape: (n_neighbors, 3)

                    # 归一化向量
                    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
                    normalized_vectors = vectors / (norms + 1e-12)  # 避免除零

                    # 计算向量点积矩阵
                    dot_products = np.dot(
                        normalized_vectors, normalized_vectors.T
                    )  # shape: (n_neighbors, n_neighbors)

                    # 取上三角部分（避免重复和自己与自己的点积）
                    upper_triangle_indices = np.triu_indices_from(dot_products, k=1)
                    cos_angles = dot_products[upper_triangle_indices]

                    # 处理数值精度问题并计算角度
                    cos_angles = np.clip(cos_angles, -1.0, 1.0)
                    angles = np.arccos(cos_angles)

                    if len(angles) > 0:
                        # 计算角分布直方图
                        hist, _ = np.histogram(angles, bins=self.angle_bins)

                        # 归一化
                        angular_distribution = hist / (np.sum(hist) + 1e-12)  # 避免除零

                        # 直接使用triplet_key作为结果的键
                        angular_results[triplet_key] = angular_distribution

                        if self.show_running_time:
                            print(
                                f"  [ANGULAR DEBUG] {triplet_key}: {len(angles)} 个角度"
                            )

            except ValueError:
                # 如果没有找到对应的截断半径配置，跳过此物种
                continue

        if self.show_running_time:
            t5 = time.perf_counter()
            print(f"  [ANGULAR DEBUG] 角分布计算: {(t5-t4)*1000:.3f} ms")
            print(
                f"  [ANGULAR DEBUG] 角分布计算总时间: {(t5-start_angular)*1000:.3f} ms"
            )

        return angular_results

    def _create_analysis_result(self, result_dict, **kwargs):
        """
        创建角分布特定的分析结果对象

        重写基类方法以返回AngularResult而不是AnalysisResult

        Parameters
        ----------
        result_dict : dict
            分析结果字典，键为原子索引，值为按triplet_type分组的角分布字典
        **kwargs : dict
            额外的元数据

        Returns
        -------
        AngularResult
            角分布分析结果对象
        """
        atom_indices = list(result_dict.keys())

        # 转换结果数据格式
        triplet_angular_data = {}

        for atom_idx in atom_indices:
            angular_dict = result_dict[atom_idx]
            triplet_data = {}

            # 提取每个triplet_type的角分布值
            for triplet_type, (angle_values, angular_values) in angular_dict.items():
                triplet_data[triplet_type] = angular_values

            triplet_angular_data[atom_idx] = triplet_data

        # 创建元数据
        metadata = {
            "angle_range": self.angle_range,
            "n_bins": self.n_bins,
            "cutoff_radii": self.cutoff_radii,
            "timestep": self.universe.timestep,
            "scheduler": self.parallel_manager.scheduler,
            "n_atoms": len(atom_indices),
        }
        metadata.update(kwargs)

        return AngularResult(
            universe=self.universe,
            atom_indices=atom_indices,
            angle_bins=self.angle_bins,
            triplet_angular_data=triplet_angular_data,
            metadata=metadata,
        )

    def compute_average_angular(self, atom_indices=None, triplet_type=None, **kwargs):
        """
        计算平均角分布

        Parameters
        ----------
        atom_indices : List[int], optional
            要分析的原子索引列表
        triplet_type : str
            三原子类型（必须指定）
        **kwargs : dict
            其他参数

        Returns
        -------
        tuple
            (angle_values, average_angular) - 角度值和平均角分布值
        """
        if triplet_type is None:
            raise ValueError(
                "Must specify triplet_type for angular distribution calculation."
            )

        # 并行分析所有原子
        result = self.analyze_parallel(atom_indices=atom_indices, **kwargs)

        # 收集特定triplet_type的角分布数据
        triplet_angular_values = []
        for atom_idx in result.atom_indices:
            if atom_idx in result.triplet_angular_data:
                atom_triplet_data = result.triplet_angular_data[atom_idx]
                if triplet_type in atom_triplet_data:
                    triplet_angular_values.append(atom_triplet_data[triplet_type])

        if not triplet_angular_values:
            raise ValueError(f"No angular data found for triplet type '{triplet_type}'")

        # 计算平均值
        triplet_angular_values = np.array(triplet_angular_values)
        average_angular = np.mean(triplet_angular_values, axis=0)

        return result.angle_bins[:-1], average_angular
