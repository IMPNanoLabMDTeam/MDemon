"""
Cylindrical selection utilities for MDemon

This module provides functionality to select atoms within a cylindrical region
or cylindrical shell defined by an axis-parallel line and radius range.

Classes:
    CylindricalSelector: Select atoms within a cylindrical region or cylindrical shell
"""

import warnings
from typing import Optional, Tuple, Union

import numpy as np

from ..lib.distance import distance_array


class CylindricalSelector:
    """
    圆柱形和圆柱壳原子选择器

    用于筛选出以某根平行于轴(x, y, 或z)的直线为圆柱轴的指定半径范围内的所有原子。
    支持圆柱选择（rmin=0）和圆柱壳选择（rmin>0）。

    Parameters
    ----------
    universe : MDemon.Universe
        分析的宇宙对象
    axis : {'x', 'y', 'z'}
        圆柱轴方向，平行于指定的坐标轴
    center_point : array_like, shape (3,)
        圆柱轴上的一个参考点，格式为 [x, y, z]
    rmax : float, optional
        圆柱外半径。如果提供了radius参数，则此参数被忽略
    rmin : float, optional
        圆柱内半径，默认为0.0（实心圆柱）
    radius : float, optional
        圆柱半径，等价于rmax，rmin=0。为保持向后兼容性保留此参数。
        如果同时提供radius和rmax，优先使用rmax/rmin组合
    axis_range : tuple, optional
        沿圆柱轴方向的选择范围 (min_val, max_val)。
        如果为None，则不限制轴向范围

    Examples
    --------
    >>> # 选择以z轴为圆柱轴，通过点(0,0,0)，半径为5的圆柱内的原子（实心圆柱）
    >>> selector = CylindricalSelector(universe, axis='z',
    ...                                center_point=[0, 0, 0], rmax=5.0)
    >>> atom_indices = selector.select()

    >>> # 选择圆柱壳：内径2，外径5
    >>> selector = CylindricalSelector(universe, axis='z',
    ...                                center_point=[0, 0, 0],
    ...                                rmin=2.0, rmax=5.0)
    >>> atom_indices = selector.select()

    >>> # 向后兼容：使用radius参数
    >>> selector = CylindricalSelector(universe, axis='z',
    ...                                center_point=[0, 0, 0], radius=5.0)
    >>> atom_indices = selector.select()

    >>> # 选择以x轴为圆柱轴，通过点(10,5,5)，内径3外径8，x范围在[8,12]的圆柱壳内的原子
    >>> selector = CylindricalSelector(universe, axis='x',
    ...                                center_point=[10, 5, 5],
    ...                                rmin=3.0, rmax=8.0,
    ...                                axis_range=(8, 12))
    >>> atom_indices = selector.select()
    """

    def __init__(
        self,
        universe,
        axis: str,
        center_point: Union[list, tuple, np.ndarray],
        rmax: Optional[float] = None,
        rmin: float = 0.0,
        radius: Optional[float] = None,
        axis_range: Optional[Tuple[float, float]] = None,
    ):
        self.universe = universe
        self.axis = axis.lower()
        self.center_point = np.array(center_point, dtype=np.float64)
        self.axis_range = axis_range

        # 处理半径参数的优先级和兼容性
        self._setup_radius_parameters(rmax, rmin, radius)

        # 验证输入参数
        self._validate_inputs()

        # 确定轴向索引和2D投影索引
        self._setup_axis_mapping()

    def _setup_radius_parameters(
        self, rmax: Optional[float], rmin: float, radius: Optional[float]
    ):
        """设置半径参数，处理向后兼容性"""
        if rmax is not None:
            # 优先使用rmax/rmin组合
            self.rmax = rmax
            self.rmin = rmin
            if radius is not None:
                warnings.warn(
                    "Both 'rmax' and 'radius' parameters provided. Using 'rmax/rmin' combination. "
                    "'radius' parameter is ignored.",
                    UserWarning,
                )
        elif radius is not None:
            # 使用radius参数（向后兼容）
            self.rmax = radius
            self.rmin = 0.0
            if rmin != 0.0:
                warnings.warn(
                    "Both 'radius' and 'rmin' parameters provided. Setting rmin=0.0 for compatibility. "
                    "Use 'rmax' and 'rmin' parameters for cylindrical shell selection.",
                    UserWarning,
                )
        else:
            raise ValueError("Either 'rmax' or 'radius' parameter must be provided")

    def _validate_inputs(self):
        """验证输入参数的有效性"""
        if self.axis not in ["x", "y", "z"]:
            raise ValueError(f"axis must be 'x', 'y', or 'z', got '{self.axis}'")

        if self.center_point.shape != (3,):
            raise ValueError(
                f"center_point must have shape (3,), got {self.center_point.shape}"
            )

        if self.rmax <= 0:
            raise ValueError(f"rmax must be positive, got {self.rmax}")

        if self.rmin < 0:
            raise ValueError(f"rmin must be non-negative, got {self.rmin}")

        if self.rmin >= self.rmax:
            raise ValueError(f"rmin ({self.rmin}) must be less than rmax ({self.rmax})")

        if self.axis_range is not None:
            if len(self.axis_range) != 2:
                raise ValueError("axis_range must be a tuple of length 2")
            if self.axis_range[0] >= self.axis_range[1]:
                raise ValueError("axis_range[0] must be less than axis_range[1]")

    def _setup_axis_mapping(self):
        """设置轴向映射"""
        # 轴向索引映射
        axis_map = {"x": 0, "y": 1, "z": 2}
        self.axis_index = axis_map[self.axis]

        print(f"Cylindrical selection setup:")
        print(f"  Axis: {self.axis} (index {self.axis_index})")
        print(f"  Center point: {self.center_point}")
        if self.rmin > 0:
            print(f"  Radial range: {self.rmin} - {self.rmax} (cylindrical shell)")
        else:
            print(f"  Radius: {self.rmax} (solid cylinder)")
        if self.axis_range is not None:
            print(f"  Axis range: {self.axis_range}")

    @property
    def is_cylindrical_shell(self) -> bool:
        """判断是否为圆柱壳选择（rmin > 0）"""
        return self.rmin > 0

    @property
    def radius(self) -> float:
        """向后兼容性属性，返回rmax"""
        return self.rmax

    def select(self, return_mask: bool = False) -> Union[np.ndarray, np.ndarray]:
        """
        执行圆柱形选择

        Parameters
        ----------
        return_mask : bool, optional
            如果为True，返回布尔掩码；如果为False，返回原子索引数组。
            默认为False

        Returns
        -------
        numpy.ndarray
            如果return_mask=True：布尔掩码，形状为(n_atoms,)
            如果return_mask=False：选中的原子索引数组
        """
        # 获取所有原子坐标
        all_coordinates = self.universe.atoms.coordinate
        n_atoms = len(all_coordinates)

        print(f"Starting cylindrical selection for {n_atoms} atoms...")

        # 1. 轴向筛选（如果指定了axis_range）
        if self.axis_range is not None:
            axis_coords = all_coordinates[:, self.axis_index]
            axis_mask = (axis_coords >= self.axis_range[0]) & (
                axis_coords <= self.axis_range[1]
            )

            if not np.any(axis_mask):
                warnings.warn("No atoms found within the specified axis range")
                if return_mask:
                    return np.zeros(n_atoms, dtype=bool)
                else:
                    return np.array([], dtype=int)

            # 应用轴向筛选
            filtered_coordinates = all_coordinates[axis_mask]
            axis_filtered_indices = np.where(axis_mask)[0]
            print(f"After axis filtering: {len(filtered_coordinates)} atoms")
        else:
            filtered_coordinates = all_coordinates
            axis_filtered_indices = np.arange(n_atoms)

        # 2. 准备中心点坐标（保持3D格式）
        center_3d = self.center_point.reshape(1, -1)

        print(f"Filtered coordinates shape: {filtered_coordinates.shape}")
        print(f"Center point: {center_3d.flatten()}")

        # 3. 计算圆柱形距离（指定轴向坐标设为0）
        distances_radial = self._calculate_cylindrical_distances(
            center_3d, filtered_coordinates
        )

        # 4. 应用半径筛选 - 支持圆柱壳
        if self.is_cylindrical_shell:
            # 圆柱壳：rmin <= distance <= rmax
            radius_mask = (distances_radial >= self.rmin) & (
                distances_radial <= self.rmax
            )
            print(f"Cylindrical shell selection: {self.rmin} ≤ r ≤ {self.rmax}")
        else:
            # 实心圆柱：distance <= rmax
            radius_mask = distances_radial <= self.rmax
            print(f"Solid cylinder selection: r ≤ {self.rmax}")

        # 获取最终选中的原子索引
        selected_indices_in_filtered = np.where(radius_mask)[0]
        final_selected_indices = axis_filtered_indices[selected_indices_in_filtered]

        print(f"Final selection: {len(final_selected_indices)} atoms")

        # 5. 返回结果
        if return_mask:
            # 创建完整的布尔掩码
            full_mask = np.zeros(n_atoms, dtype=bool)
            full_mask[final_selected_indices] = True
            return full_mask
        else:
            return final_selected_indices

    def _calculate_cylindrical_distances(
        self, center_3d: np.ndarray, coordinates_3d: np.ndarray
    ) -> np.ndarray:
        """
        计算圆柱形距离

        通过将指定轴向的坐标设为0来计算径向距离，
        保持3D坐标格式以确保后端兼容性和周期性边界条件支持

        Parameters
        ----------
        center_3d : numpy.ndarray, shape (1, 3)
            3D中心点坐标
        coordinates_3d : numpy.ndarray, shape (n, 3)
            3D原子坐标数组

        Returns
        -------
        numpy.ndarray, shape (n,)
            径向距离数组
        """
        # 复制坐标以避免修改原始数据
        modified_center = center_3d.copy()
        modified_coordinates = coordinates_3d.copy()

        # 将指定轴向的坐标设为0
        modified_center[:, self.axis_index] = 0.0
        modified_coordinates[:, self.axis_index] = 0.0

        # 获取盒子信息用于周期性边界条件
        box_dims = self.universe.box[:6] if hasattr(self.universe, "box") else None

        # 使用distance_array计算距离（考虑周期性边界条件）
        distance_matrix = distance_array(
            reference=modified_center, configuration=modified_coordinates, box=box_dims
        )

        # 返回一维距离数组
        return distance_matrix.flatten()

    def get_selection_info(self) -> dict:
        """
        获取选择器的配置信息

        Returns
        -------
        dict
            包含选择器配置的字典
        """
        info = {
            "axis": self.axis,
            "axis_index": self.axis_index,
            "center_point": self.center_point.tolist(),
            "rmax": self.rmax,
            "rmin": self.rmin,
            "is_cylindrical_shell": self.is_cylindrical_shell,
            "axis_range": self.axis_range,
            "total_atoms": len(self.universe.atoms),
        }
        return info

    def print_selection_info(self):
        """打印选择器的配置信息"""
        info = self.get_selection_info()

        print("\n=== Cylindrical Selector Configuration ===")
        print(
            f"Cylindrical axis: {info['axis']} (coordinate index {info['axis_index']})"
        )
        print(f"Center point: {info['center_point']}")

        if info["is_cylindrical_shell"]:
            print(f"Cylindrical shell: rmin={info['rmin']}, rmax={info['rmax']}")
        else:
            print(f"Solid cylinder: radius={info['rmax']}")

        if info["axis_range"] is not None:
            print(f"Axis range: {info['axis_range']}")
        else:
            print("Axis range: No limit")
        print(f"Total atoms in universe: {info['total_atoms']}")
        print("=" * 45)

    def select_with_statistics(
        self, return_mask: bool = False
    ) -> Tuple[np.ndarray, dict]:
        """
        执行选择并返回统计信息

        Parameters
        ----------
        return_mask : bool, optional
            如果为True，返回布尔掩码；如果为False，返回原子索引数组

        Returns
        -------
        Tuple[numpy.ndarray, dict]
            (selected_atoms, statistics) - 选中的原子和统计信息
        """
        selected_atoms = self.select(return_mask=return_mask)

        # 计算统计信息
        if return_mask:
            n_selected = np.sum(selected_atoms)
        else:
            n_selected = len(selected_atoms)

        total_atoms = len(self.universe.atoms)
        selection_ratio = n_selected / total_atoms if total_atoms > 0 else 0

        statistics = {
            "total_atoms": total_atoms,
            "selected_atoms": n_selected,
            "selection_ratio": selection_ratio,
            "axis": self.axis,
            "rmax": self.rmax,
            "rmin": self.rmin,
            "is_cylindrical_shell": self.is_cylindrical_shell,
            "center_point": self.center_point.tolist(),
            "axis_range": self.axis_range,
        }

        return selected_atoms, statistics

    def print_selection_statistics(
        self, selected_atoms: Union[np.ndarray, None] = None
    ):
        """
        打印选择统计信息

        Parameters
        ----------
        selected_atoms : numpy.ndarray, optional
            选中的原子（索引数组或布尔掩码）。如果为None，会重新执行选择
        """
        if selected_atoms is None:
            selected_atoms, stats = self.select_with_statistics()
        else:
            # 计算统计信息
            if selected_atoms.dtype == bool:
                n_selected = np.sum(selected_atoms)
            else:
                n_selected = len(selected_atoms)

            total_atoms = len(self.universe.atoms)
            stats = {
                "total_atoms": total_atoms,
                "selected_atoms": n_selected,
                "selection_ratio": n_selected / total_atoms if total_atoms > 0 else 0,
                "axis": self.axis,
                "rmax": self.rmax,
                "rmin": self.rmin,
                "is_cylindrical_shell": self.is_cylindrical_shell,
                "center_point": self.center_point.tolist(),
                "axis_range": self.axis_range,
            }

        print("\n=== Cylindrical Selection Statistics ===")
        print(f"Total atoms: {stats['total_atoms']}")
        print(f"Selected atoms: {stats['selected_atoms']}")
        print(
            f"Selection ratio: {stats['selection_ratio']:.4f} ({stats['selection_ratio']*100:.2f}%)"
        )
        print(f"Cylinder axis: {stats['axis']}")

        if stats["is_cylindrical_shell"]:
            print(f"Cylindrical shell: rmin={stats['rmin']}, rmax={stats['rmax']}")
        else:
            print(f"Solid cylinder: radius={stats['rmax']}")

        print(f"Center point: {stats['center_point']}")
        if stats["axis_range"] is not None:
            print(f"Axis range: {stats['axis_range']}")
        print("=" * 42)
