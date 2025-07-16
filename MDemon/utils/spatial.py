"""
Spatial subdivision utilities for efficient neighbor searching

This module provides spatial subdivision functionality to accelerate RDF and other
pair-based calculations by dividing the simulation box into sub-regions and only
considering neighboring regions for each calculation.

Classes:
    SpatialSubdivision: Main class for spatial subdivision of atoms
"""

import warnings
from typing import Dict, List, Optional, Tuple, Union

import numpy as np


class SpatialSubdivision:
    """
    空间预分割类，用于高效的邻近原子搜索

    这个类将模拟盒子划分为子空间网格，每个子空间的尺寸至少为截断半径，
    且是盒子边长的整数倍。通过这种预分割，可以显著加速RDF等需要邻近原子
    搜索的计算。

    Parameters
    ----------
    universe : MDemon.Universe
        分析的宇宙对象
    r_cutoff : float
        截断半径，每个子空间的边长至少为此值
    use_masks : bool, optional
        是否生成和存储布尔掩码，默认为False。对于大体系，建议设为False以节省内存
    """

    def __init__(self, universe, r_cutoff: float, use_masks: bool = False):
        self.universe = universe
        self.r_cutoff = r_cutoff
        self.use_masks = use_masks

        # 获取盒子参数 - 参考rdf.py中的用法
        self.box = universe.box[:6]  # [a, b, c, alpha, beta, gamma]

        # 检查是否为正交盒子
        if not self._is_orthogonal_box():
            raise NotImplementedError("Currently only orthogonal boxes are supported")

        # 盒子边长
        self.box_lengths = self.box[:3]

        # 计算子空间划分
        self._calculate_subdivision()

        # 生成原子到子空间的映射
        self._generate_atom_subdivision_mapping()

    def _is_orthogonal_box(self) -> bool:
        """检查是否为正交盒子"""
        angles = self.box[3:6]
        return np.allclose(angles, [90.0, 90.0, 90.0], atol=1e-6)

    def _calculate_subdivision(self):
        """
        计算最优的子空间划分

        为每个方向选择满足条件的最小子空间数量：
        - 子空间尺寸 >= r_cutoff
        - 子空间尺寸 = box_length / n_divisions (n_divisions为整数)
        """
        self.n_subdivisions = np.zeros(3, dtype=int)
        self.subdivision_lengths = np.zeros(3)

        for i, box_length in enumerate(self.box_lengths):
            if box_length <= 0:
                raise ValueError(f"Invalid box length: {box_length}")

            # 计算满足条件的最大子空间数量
            max_n_div = int(np.floor(box_length / self.r_cutoff))

            # 确保至少有一个子空间
            if max_n_div < 1:
                max_n_div = 1
                warnings.warn(
                    f"Box length ({box_length:.3f}) is smaller than cutoff "
                    f"({self.r_cutoff:.3f}) in dimension {i}. "
                    "Using single subdivision."
                )

            self.n_subdivisions[i] = max_n_div
            self.subdivision_lengths[i] = box_length / max_n_div

        # 总的子空间数量
        self.total_subdivisions = np.prod(self.n_subdivisions)

        print(f"Spatial subdivision info:")
        print(f"  Box lengths: {self.box_lengths}")
        print(f"  Cutoff radius: {self.r_cutoff}")
        print(f"  Subdivisions: {self.n_subdivisions}")
        print(f"  Subdivision lengths: {self.subdivision_lengths}")
        print(f"  Total subdivisions: {self.total_subdivisions}")

    def _generate_atom_subdivision_mapping(self):
        """生成原子到子空间的映射"""
        # 获取原子坐标 - 参考rdf.py中的用法
        coordinates = self.universe.atoms.coordinate
        n_atoms = len(coordinates)

        # 计算每个原子所属的子空间索引
        self.atom_subdivision_indices = self._get_subdivision_indices(coordinates)

        # 根据use_masks参数决定是否生成掩码
        if self.use_masks:
            self.subdivision_masks = {}
        else:
            self.subdivision_masks = None

        self.subdivision_atom_lists = {}

        for subdivision_idx in range(self.total_subdivisions):
            mask = self.atom_subdivision_indices == subdivision_idx
            atom_indices = np.where(mask)[0]

            if len(atom_indices) > 0:  # 只存储有原子的子空间
                if self.use_masks:
                    self.subdivision_masks[subdivision_idx] = mask
                self.subdivision_atom_lists[subdivision_idx] = atom_indices

    def _get_subdivision_indices(self, coordinates: np.ndarray) -> np.ndarray:
        """
        计算每个原子所属的子空间索引

        Parameters
        ----------
        coordinates : np.ndarray
            原子坐标数组，形状为 (n_atoms, 3)

        Returns
        -------
        np.ndarray
            每个原子对应的子空间索引，形状为 (n_atoms,)
        """
        # 将坐标标准化到[0, n_subdivisions)范围
        normalized_coords = coordinates / self.subdivision_lengths

        # 转换为整数坐标（子空间索引）
        subdivision_coords = np.floor(normalized_coords).astype(int)

        # 处理边界情况（确保在有效范围内）
        for i in range(3):
            subdivision_coords[:, i] = np.clip(
                subdivision_coords[:, i], 0, self.n_subdivisions[i] - 1
            )

        # 将3D索引转换为1D索引
        indices = self._coords_3d_to_1d(
            subdivision_coords[:, 0], subdivision_coords[:, 1], subdivision_coords[:, 2]
        )

        return indices

    def _coords_3d_to_1d(
        self, i: np.ndarray, j: np.ndarray, k: np.ndarray
    ) -> np.ndarray:
        """
        将3D子空间坐标转换为1D索引（向量化版本）

        Parameters
        ----------
        i, j, k : np.ndarray
            3D坐标数组

        Returns
        -------
        np.ndarray
            1D索引数组
        """
        return (
            i * self.n_subdivisions[1] * self.n_subdivisions[2]
            + j * self.n_subdivisions[2]
            + k
        )

    def _index_1d_to_3d(self, index: int) -> Tuple[int, int, int]:
        """
        将1D索引转换为3D子空间坐标

        Parameters
        ----------
        index : int
            1D索引

        Returns
        -------
        Tuple[int, int, int]
            (i, j, k) 3D坐标
        """
        k = index % self.n_subdivisions[2]
        temp = index // self.n_subdivisions[2]
        j = temp % self.n_subdivisions[1]
        i = temp // self.n_subdivisions[1]
        return i, j, k

    def get_neighboring_subdivisions(
        self, subdivision_index: int, include_self: bool = True
    ) -> List[int]:
        """
        获取相邻的子空间索引（考虑周期性边界条件）

        Parameters
        ----------
        subdivision_index : int
            中心子空间的索引
        include_self : bool, optional
            是否包含自身，默认为True

        Returns
        -------
        List[int]
            相邻子空间的索引列表（包括中心子空间本身，如果include_self=True）
        """
        # 将1D索引转换为3D坐标
        i, j, k = self._index_1d_to_3d(subdivision_index)

        neighbors = []

        # 遍历所有相邻的子空间（包括对角相邻）
        for di in [-1, 0, 1]:
            for dj in [-1, 0, 1]:
                for dk in [-1, 0, 1]:
                    # 如果不包含自身，跳过(0,0,0)
                    if not include_self and di == 0 and dj == 0 and dk == 0:
                        continue

                    # 考虑周期性边界条件
                    ni = (i + di) % self.n_subdivisions[0]
                    nj = (j + dj) % self.n_subdivisions[1]
                    nk = (k + dk) % self.n_subdivisions[2]

                    neighbor_index = self._coords_3d_to_1d(
                        np.array([ni]), np.array([nj]), np.array([nk])
                    )[0]
                    neighbors.append(neighbor_index)

        return neighbors

    def get_relevant_atoms_for_subdivision(
        self, subdivision_index: int
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """
        获取计算某个子空间RDF时需要考虑的所有原子

        Parameters
        ----------
        subdivision_index : int
            中心子空间的索引

        Returns
        -------
        Tuple[np.ndarray, Optional[np.ndarray]]
            (atom_indices, mask) - 相关原子的索引数组和对应的布尔掩码
            如果use_masks=False，mask为None
        """
        neighboring_indices = self.get_neighboring_subdivisions(
            subdivision_index, include_self=True
        )

        # 收集所有相邻子空间的原子
        relevant_atom_indices = []

        if self.use_masks:
            n_atoms = len(self.universe.atoms)
            combined_mask = np.zeros(n_atoms, dtype=bool)
        else:
            combined_mask = None

        for neighbor_idx in neighboring_indices:
            if neighbor_idx in self.subdivision_atom_lists:
                atom_indices = self.subdivision_atom_lists[neighbor_idx]
                relevant_atom_indices.extend(atom_indices)
                if self.use_masks:
                    combined_mask[atom_indices] = True

        relevant_atom_indices = np.array(relevant_atom_indices, dtype=int)

        return relevant_atom_indices, combined_mask

    def get_subdivision_for_atom(self, atom_index: int) -> int:
        """
        获取某个原子所属的子空间索引

        Parameters
        ----------
        atom_index : int
            原子索引

        Returns
        -------
        int
            子空间索引
        """
        if atom_index < 0 or atom_index >= len(self.atom_subdivision_indices):
            raise ValueError(f"Invalid atom index: {atom_index}")

        return self.atom_subdivision_indices[atom_index]

    def get_atoms_in_subdivision(self, subdivision_index: int) -> np.ndarray:
        """
        获取某个子空间中的所有原子索引

        Parameters
        ----------
        subdivision_index : int
            子空间索引

        Returns
        -------
        np.ndarray
            该子空间中的原子索引数组
        """
        if subdivision_index in self.subdivision_atom_lists:
            return self.subdivision_atom_lists[subdivision_index].copy()
        else:
            return np.array([], dtype=int)

    def get_subdivision_center(self, subdivision_index: int) -> np.ndarray:
        """
        获取子空间的几何中心位置

        Parameters
        ----------
        subdivision_index : int
            子空间索引

        Returns
        -------
        np.ndarray
            子空间中心的坐标，形状为(3,)
        """
        i, j, k = self._index_1d_to_3d(subdivision_index)

        center = np.array(
            [
                (i + 0.5) * self.subdivision_lengths[0],
                (j + 0.5) * self.subdivision_lengths[1],
                (k + 0.5) * self.subdivision_lengths[2],
            ]
        )

        return center

    def get_subdivision_bounds(
        self, subdivision_index: int
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        获取子空间的边界

        Parameters
        ----------
        subdivision_index : int
            子空间索引

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            (min_bounds, max_bounds) - 最小和最大边界坐标
        """
        i, j, k = self._index_1d_to_3d(subdivision_index)

        min_bounds = np.array(
            [
                i * self.subdivision_lengths[0],
                j * self.subdivision_lengths[1],
                k * self.subdivision_lengths[2],
            ]
        )

        max_bounds = np.array(
            [
                (i + 1) * self.subdivision_lengths[0],
                (j + 1) * self.subdivision_lengths[1],
                (k + 1) * self.subdivision_lengths[2],
            ]
        )

        return min_bounds, max_bounds

    def get_statistics(self) -> Dict:
        """
        获取空间分割的统计信息

        Returns
        -------
        Dict
            包含分割统计信息的字典
        """
        # 计算每个子空间的原子数量
        atom_counts = []
        occupied_subdivisions = 0

        for subdivision_idx in range(self.total_subdivisions):
            if subdivision_idx in self.subdivision_atom_lists:
                count = len(self.subdivision_atom_lists[subdivision_idx])
                atom_counts.append(count)
                occupied_subdivisions += 1
            else:
                atom_counts.append(0)

        atom_counts = np.array(atom_counts)

        stats = {
            "total_subdivisions": self.total_subdivisions,
            "occupied_subdivisions": occupied_subdivisions,
            "empty_subdivisions": self.total_subdivisions - occupied_subdivisions,
            "total_atoms": len(self.universe.atoms),
            "atoms_per_subdivision": {
                "mean": (
                    np.mean(atom_counts[atom_counts > 0])
                    if occupied_subdivisions > 0
                    else 0
                ),
                "std": (
                    np.std(atom_counts[atom_counts > 0])
                    if occupied_subdivisions > 0
                    else 0
                ),
                "min": (
                    np.min(atom_counts[atom_counts > 0])
                    if occupied_subdivisions > 0
                    else 0
                ),
                "max": np.max(atom_counts),
            },
            "subdivision_dimensions": self.n_subdivisions.tolist(),
            "subdivision_lengths": self.subdivision_lengths.tolist(),
            "cutoff_radius": self.r_cutoff,
            "box_lengths": self.box_lengths.tolist(),
        }

        return stats

    def print_statistics(self):
        """打印空间分割的统计信息"""
        stats = self.get_statistics()

        print("\n=== Spatial Subdivision Statistics ===")
        print(f"Box dimensions: {stats['box_lengths']}")
        print(f"Cutoff radius: {stats['cutoff_radius']:.3f}")
        print(f"Subdivision grid: {stats['subdivision_dimensions']}")
        print(
            f"Subdivision lengths: {[f'{x:.3f}' for x in stats['subdivision_lengths']]}"
        )
        print(f"Total subdivisions: {stats['total_subdivisions']}")
        print(f"Occupied subdivisions: {stats['occupied_subdivisions']}")
        print(f"Empty subdivisions: {stats['empty_subdivisions']}")
        print(f"Total atoms: {stats['total_atoms']}")
        print("Atoms per subdivision:")
        print(f"  Mean: {stats['atoms_per_subdivision']['mean']:.1f}")
        print(f"  Std:  {stats['atoms_per_subdivision']['std']:.1f}")
        print(f"  Min:  {stats['atoms_per_subdivision']['min']}")
        print(f"  Max:  {stats['atoms_per_subdivision']['max']}")
        print("=" * 40)
