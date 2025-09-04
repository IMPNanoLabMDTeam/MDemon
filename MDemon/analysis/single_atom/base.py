"""
Base classes for single atom analysis in MDemon

This module provides the foundation for single atom analysis, including abstract
base classes, parallel processing integration, and atom selection utilities.

Classes:
    SingleAtomAnalyzer: Abstract base class for single atom analysis
    AnalysisConfig: Configuration management for analysis parameters
    AnalysisResult: Base class for analysis results
"""

import warnings
from abc import ABC, abstractmethod
from datetime import datetime
from typing import Any, Dict, List, Optional, Union

import numpy as np
from dask import delayed

from ...utils.parallel import DaskParallelManager
from ...utils.spatial import SpatialSubdivision


class SingleAtomAnalyzer(ABC):
    """
    单原子分析器基类，基于Dask的简化设计

    这个抽象基类为所有单原子分析提供了统一的接口和并行处理能力。
    所有具体的分析器（如RDF、扩散、配位数分析）都应该继承此类。

    Parameters
    ----------
    universe : MDemon.Universe
        分析的宇宙对象，包含原子信息和轨迹数据
    atom_selection : numpy.ndarray, optional
        原子选择掩码，长度为len(universe.atoms)的布尔数组
        默认为None（等价于全选所有原子）
    scheduler : str, optional
        Dask调度器：'threads', 'processes', 'distributed'
        Default: 'threads'
    **kwargs : dict
        传递给DaskParallelManager的其他参数

    Attributes
    ----------
    universe : MDemon.Universe
        分析的宇宙对象
    atom_selection : numpy.ndarray or None
        原子选择掩码
    parallel_manager : DaskParallelManager
        并行处理管理器
    n_selected_atoms : int
        选中的原子数量
    """

    def __init__(
        self,
        universe,
        atom_selection=None,
        scheduler="threads",
        enable_spatial_subdivision=True,
        r_cutoff=None,
        use_subdivision_masks=False,
        show_running_time=False,
        **kwargs,
    ):
        """
        初始化单原子分析器

        Parameters
        ----------
        universe : MDemon.Universe
            分析的宇宙对象
        atom_selection : numpy.ndarray, optional
            原子选择掩码，长度为len(universe.atoms)的布尔数组
            默认为None（等价于全选所有原子）
        scheduler : str, optional
            Dask调度器：'threads', 'processes', 'distributed'
        enable_spatial_subdivision : bool, optional
            是否启用空间预分割优化，默认为True
        r_cutoff : float, optional
            空间预分割的截断半径，如果为None则由子类确定
        use_subdivision_masks : bool, optional
            是否在空间预分割中使用掩码，默认为False。对于大体系建议保持False
        show_running_time : bool, optional
            是否显示运行时间调试信息，默认为False
        **kwargs : dict
            传递给DaskParallelManager的其他参数
        """
        self.universe = universe
        self.atom_selection = self._validate_atom_selection(atom_selection)

        # 直接使用DaskParallelManager，自动内存管理
        self.parallel_manager = DaskParallelManager(scheduler=scheduler, **kwargs)

        # 缓存选中的原子数量
        self.n_selected_atoms = len(self._get_selected_atoms())

        # 分析配置
        self.config = AnalysisConfig()

        # 空间预分割设置
        self.enable_spatial_subdivision = enable_spatial_subdivision
        self.r_cutoff = r_cutoff
        self.use_subdivision_masks = use_subdivision_masks
        self.spatial_subdivision = None

        # 时间调试设置
        self.show_running_time = show_running_time

        # 如果启用空间预分割且有截断半径，立即初始化
        if self.enable_spatial_subdivision and self.r_cutoff is not None:
            self._initialize_spatial_subdivision()

        # 结果缓存
        self._results_cache = {}

    def _validate_atom_selection(self, atom_selection):
        """
        验证原子选择掩码

        Parameters
        ----------
        atom_selection : numpy.ndarray or None
            原子选择掩码

        Returns
        -------
        numpy.ndarray or None
            验证后的原子选择掩码

        Raises
        ------
        ValueError
            如果掩码长度不匹配原子数量
        TypeError
            如果掩码类型不正确
        """
        if atom_selection is None:
            return None

        if not isinstance(atom_selection, np.ndarray):
            try:
                atom_selection = np.array(atom_selection)
            except Exception as e:
                raise TypeError(
                    f"atom_selection must be numpy array or convertible to array: {e}"
                )

        if atom_selection.dtype != bool:
            if np.issubdtype(atom_selection.dtype, np.integer):
                # 如果是整数数组，转换为布尔掩码
                mask = np.zeros(len(self.universe.atoms), dtype=bool)
                mask[atom_selection] = True
                atom_selection = mask
            else:
                raise TypeError(
                    "atom_selection must be boolean array or integer indices"
                )

        if len(atom_selection) != len(self.universe.atoms):
            raise ValueError(
                f"atom_selection mask length ({len(atom_selection)}) "
                f"must match number of atoms ({len(self.universe.atoms)})"
            )

        return atom_selection

    def _get_selected_atoms(self):
        """
        获取选中的原子索引

        Returns
        -------
        list
            选中的原子索引列表
        """
        if self.atom_selection is None:
            # 默认全选，返回所有原子索引
            return list(range(len(self.universe.atoms)))
        else:
            # atom_selection是一个掩码，长度为len(self.universe.atoms)
            return np.where(self.atom_selection)[0].tolist()

    def _initialize_spatial_subdivision(self):
        """初始化空间预分割"""
        if self.r_cutoff is None:
            warnings.warn("Cannot initialize spatial subdivision without r_cutoff")
            return

        try:
            self.spatial_subdivision = SpatialSubdivision(
                self.universe, self.r_cutoff, use_masks=self.use_subdivision_masks
            )
            mask_status = (
                "with masks" if self.use_subdivision_masks else "without masks"
            )
            print(
                f"Spatial subdivision initialized with cutoff {self.r_cutoff:.3f} Å ({mask_status})"
            )
        except Exception as e:
            warnings.warn(f"Failed to initialize spatial subdivision: {e}")
            self.enable_spatial_subdivision = False

    def _get_relevant_atoms_for_analysis(self, atom_index: int):
        """
        获取分析某个原子时需要考虑的所有相关原子

        Parameters
        ----------
        atom_index : int
            中心原子索引

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            (relevant_atom_indices, combined_mask) - 相关原子索引和掩码
            如果未启用空间预分割，返回所有原子
        """
        if self.enable_spatial_subdivision and self.spatial_subdivision is not None:
            # 使用空间预分割获取相关原子
            atom_subdivision = self.spatial_subdivision.get_subdivision_for_atom(
                atom_index
            )
            relevant_atoms, mask = (
                self.spatial_subdivision.get_relevant_atoms_for_subdivision(
                    atom_subdivision
                )
            )

            return relevant_atoms, mask
        else:
            # 如果未启用空间预分割，返回所有原子
            n_atoms = len(self.universe.atoms)
            all_atoms = np.arange(n_atoms)
            all_mask = np.ones(n_atoms, dtype=bool)
            return all_atoms, all_mask

    def set_spatial_subdivision_cutoff(self, r_cutoff: float):
        """
        设置空间预分割的截断半径并重新初始化

        Parameters
        ----------
        r_cutoff : float
            新的截断半径
        """
        self.r_cutoff = r_cutoff
        if self.enable_spatial_subdivision:
            self._initialize_spatial_subdivision()

    @abstractmethod
    def analyze_single_atom(self, atom_index: int, **kwargs) -> Any:
        """
        分析单个原子的抽象方法

        所有子类必须实现此方法来定义具体的分析逻辑。

        Parameters
        ----------
        atom_index : int
            要分析的原子索引
        **kwargs : dict
            分析参数

        Returns
        -------
        Any
            分析结果，具体类型由子类定义
        """
        pass

    def analyze_parallel(self, atom_indices=None, **kwargs):
        """
        并行分析多个原子，Dask自动处理内存管理

        Parameters
        ----------
        atom_indices : list, optional
            要分析的原子索引列表，None表示使用所有选中的原子
        **kwargs : dict
            传递给analyze_single_atom的参数

        Returns
        -------
        dict
            分析结果字典，键为原子索引，值为分析结果
        """
        if atom_indices is None:
            atom_indices = self._get_selected_atoms()
        elif not isinstance(atom_indices, (list, tuple, np.ndarray)):
            atom_indices = [atom_indices]

        # 验证原子索引
        atom_indices = self._validate_atom_indices(atom_indices)

        # 如果启用时间调试，限制原子数量不超过3
        if self.show_running_time and len(atom_indices) > 3:
            atom_indices = atom_indices[:3]
            print(
                f"[DEBUG] show_running_time enabled, limiting analysis to first 3 atoms: {atom_indices}"
            )

        # Dask自动优化内存使用和任务调度
        results = self.parallel_manager.parallel_apply(
            lambda idx: self.analyze_single_atom(idx, **kwargs), atom_indices
        )

        # 将结果组织为字典
        result_dict = dict(zip(atom_indices, results))

        # 创建结果对象
        return self._create_analysis_result(result_dict, **kwargs)

    def analyze_large_system(self, batch_size="auto", **kwargs):
        """
        分析大规模系统，Dask自动分块和内存管理

        Parameters
        ----------
        batch_size : str or int, optional
            批处理大小，'auto'表示自动优化
        **kwargs : dict
            传递给analyze_single_atom的参数

        Returns
        -------
        AnalysisResult
            分析结果对象
        """
        atom_indices = self._get_selected_atoms()

        # 对于大规模系统，使用batch_compute自动优化
        tasks = [
            delayed(self.analyze_single_atom)(idx, **kwargs) for idx in atom_indices
        ]

        results = self.parallel_manager.batch_compute(tasks, batch_size=batch_size)

        # 将结果组织为字典
        result_dict = dict(zip(atom_indices, results))

        # 创建结果对象
        return self._create_analysis_result(result_dict, **kwargs)

    def _validate_atom_indices(self, atom_indices):
        """
        验证原子索引的有效性

        Parameters
        ----------
        atom_indices : list
            原子索引列表

        Returns
        -------
        list
            验证后的原子索引列表

        Raises
        ------
        ValueError
            如果原子索引超出范围
        """
        atom_indices = np.asarray(atom_indices)

        if np.any(atom_indices < 0) or np.any(atom_indices >= len(self.universe.atoms)):
            raise ValueError(
                f"Atom indices must be in range [0, {len(self.universe.atoms)-1}]"
            )

        return atom_indices.tolist()

    def _create_analysis_result(self, result_dict, **kwargs):
        """
        创建分析结果对象

        Parameters
        ----------
        result_dict : dict
            分析结果字典
        **kwargs : dict
            额外的元数据

        Returns
        -------
        AnalysisResult
            分析结果对象
        """
        analysis_type = self.__class__.__name__.replace("Analyzer", "").lower()

        metadata = {
            "analysis_type": analysis_type,
            "scheduler": self.parallel_manager.scheduler,
            "n_atoms": len(result_dict),
            "universe_info": {
                "n_atoms": len(self.universe.atoms),
                "n_frames": getattr(self.universe, "n_frames", 1),
            },
        }
        metadata.update(kwargs)

        return AnalysisResult(
            universe=self.universe,
            analysis_type=analysis_type,
            atom_indices=list(result_dict.keys()),
            data=result_dict,
            metadata=metadata,
        )

    def get_memory_usage(self):
        """
        获取内存使用情况

        Returns
        -------
        dict
            内存使用信息
        """
        return self.parallel_manager.get_memory_usage()

    def get_system_info(self):
        """
        获取系统信息

        Returns
        -------
        dict
            系统信息
        """
        return self.parallel_manager.get_system_info()

    def set_config(self, **kwargs):
        """
        设置分析配置

        Parameters
        ----------
        **kwargs : dict
            配置参数
        """
        self.config.update(**kwargs)

    def close(self):
        """清理资源"""
        if hasattr(self, "parallel_manager"):
            self.parallel_manager.close()

    def __enter__(self):
        """上下文管理器支持"""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """上下文管理器支持"""
        self.close()

    def __del__(self):
        """析构函数，确保资源清理"""
        self.close()


class AnalysisConfig:
    """
    分析配置管理

    这个类管理分析过程中的各种配置参数，包括并行处理、内存管理等。

    Attributes
    ----------
    config : dict
        配置参数字典
    """

    def __init__(self, config_dict=None):
        """
        初始化配置管理器

        Parameters
        ----------
        config_dict : dict, optional
            初始配置字典
        """
        self.config = self._load_default_config()
        if config_dict:
            self.config.update(config_dict)

    def _load_default_config(self):
        """加载默认配置"""
        return {
            "parallel": {
                "scheduler": "threads",
                "n_workers": None,
                "memory_limit": "auto",
            },
            "memory": {
                "chunk_size": "auto",
                "batch_size": "auto",
                "cache_size": "1GB",
            },
            "analysis": {
                "precision": "float64",
                "error_handling": "warn",
                "progress_bar": True,
            },
        }

    def get_parallel_config(self):
        """获取并行配置"""
        return self.config.get("parallel", {})

    def get_memory_config(self):
        """获取内存配置"""
        return self.config.get("memory", {})

    def get_analysis_config(self):
        """获取分析配置"""
        return self.config.get("analysis", {})

    def update(self, **kwargs):
        """更新配置"""
        for key, value in kwargs.items():
            if key in self.config:
                if isinstance(self.config[key], dict) and isinstance(value, dict):
                    self.config[key].update(value)
                else:
                    self.config[key] = value
            else:
                self.config[key] = value

    def optimize_for_system(self, n_atoms, n_frames=1, system_memory=None):
        """
        根据系统参数优化配置

        Parameters
        ----------
        n_atoms : int
            原子数量
        n_frames : int, optional
            帧数量
        system_memory : int, optional
            系统内存（字节），None表示自动检测
        """
        import psutil

        if system_memory is None:
            system_memory = psutil.virtual_memory().available

        # 估算内存需求
        estimated_memory_per_atom = 1024  # 每个原子约1KB
        total_memory_needed = n_atoms * n_frames * estimated_memory_per_atom

        # 优化批处理大小
        if total_memory_needed > system_memory * 0.8:
            # 如果内存不足，使用批处理
            batch_size = max(
                1, int(system_memory * 0.8 / (n_frames * estimated_memory_per_atom))
            )
            self.config["memory"]["batch_size"] = batch_size

        # 优化调度器选择
        if n_atoms > 10000:
            self.config["parallel"]["scheduler"] = "processes"
        elif n_atoms > 1000:
            self.config["parallel"]["scheduler"] = "threads"


class AnalysisResult:
    """
    分析结果基类

    这个类封装了分析结果，提供统一的接口来访问、保存和可视化结果。

    Parameters
    ----------
    analysis_type : str
        分析类型
    atom_indices : list
        分析的原子索引列表
    data : dict
        分析结果数据
    metadata : dict, optional
        元数据信息

    Attributes
    ----------
    analysis_type : str
        分析类型
    atom_indices : list
        分析的原子索引列表
    data : dict
        分析结果数据
    metadata : dict
        元数据信息
    timestamp : datetime
        创建时间戳
    """

    def __init__(
        self,
        analysis_type,
        atom_indices,
        data,
        metadata=None,
        universe=None,
        averaged_data=None,
        distributions=None,
    ):
        """
        初始化分析结果

        Parameters
        ----------
        universe : Universe
            宇宙对象
        analysis_type : str
            分析类型
        atom_indices : list
            分析的原子索引列表
        data : dict
            分析结果数据
        metadata : dict, optional
            元数据信息
        """
        self.universe = universe
        self.analysis_type = analysis_type
        self.atom_indices = atom_indices
        self.data = data
        self.metadata = metadata or {}
        self.timestamp = datetime.now()
        self.averaged_data = averaged_data
        self.distributions = distributions

    def get_result(self, atom_index):
        """
        获取特定原子的结果

        Parameters
        ----------
        atom_index : int
            原子索引

        Returns
        -------
        Any
            该原子的分析结果
        """
        return self.data.get(atom_index)

    def get_all_results(self):
        """
        获取所有结果

        Returns
        -------
        dict
            所有分析结果
        """
        return self.data

    def to_dict(self):
        """
        转换为字典格式

        Returns
        -------
        dict
            字典格式的结果
        """
        return {
            "analysis_type": self.analysis_type,
            "atom_indices": self.atom_indices,
            "data": self.data,
            "metadata": self.metadata,
            "timestamp": self.timestamp.isoformat(),
            "averaged_data": self.averaged_data,
            "distributions": self.distributions,
        }

    def save(self, filename, format="pickle"):
        """
        保存结果到文件

        Parameters
        ----------
        filename : str
            文件名
        format : str, optional
            文件格式：'pickle', 'json', 'hdf5'
        """
        if format == "pickle":
            import pickle

            with open(filename, "wb") as f:
                pickle.dump(self.to_dict(), f)
        elif format == "json":
            import json

            with open(filename, "w") as f:
                json.dump(self.to_dict(), f, indent=2, default=str)
        elif format == "hdf5":
            try:
                import h5py

                with h5py.File(filename, "w") as f:
                    f.attrs["analysis_type"] = self.analysis_type
                    f.attrs["timestamp"] = self.timestamp.isoformat()

                    # 保存原子索引
                    f.create_dataset("atom_indices", data=self.atom_indices)

                    # 保存数据
                    data_group = f.create_group("data")
                    for atom_idx, result in self.data.items():
                        data_group.create_dataset(str(atom_idx), data=result)

            except ImportError:
                warnings.warn("h5py not available, falling back to pickle format")
                self.save(filename, format="pickle")
        else:
            raise ValueError(f"Unsupported format: {format}")

    @classmethod
    def load(cls, filename, format="pickle"):
        """
        从文件加载结果

        Parameters
        ----------
        filename : str
            文件名
        format : str, optional
            文件格式：'pickle', 'json', 'hdf5'

        Returns
        -------
        AnalysisResult
            加载的结果对象
        """
        if format == "pickle":
            import pickle

            with open(filename, "rb") as f:
                data = pickle.load(f)
        elif format == "json":
            import json

            with open(filename, "r") as f:
                data = json.load(f)
        elif format == "hdf5":
            try:
                import h5py

                with h5py.File(filename, "r") as f:
                    analysis_type = f.attrs["analysis_type"]
                    atom_indices = f["atom_indices"][:].tolist()

                    # 加载数据
                    data_dict = {}
                    for atom_idx in f["data"].keys():
                        data_dict[int(atom_idx)] = f["data"][atom_idx][:]

                    data = {
                        "analysis_type": analysis_type,
                        "atom_indices": atom_indices,
                        "data": data_dict,
                        "metadata": {},
                        "timestamp": f.attrs["timestamp"],
                    }
            except ImportError:
                warnings.warn("h5py not available, falling back to pickle format")
                return cls.load(filename, format="pickle")
        else:
            raise ValueError(f"Unsupported format: {format}")

        return cls(
            analysis_type=data["analysis_type"],
            atom_indices=data["atom_indices"],
            data=data["data"],
            metadata=data["metadata"],
            averaged_data=data["averaged_data"],
            distributions=data["distributions"] if "distributions" in data else None,
        )

    def __repr__(self):
        return (
            f"AnalysisResult(analysis_type='{self.analysis_type}', "
            f"n_atoms={len(self.atom_indices)}, "
            f"timestamp='{self.timestamp.strftime('%Y-%m-%d %H:%M:%S')}')"
        )

    def __len__(self):
        return len(self.atom_indices)

    def __getitem__(self, atom_index):
        return self.get_result(atom_index)


# 辅助函数
def create_atom_selection_mask(universe, **criteria):
    """
    创建原子选择掩码的便利函数

    Parameters
    ----------
    universe : MDemon.Universe
        宇宙对象
    **criteria : dict
        选择条件，如 species=1, element=6, charge_range=(0, 1)

    Returns
    -------
    numpy.ndarray
        布尔掩码数组
    """
    mask = np.ones(len(universe.atoms), dtype=bool)

    for key, value in criteria.items():
        if key == "species":
            atom_mask = np.array([atom.species == value for atom in universe.atoms])
        elif key == "element":
            atom_mask = np.array([atom.element == value for atom in universe.atoms])
        elif key == "charge_range":
            min_charge, max_charge = value
            atom_mask = np.array(
                [min_charge <= atom.charge <= max_charge for atom in universe.atoms]
            )
        elif key == "mass_range":
            min_mass, max_mass = value
            atom_mask = np.array(
                [min_mass <= atom.mass <= max_mass for atom in universe.atoms]
            )
        elif key == "position_range":
            axis, min_pos, max_pos = value
            atom_mask = np.array(
                [min_pos <= atom.coordinate[axis] <= max_pos for atom in universe.atoms]
            )
        elif key == "indices":
            atom_mask = np.zeros(len(universe.atoms), dtype=bool)
            atom_mask[value] = True
        else:
            warnings.warn(f"Unknown selection criterion: {key}")
            continue

        mask = mask & atom_mask

    return mask


def validate_universe(universe):
    """
    Validate Universe objects for single atom analysis

    This function checks if a Universe object has the required attributes
    and structure for single atom analysis.

    Parameters
    ----------
    universe : MDemon.Universe
        The Universe object to validate

    Returns
    -------
    bool
        True if the universe is valid for single atom analysis

    Raises
    ------
    ValueError
        If the universe is not suitable for single atom analysis
    """
    if universe is None:
        raise ValueError("Universe cannot be None")

    # Check if universe has atoms
    if not hasattr(universe, "atoms"):
        raise ValueError("Universe must have atoms attribute")

    if len(universe.atoms) == 0:
        raise ValueError("Universe contains no atoms")

    # Check if atoms have required attributes
    for i, atom in enumerate(universe.atoms):
        if not hasattr(atom, "coordinate"):
            raise ValueError(f"Atom {i} is missing coordinate attribute")

        if not hasattr(atom, "species"):
            warnings.warn(f"Atom {i} is missing species attribute")

        if not hasattr(atom, "element"):
            warnings.warn(f"Atom {i} is missing element attribute")

    # Check if universe has n_atoms attribute
    if hasattr(universe, "n_atoms"):
        if universe.n_atoms != len(universe.atoms):
            warnings.warn(
                f"Universe.n_atoms ({universe.n_atoms}) does not match "
                f"actual number of atoms ({len(universe.atoms)})"
            )

    return True
