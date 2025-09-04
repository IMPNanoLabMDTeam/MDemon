"""
Main interface for single atom analysis in MDemon

This module provides a unified interface for all single atom analysis capabilities,
making it easy to perform multiple types of analysis on the same system.

Classes:
    SingleAtomAnalysis: Main unified interface for single atom analysis
"""

import warnings
from typing import Any, Dict, List, Optional, Union

import numpy as np

from .angular import AngularAnalyzer
from .base import AnalysisConfig, validate_universe
from .coordination import CoordinationAnalyzer
from .rdf import RDFAnalyzer


class SingleAtomAnalysis:
    """
    单原子分析主接口

    这个类提供了一个统一的接口来访问所有单原子分析功能，包括RDF、角分布、配位数等分析。
    它自动管理不同分析器的实例化和配置，简化用户的使用体验。

    Parameters
    ----------
    universe : MDemon.Universe
        分析的宇宙对象
    config : AnalysisConfig or dict, optional
        分析配置对象或配置字典
    **kwargs : dict
        传递给各个分析器的通用参数

    Attributes
    ----------
    universe : MDemon.Universe
        分析的宇宙对象
    config : AnalysisConfig
        分析配置
    _analyzers : dict
        缓存的分析器实例

    Examples
    --------
    >>> from MDemon.analysis.single_atom import SingleAtomAnalysis
    >>> analysis = SingleAtomAnalysis(universe)
    >>>
    >>> # RDF分析
    >>> rdf_result = analysis.rdf(atom_indices=[0, 1, 2])
    >>>
    >>> # 配位数分析
    >>> cutoff_radii = {'1-1': 2.5, '1-2': 3.0}
    >>> coord_result = analysis.coordination(cutoff_radii=cutoff_radii, atom_indices=[0, 1, 2])
    >>>
    >>> # 角分布分析
    >>> angular_result = analysis.angular(cutoff_radii=cutoff_radii, atom_indices=[0, 1, 2])
    >>>
    >>> # 执行所有分析
    >>> all_results = analysis.analyze_all(cutoff_radii=cutoff_radii, atom_indices=[0, 1, 2])
    """

    def __init__(self, universe, config=None, **kwargs):
        """
        初始化单原子分析主接口

        Parameters
        ----------
        universe : MDemon.Universe
            分析的宇宙对象
        config : AnalysisConfig or dict, optional
            分析配置对象或配置字典
        **kwargs : dict
            传递给各个分析器的通用参数
        """
        # 验证universe对象
        validate_universe(universe)

        self.universe = universe
        self.config = config if config is not None else AnalysisConfig()
        self._analyzers = {}
        self._default_kwargs = kwargs

        # 系统信息
        self.n_atoms = len(universe.atoms)

    def _get_analyzer(self, analyzer_type, **kwargs):
        """
        获取或创建分析器实例

        Parameters
        ----------
        analyzer_type : str
            分析器类型
        **kwargs : dict
            传递给分析器的参数

        Returns
        -------
        analyzer
            分析器实例
        """
        # 创建缓存键，处理字典类型的参数
        hashable_kwargs = []
        for key, value in sorted(kwargs.items()):
            if isinstance(value, dict):
                # 将字典转换为可哈希的tuple of tuples
                hashable_value = tuple(sorted(value.items()))
            elif isinstance(value, list):
                # 将列表转换为tuple
                hashable_value = tuple(value)
            elif isinstance(value, np.ndarray):
                # 将numpy数组转换为tuple
                hashable_value = tuple(value.flatten())
            else:
                hashable_value = value
            hashable_kwargs.append((key, hashable_value))

        cache_key = (analyzer_type, tuple(hashable_kwargs))

        if cache_key not in self._analyzers:
            # 合并默认参数和用户参数
            merged_kwargs = {**self._default_kwargs, **kwargs}

            # 创建分析器实例
            if analyzer_type == "rdf":
                self._analyzers[cache_key] = RDFAnalyzer(
                    universe=self.universe, **merged_kwargs
                )
            elif analyzer_type == "angular":
                self._analyzers[cache_key] = AngularAnalyzer(
                    universe=self.universe, **merged_kwargs
                )
            elif analyzer_type == "coordination":
                self._analyzers[cache_key] = CoordinationAnalyzer(
                    universe=self.universe, **merged_kwargs
                )
            else:
                raise ValueError(f"Unknown analyzer type: {analyzer_type}")

        return self._analyzers[cache_key]

    def rdf(
        self,
        atom_selection=None,
        atom_indices=None,
        r_range=(0.0, 10.0),
        n_bins=100,
        reference_atoms=None,
        scheduler="threads",
        **kwargs,
    ):
        """
        径向分布函数分析

        Parameters
        ----------
        atom_selection : numpy.ndarray, optional
            原子选择掩码，长度为len(universe.atoms)的布尔数组
        atom_indices : List[int], optional
            要分析的原子索引列表，如果提供则覆盖atom_selection
        r_range : tuple, optional
            径向距离范围 (r_min, r_max)，Default: (0.0, 10.0)
        n_bins : int, optional
            距离分箱数量，Default: 100
        reference_atoms : numpy.ndarray, optional
            参考原子的索引数组，如果为None则使用所有原子
        scheduler : str, optional
            Dask调度器类型，Default: 'threads'
        **kwargs : dict
            传递给RDFAnalyzer的其他参数

        Returns
        -------
        RDFResult
            RDF分析结果对象
        """
        # 获取分析器
        analyzer = self._get_analyzer(
            "rdf",
            atom_selection=atom_selection,
            r_range=r_range,
            n_bins=n_bins,
            scheduler=scheduler,
            **kwargs,
        )

        # 执行分析
        return analyzer.analyze_parallel(
            atom_indices=atom_indices, reference_atoms=reference_atoms, **kwargs
        )

    def angular(
        self,
        cutoff_radii,
        atom_selection=None,
        atom_indices=None,
        angle_range=(0.0, np.pi),
        n_bins=90,
        scheduler="threads",
        **kwargs,
    ):
        """
        角分布分析

        Parameters
        ----------
        atom_selection : numpy.ndarray, optional
            原子选择掩码，长度为len(universe.atoms)的布尔数组
        atom_indices : List[int], optional
            要分析的原子索引列表，如果提供则覆盖atom_selection
        cutoff_radii : dict
            三原子类型的截断半径字典，格式为 {'1-1': 2.5, '1-2': 3.0, ...}
            其中键表示三原子类型：
            - '1-1' 表示 1-1-1 角分布（同类型原子）
            - '1-2' 表示 2-1-2 角分布（不同类型原子）
            此参数是必需的，因为截断半径是高度特异的
        angle_range : tuple, optional
            角度范围（弧度） (angle_min, angle_max)，Default: (0.0, π)
        n_bins : int, optional
            角度分箱数量，Default: 90
        scheduler : str, optional
            Dask调度器类型，Default: 'threads'
        **kwargs : dict
            传递给AngularAnalyzer的其他参数

        Returns
        -------
        AngularResult
            角分布分析结果对象

        Notes
        -----
        角分布分析计算中心原子与两个邻近原子之间的夹角分布。
        支持的三原子类型：
        - a-a-a类型：同类型原子组合
        - a-b-a类型：不同类型原子组合

        截断半径用于限制邻近原子的搜索范围，不同原子对可以设置不同的截断半径。
        """
        # 获取分析器
        analyzer = self._get_analyzer(
            "angular",
            atom_selection=atom_selection,
            cutoff_radii=cutoff_radii,
            angle_range=angle_range,
            n_bins=n_bins,
            scheduler=scheduler,
            **kwargs,
        )

        # 执行分析
        return analyzer.analyze_parallel(atom_indices=atom_indices, **kwargs)

    def coordination(
        self,
        cutoff_radii,
        atom_selection=None,
        atom_indices=None,
        scheduler="threads",
        **kwargs,
    ):
        """
        配位数分析

        Parameters
        ----------
        cutoff_radii : dict
            中心-邻近原子对的截断半径字典，格式为 {'1-1': 2.5, '1-2': 3.0, ...}
            其中键表示中心-邻近原子对：
            - '1-1' 表示中心原子类型1，邻近原子类型1
            - '1-2' 表示中心原子类型1，邻近原子类型2
            此参数是必需的，因为截断半径是高度特异的
        atom_selection : numpy.ndarray, optional
            原子选择掩码，长度为len(universe.atoms)的布尔数组
        atom_indices : List[int], optional
            要分析的原子索引列表，如果提供则覆盖atom_selection
        scheduler : str, optional
            Dask调度器类型，Default: 'threads'
        **kwargs : dict
            传递给CoordinationAnalyzer的其他参数

        Returns
        -------
        CoordinationResult
            配位数分析结果对象

        Notes
        -----
        配位数分析计算每个原子的配位环境，与RDF和角分布分析不同，
        配位数分析重点保留每个原子的个体配位数值，而不是平均分布。

        这使得配位数分析特别适用于：
        - 缺陷分析：识别配位数异常的原子
        - 界面研究：分析界面区域的配位环境
        - 结构多样性：研究局部结构的多样性

        截断半径用于定义配位邻近的范围，不同原子对可以设置不同的截断半径。
        """
        # 获取分析器
        analyzer = self._get_analyzer(
            "coordination",
            atom_selection=atom_selection,
            cutoff_radii=cutoff_radii,
            scheduler=scheduler,
            **kwargs,
        )

        # 执行分析
        return analyzer.analyze_parallel(atom_indices=atom_indices, **kwargs)

    def analyze_all(self, atom_selection=None, atom_indices=None, **kwargs):
        """
        执行所有可用的分析

        Parameters
        ----------
        atom_selection : numpy.ndarray, optional
            原子选择掩码
        atom_indices : List[int], optional
            要分析的原子索引列表
        **kwargs : dict
            传递给各个分析器的参数

        Returns
        -------
        dict
            包含所有分析结果的字典
        """
        results = {}

        # RDF分析
        try:
            results["rdf"] = self.rdf(
                atom_selection=atom_selection, atom_indices=atom_indices, **kwargs
            )
        except Exception as e:
            warnings.warn(f"RDF analysis failed: {e}")
            results["rdf"] = None

        # 角分布分析
        try:
            # 角分布分析需要cutoff_radii参数
            if "cutoff_radii" in kwargs:
                results["angular"] = self.angular(
                    atom_selection=atom_selection, atom_indices=atom_indices, **kwargs
                )
            else:
                warnings.warn("Angular analysis requires 'cutoff_radii' parameter")
                results["angular"] = None
        except Exception as e:
            warnings.warn(f"Angular analysis failed: {e}")
            results["angular"] = None

        # 配位数分析
        try:
            # 配位数分析需要cutoff_radii参数
            if "cutoff_radii" in kwargs:
                results["coordination"] = self.coordination(
                    atom_selection=atom_selection, atom_indices=atom_indices, **kwargs
                )
            else:
                warnings.warn("Coordination analysis requires 'cutoff_radii' parameter")
                results["coordination"] = None
        except Exception as e:
            warnings.warn(f"Coordination analysis failed: {e}")
            results["coordination"] = None

        # 扩散分析（待实现）
        results["diffusion"] = None

        return results

    def get_available_analyses(self):
        """
        获取可用的分析方法列表

        Returns
        -------
        dict
            可用分析方法及其状态
        """
        return {
            "rdf": "available",
            "angular": "available",
            "coordination": "available",
            "diffusion": "not_implemented",
        }

    def get_system_info(self):
        """
        获取系统信息

        Returns
        -------
        dict
            系统信息字典
        """
        info = {
            "n_atoms": len(self.universe.atoms),
            "config": (
                self.config.to_dict()
                if hasattr(self.config, "to_dict")
                else str(self.config)
            ),
            "available_analyses": self.get_available_analyses(),
            "current_timestep": self.universe.timestep,
        }

        return info

    def close(self):
        """关闭所有分析器并清理资源"""
        for analyzer in self._analyzers.values():
            if hasattr(analyzer, "close"):
                analyzer.close()
        self._analyzers.clear()
