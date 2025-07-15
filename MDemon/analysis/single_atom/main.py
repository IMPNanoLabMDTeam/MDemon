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

from .base import AnalysisConfig, validate_universe
from .rdf import RDFAnalyzer


class SingleAtomAnalysis:
    """
    单原子分析主接口

    这个类提供了一个统一的接口来访问所有单原子分析功能，包括RDF、扩散、配位数等分析。
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
        # 创建缓存键
        cache_key = (analyzer_type, tuple(sorted(kwargs.items())))

        if cache_key not in self._analyzers:
            # 合并默认参数和用户参数
            merged_kwargs = {**self._default_kwargs, **kwargs}

            # 创建分析器实例
            if analyzer_type == "rdf":
                self._analyzers[cache_key] = RDFAnalyzer(
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

        # 扩散分析（待实现）
        results["diffusion"] = None

        # 配位数分析（待实现）
        results["coordination"] = None

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
            "diffusion": "not_implemented",
            "coordination": "not_implemented",
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
