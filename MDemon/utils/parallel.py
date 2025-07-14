"""
Parallel processing utilities for MDemon using Dask

This module provides Dask-based parallel processing capabilities with automatic
memory management and support for different schedulers (threads, processes, distributed).

Classes:
    DaskParallelManager: Main parallel processing manager with automatic memory management
"""

import warnings
from typing import Any, Callable, List, Optional, Union

import dask
import dask.array as da
import numpy as np
import psutil
from dask import compute, delayed
from dask.distributed import Client, as_completed


class DaskParallelManager:
    """
    基于Dask的并行处理管理器，自动内存管理

    这个类提供了一个统一的接口来处理并行计算任务，支持多种调度器类型，
    并自动管理内存使用和任务调度。

    Parameters
    ----------
    scheduler : str, optional
        Dask调度器类型：'threads', 'processes', 'distributed'
        Default: 'threads'
    memory_limit : str or int, optional
        内存限制，'auto'表示自动检测，或指定GB数
        Default: 'auto'
    n_workers : int, optional
        工作进程/线程数量，None表示自动检测
        Default: None
    **kwargs : dict
        传递给Dask客户端的其他参数

    Attributes
    ----------
    scheduler : str
        当前使用的调度器类型
    memory_limit : str or int
        内存限制设置
    n_workers : int
        工作进程/线程数量
    _client : dask.distributed.Client or None
        分布式客户端实例（仅在distributed模式下使用）
    """

    def __init__(
        self, scheduler="threads", memory_limit="auto", n_workers=None, **kwargs
    ):
        """
        初始化Dask并行处理管理器

        Parameters
        ----------
        scheduler : str, optional
            Dask调度器：'threads', 'processes', 'distributed'
        memory_limit : str or int, optional
            内存限制，'auto'表示自动检测，或指定GB数
        n_workers : int, optional
            工作进程/线程数量，None表示自动检测
        **kwargs : dict
            传递给Dask客户端的其他参数
        """
        self.scheduler = scheduler
        self.memory_limit = self._parse_memory_limit(memory_limit)
        self.n_workers = n_workers or self._get_optimal_workers()
        self._client = None

        # 设置Dask配置
        self._configure_dask()

        # 如果使用分布式调度器，创建客户端
        if scheduler == "distributed":
            try:
                from dask.distributed import Client, LocalCluster

                cluster = LocalCluster(
                    n_workers=self.n_workers, memory_limit=self.memory_limit, **kwargs
                )
                self._client = Client(cluster)

            except ImportError:
                warnings.warn(
                    "dask.distributed not available, falling back to 'threads' scheduler",
                    RuntimeWarning,
                )
                self.scheduler = "threads"

    def _parse_memory_limit(self, memory_limit):
        """解析内存限制参数"""
        if memory_limit == "auto":
            # 自动检测可用内存，使用80%作为安全限制
            available_memory = psutil.virtual_memory().available
            return int(available_memory * 0.8)
        elif isinstance(memory_limit, str) and memory_limit.endswith("GB"):
            # 解析GB格式
            return int(float(memory_limit[:-2]) * 1024**3)
        elif isinstance(memory_limit, (int, float)):
            return int(memory_limit)
        else:
            raise ValueError(f"Invalid memory_limit format: {memory_limit}")

    def _get_optimal_workers(self):
        """获取最优的工作进程数量"""
        cpu_count = psutil.cpu_count(logical=False)  # 物理核心数
        if self.scheduler == "threads":
            # 线程调度器可以使用更多工作线程
            return min(cpu_count * 2, 16)
        else:
            # 进程调度器使用物理核心数
            return cpu_count

    def _configure_dask(self):
        """配置Dask全局设置"""
        # 设置Dask配置以优化性能
        dask.config.set(
            {
                "array.chunk-size": "128MB",
                "array.slicing.split_large_chunks": True,
                "optimization.inline-functions": True,
                "optimization.fuse.active": True,
            }
        )

    def parallel_apply(
        self, func: Callable, data_list: List[Any], chunk_size: Union[str, int] = "auto"
    ) -> List[Any]:
        """
        并行应用函数到数据列表，自动内存管理

        Parameters
        ----------
        func : callable
            要应用的函数
        data_list : list
            输入数据列表
        chunk_size : str or int, optional
            分块大小，'auto'表示自动确定

        Returns
        -------
        list
            计算结果列表
        """
        # 创建delayed任务
        tasks = [delayed(func)(item) for item in data_list]

        # Dask自动处理内存管理和任务调度
        return compute(*tasks, scheduler=self.scheduler)

    def parallel_map_blocks(
        self,
        func: Callable,
        data_array: np.ndarray,
        chunks: Union[str, int] = "auto",
        **kwargs,
    ) -> da.Array:
        """
        对大数组并行应用函数，支持大于内存的数据

        Parameters
        ----------
        func : callable
            要应用的函数
        data_array : numpy.ndarray or dask.array.Array
            输入数组
        chunks : str or int, optional
            分块大小，'auto'表示自动确定
        **kwargs : dict
            传递给map_blocks的其他参数

        Returns
        -------
        dask.array.Array
            计算结果数组
        """
        if not isinstance(data_array, da.Array):
            # 转换为dask数组，自动分块
            data_array = da.from_array(data_array, chunks=chunks)

        return data_array.map_blocks(func, **kwargs)

    def create_dask_array(
        self,
        shape: tuple,
        dtype: np.dtype = np.float64,
        chunks: Union[str, int] = "auto",
    ) -> da.Array:
        """
        创建dask数组，自动内存优化

        Parameters
        ----------
        shape : tuple
            数组形状
        dtype : numpy.dtype, optional
            数组数据类型
        chunks : str or int, optional
            分块大小，'auto'表示自动确定

        Returns
        -------
        dask.array.Array
            创建的dask数组
        """
        return da.zeros(shape, dtype=dtype, chunks=chunks)

    def batch_compute(
        self, tasks: List[Any], batch_size: Union[str, int] = "auto"
    ) -> List[Any]:
        """
        批量计算任务，防止内存溢出

        Parameters
        ----------
        tasks : list
            delayed任务列表
        batch_size : str or int, optional
            批处理大小，'auto'表示自动确定最优批大小

        Returns
        -------
        list
            计算结果列表
        """
        if batch_size == "auto":
            # Dask自动确定最优批大小
            return compute(*tasks, scheduler=self.scheduler)
        else:
            # 手动批处理
            results = []
            for i in range(0, len(tasks), batch_size):
                batch = tasks[i : i + batch_size]
                batch_results = compute(*batch, scheduler=self.scheduler)
                results.extend(batch_results)
            return results

    def map_partitions(
        self, func: Callable, data: da.Array, *args, **kwargs
    ) -> da.Array:
        """
        对数组分区并行应用函数

        Parameters
        ----------
        func : callable
            要应用的函数
        data : dask.array.Array
            输入数组
        *args : tuple
            传递给函数的位置参数
        **kwargs : dict
            传递给函数的关键字参数

        Returns
        -------
        dask.array.Array
            计算结果数组
        """
        return data.map_partitions(func, *args, **kwargs)

    def persist(self, *args, **kwargs):
        """
        持久化数据到内存中，用于重复访问

        Parameters
        ----------
        *args : tuple
            要持久化的对象
        **kwargs : dict
            传递给persist的其他参数

        Returns
        -------
        tuple or single object
            持久化后的对象
        """
        return dask.persist(*args, **kwargs)

    def visualize(self, *args, filename=None, **kwargs):
        """
        可视化计算图

        Parameters
        ----------
        *args : tuple
            要可视化的对象
        filename : str, optional
            输出文件名
        **kwargs : dict
            传递给visualize的其他参数
        """
        return dask.visualize(*args, filename=filename, **kwargs)

    def get_memory_usage(self):
        """
        获取当前内存使用情况

        Returns
        -------
        dict
            内存使用信息
        """
        memory_info = psutil.virtual_memory()
        return {
            "total": memory_info.total,
            "available": memory_info.available,
            "used": memory_info.used,
            "percent": memory_info.percent,
            "limit": self.memory_limit,
        }

    def get_system_info(self):
        """
        获取系统信息

        Returns
        -------
        dict
            系统信息
        """
        return {
            "cpu_count": psutil.cpu_count(logical=False),
            "logical_cpu_count": psutil.cpu_count(logical=True),
            "memory_total": psutil.virtual_memory().total,
            "scheduler": self.scheduler,
            "n_workers": self.n_workers,
            "dask_version": dask.__version__,
        }

    def close(self):
        """关闭分布式客户端，释放资源"""
        if self._client:
            self._client.close()
            self._client = None

    def __enter__(self):
        """上下文管理器支持"""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """上下文管理器支持"""
        self.close()

    def __del__(self):
        """析构函数，确保资源清理"""
        self.close()


# 便利函数，提供简化的接口
def parallel_apply(
    func: Callable, data_list: List[Any], scheduler: str = "threads", **kwargs
) -> List[Any]:
    """
    便利函数：并行应用函数到数据列表

    Parameters
    ----------
    func : callable
        要应用的函数
    data_list : list
        输入数据列表
    scheduler : str, optional
        调度器类型：'threads', 'processes', 'distributed'
    **kwargs : dict
        传递给DaskParallelManager的其他参数

    Returns
    -------
    list
        计算结果列表
    """
    with DaskParallelManager(scheduler=scheduler, **kwargs) as manager:
        return manager.parallel_apply(func, data_list)


def create_delayed_tasks(func: Callable, data_list: List[Any]) -> List[Any]:
    """
    便利函数：创建延迟任务列表

    Parameters
    ----------
    func : callable
        要应用的函数
    data_list : list
        输入数据列表

    Returns
    -------
    list
        延迟任务列表
    """
    return [delayed(func)(item) for item in data_list]


def compute_tasks(tasks: List[Any], scheduler: str = "threads", **kwargs) -> List[Any]:
    """
    便利函数：计算延迟任务

    Parameters
    ----------
    tasks : list
        延迟任务列表
    scheduler : str, optional
        调度器类型：'threads', 'processes', 'distributed'
    **kwargs : dict
        传递给compute的其他参数

    Returns
    -------
    list
        计算结果列表
    """
    return compute(*tasks, scheduler=scheduler, **kwargs)
