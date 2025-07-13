import os

import numpy as np

import MDemon as md


def test_create_universe_from_lammpsdatafile():
    data_dir = os.path.join("tests", "data", "lammps", "CNT")
    filename = "SWNT-8-8-graphene_hole_h2o_random-66666-1-10-30-3.data"
    filepath = os.path.join(data_dir, filename)
    u = md.Universe(filepath)
    a = u.atoms[0]
    assert len(a.coordinate) == 3
    assert list(a.mle_top.keys())[0] >= 0
    a.detailed = False
    assert a.mle_top >= 0
    box = u.box
    assert len(box) == 12
    m = u.molecules[0]
    assert len(m.atms) > 0
    ag = u.angles[0]
    assert len(ag.atms) == 3

    # test distance function
    from MDemon.lib.distance import distance_array, self_distance_array

    d1 = distance_array(u.atoms[0], u.atoms)
    assert d1.shape == tuple([1, len(u.atoms)])
    d2 = self_distance_array(u.atoms[0:4])
    assert len(d2) == 6


def test_create_universe_from_LammpsReaxff():
    data_dir = os.path.join("tests", "data", "lammps", "IrradiatedKapton")
    modelpri = "KAPTON5_504-33r_irradiated"
    datafile = os.path.join(data_dir, modelpri + ".data")
    bondfile = os.path.join(data_dir, modelpri + ".reaxff")

    u = md.Universe(datafile, bondfile)
    b1_ix = min(u.atoms[0].bnds)
    a1_ix = min(u.bonds[b1_ix].atm_top)
    assert a1_ix == 0

    mol = u.molecules[10]
    assert type(mol.mass) is np.float32 or type(mol.mass) is np.float64
    mol.create_rings(multiring=True)
    assert u.rings[0].atms
    assert u.multirings[0].rngs and u.multirings[0].atms


def test_atomic_style_auto_detection():
    """测试atomic格式的自动检测"""
    data_dir = os.path.join("tests", "data", "lammps", "atomic_test")
    filename = "simple_atomic.data"
    filepath = os.path.join(data_dir, filename)

    u = md.Universe(filepath)
    assert len(u.atoms) == 4
    assert len(set(atom.species for atom in u.atoms)) == 2

    # 验证第一个原子的属性
    atom0 = u.atoms[0]
    coord = atom0.coordinate
    assert np.allclose(coord, [0.0, 0.0, 0.0], atol=1e-6)
    assert atom0.species == 1


def test_atomic_style_explicit():
    """测试明确指定atomic格式"""
    data_dir = os.path.join("tests", "data", "lammps", "atomic_test")
    filename = "simple_atomic.data"
    filepath = os.path.join(data_dir, filename)

    u = md.Universe(filepath, atom_style="id type x y z")
    assert len(u.atoms) == 4

    # 验证原子坐标和类型
    expected_coords = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
    ]
    expected_types = [1, 1, 1, 2]

    for i, atom in enumerate(u.atoms):
        assert np.allclose(atom.coordinate, expected_coords[i], atol=1e-6)
        assert atom.species == expected_types[i]


def test_atomic_style_with_flags():
    """测试包含flag信息的atomic格式"""
    data_dir = os.path.join("tests", "data", "lammps", "atomic_test")
    filename = "atomic_with_flags.data"
    filepath = os.path.join(data_dir, filename)

    # 自动检测8字段格式
    u1 = md.Universe(filepath)
    # 明确指定格式
    u2 = md.Universe(filepath, atom_style="id type x y z")

    # 两种方式结果应该一致
    assert len(u1.atoms) == len(u2.atoms) == 4

    for i in range(len(u1.atoms)):
        assert np.allclose(u1.atoms[i].coordinate, u2.atoms[i].coordinate, atol=1e-6)
        assert u1.atoms[i].species == u2.atoms[i].species
