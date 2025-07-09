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
