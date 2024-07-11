import numpy as np

import MDemon as md


def test_create_universe_from_lammpsdatafile():
    dir_ = r".\tests\data\lammps\CNT\\"
    filename = r"SWNT-8-8-graphene_hole_h2o_random-66666-1-10-30-3.data"
    u = md.Universe(dir_ + filename)
    a = u.atoms[0]
    assert len(a.coordinate) == 3
    assert list(a.mle_ix.keys())[0] >= 0
    a.detailed = False
    assert a.mle_ix >= 0
    box = u.box
    assert len(box) == 6
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
    dir_ = r".\tests\data\lammps\IrradiatedKapton\\"
    modelpri = r"KAPTON5_504-33r_irradiated"
    datafile = dir_ + modelpri + ".data"
    bondfile = dir_ + modelpri + ".reaxff"

    u = md.Universe(datafile, bondfile)
    b1_ix = min(u.atoms[0].bnd_ix)
    a1_ix = min(u.bonds[b1_ix].atms)
    assert a1_ix == 0

    mol = u.molecules[10]
    assert type(mol.mass) is np.float32
    mol.create_rings(multiring=False)
    mol.draw_atoms("graph")
