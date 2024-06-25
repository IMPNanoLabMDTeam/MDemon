import MDemon as md


def test_create_universe_from_lammpsdatafile():
    dir_ = r"C:\Users\amphi\MDemon_DatabaseUpdate\tests\data\lammps\CNT\\"
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
