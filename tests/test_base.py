import MDemon as md


def test_create_universe_from_lammpsdatafile():
    dir_ = r"D:\MDemon_dev\tests\data\lammps\CNT\\"
    filename = r"SWNT-8-8-graphene_hole_h2o_random-66666-1-10-30-3.data"
    u = md.Universe(dir_ + filename)
    a = u.atoms[0]
    assert len(a.coordinate) == 3
    box = u.box
    assert len(box) == 6
