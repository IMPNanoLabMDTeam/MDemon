import numpy as np

from ..core.structureattr import BondOrder, Composition, Connection, LonePair, Valence
from .base import ReaderBase


def fconnect2bonds(fconnect):
    row = np.zeros(0, dtype=np.int32)
    col = np.zeros(0, dtype=np.int32)
    row_connec = np.zeros(0, dtype=np.int32)
    col_connec = np.zeros(0, dtype=np.int32)
    data_connec = np.zeros(0, dtype=np.int32)
    nbd = 0
    crsmtrx = fconnect._source_register[("Atom_Base", "Atom_Base")].values
    N = crsmtrx.shape[0]
    for i in range(N - 1):
        for j in crsmtrx.indices:
            if j > i:
                row = np.concatenate((row, [nbd, nbd]))
                col = np.concatenate((col, [i, j]))
                row_connec = np.concatenate((row_connec, [i, j]))
                col_connec = np.concatenate((col_connec, [j, i]))
                data_connec = np.concatenate((data_connec, [nbd, nbd]))
                nbd += 1

    data = np.ones(len(col), dtype=np.float32)
    dbase = fconnect._database
    Composition(np.array([row, col, data]), sid=("Bond_Base", "Atom_Base"), N=nbd, M=N)

    Connection(
        np.array([row_connec, col_connec, data_connec]),
        sid=("Atom_Base", "Atom_Base"),
        database=dbase,
        N=N,
        M=N,
    )

    bond_ids = np.arange(1, nbd + 1, dtype=np.int32)
    dbase.id._update_source(bond_ids, ("Bond_Base",))
    dbase.ix._update_source(np.arange(nbd, dtype=np.int32), ("Bond_Base",))


# TODO: Develop a dynamic-reader version.
class REAXFFReader(ReaderBase):
    """ """

    format = "REAXFF"

    def iterdata(self):
        with open(self.filename) as f:
            for line in f:
                if "#" not in line:
                    yield line

    def parse(self, **kwargs):
        """
        Parses a REAXFF_ BOND file.

        Returns
        -------
        BOND
        """

        datalines = list(self.iterdata())

        dbase = kwargs["database"]
        n_atoms = dbase.n_atoms

        atom_ids = np.zeros(n_atoms, dtype=np.int32)
        types = np.zeros(n_atoms, dtype=np.int32)
        charges = np.zeros(n_atoms, dtype=np.float32)
        lonepairs = np.zeros(n_atoms, dtype=np.float32)

        n_bonds = 0

        for i, line in enumerate(datalines):
            line = line.split()
            atom_ids[i] = line[0]
            types[i] = line[1]
            nb = np.int32(line[2])
            charges[i] = line[6 + nb * 2]
            lonepairs[i] = line[5 + nb * 2]

            n_bonds += nb

        order = np.argsort(atom_ids)
        charges = charges[order]
        dbase.charge._update_source(
            charges, ("Atom_Base",), ix=np.arange(n_atoms, dtype=np.int32)
        )

        lonepairs = lonepairs[order]
        LonePair(lonepairs, sid=("Atom_Base",), database=dbase)

        n_bonds /= 2
        n_bonds = np.int32(n_bonds)
        bond_ids = np.arange(1, n_bonds + 1, dtype=np.int32)
        bos = np.zeros(n_bonds, dtype=np.float32)
        row = np.zeros(2 * n_bonds, dtype=np.int32)
        col = np.zeros(2 * n_bonds, dtype=np.int32)
        data = np.ones(2 * n_bonds, dtype=np.int32)
        col_b = np.zeros(2 * n_bonds, dtype=np.int32)
        m = 0
        n = 0
        combi_sets = {}
        for i, line in enumerate(datalines):
            line = line.split()
            nb = np.int32(line[2])
            a1_id = atom_ids[i]
            row[m : m + nb] = np.full(nb, a1_id, dtype=np.int32)
            col[m : m + nb] = line[3 : 3 + nb]
            for j in range(nb):
                a2_id = np.int32(line[3 + j])
                combi = frozenset((a1_id, a2_id))
                if combi not in combi_sets:
                    combi_sets[combi] = n
                    b_id = bond_ids[n]
                    bos[n] = line[4 + nb + j]
                    n += 1
                else:
                    b_id = bond_ids[combi_sets[combi]]
                    combi_sets.pop(combi)

                col_b[m + j] = a2_id
                data[m + j] = b_id - 1
            m += nb

        BondOrder(bos, sid=("Bond_Base",), database=dbase)
        # "Minus 1" scheme is only a temporal plan before cache function implementation,
        # which is essential to construct a reverse table(from id to ix).
        row -= 1
        col -= 1
        col_b -= 1

        dbase.composition._update_source(
            np.array([data, col_b, np.zeros(2 * n_bonds, dtype=np.float32)]),
            ("Bond_Base", "Atom_Base"),
            N=n_bonds,
            M=n_atoms,
            reverse=True,
        )

        Valence.start_from_bondorders(bos, database=dbase)

        Connection(
            np.array([row, col, data]),
            sid=("Atom_Base", "Atom_Base"),
            database=dbase,
            N=n_atoms,
            M=n_atoms,
        )

        dbase.id._update_source(bond_ids, ("Bond_Base",))
        dbase.ix._update_source(np.arange(n_bonds, dtype=np.int32), ("Bond_Base",))
        dbase.connection.union("Atom_Base", "Molecule_Base", process_attr=False)
        return dbase
