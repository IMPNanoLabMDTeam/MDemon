import numpy as np

from .base import DynamicReaderBase, ReaderBase, squash_by


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

        dbase = kwargs["database"]
        n_atoms = dbase.n_atoms

        atom_ids = np.zeros(n_atoms, dtype=np.int32)
        types = np.zeros(n_atoms, dtype=np.int32)
        nbonds = np.zeros(n_atoms, dtype=np.int32)
        charges = np.zeros(n_atoms, dtype=np.float32)

        datalines = list(self.iterdata())

        for i, line in enumerate(datalines):
            line = line.split()
            atom_ids[i] = line[0]
            types[i] = line[1]
            nbonds[i] = line[2]
            charges[i] = line[6 + nbonds[i] * 2]
