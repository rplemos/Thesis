class Protein:
    def __init__(self):
        self.title = None
        self.id = None
        self.chains = []

    def set_title(self, title):
        if self.title is None:
            self.title = title
        else:
            self.title += " " + title.strip()

    def get_chains(self):
        for chain in self.chains:
            yield chain

    def get_residues(self):
        for chain in self.get_chains():
            for residue in chain.residues:
                yield residue

    def get_atoms(self):
        for residue in self.get_residues():
            for atom in residue.atoms:
                yield atom


class Chain:
    def __init__(self):
        self.id = None
        self.residues = []

    def get_residues(self):
        for residue in self.residues:
            yield residue


class Residue:
    def __init__(self):
        self.resnum = None
        self.resname = None
        self.atoms = []
        self.chain = None
        self.ring = False
        self.normal_vector = None

    def get_atoms(self):
        for atom in self.atoms:
            yield atom


class Atom:
    def __init__(self):
        self.atomname = None
        self.x = None
        self.y = None
        self.z = None
        self.residue = None

    def set_atom_info(self, atomname, x, y, z):
        self.atomname = atomname
        self.x = x
        self.y = y
        self.z = z
