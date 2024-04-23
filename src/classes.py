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


class Contact:
    def __init__(self, idchain1, residueatom1, idchain2, residueatom2, distance, type, atom1, atom2):
        self.idchain1 = idchain1
        self.residueatom1 = residueatom1
        self.idchain2 = idchain2
        self.residueatom2 = residueatom2
        self.distance = distance
        self.type = type
        self.atom1 = atom1
        self.atom2 = atom2

    def to_list(self): # all info besides the atom objects
        return [self.idchain1, self.residueatom1, self.idchain2, self.residueatom2,
            self.distance, self.type]
        

class Match:
    def __init__(self, avd, contact1, contact2):
        self.avd = avd
        self.contact1 = contact1
        self.contact2 = contact2
        
    def __eq__(self, other):
        if isinstance(other, Match):
            return (self.avd == other.avd and
                    self.contact1 == other.contact1 and
                    self.contact2 == other.contact2)
        return False