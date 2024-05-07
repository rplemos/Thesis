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
                
    def count_residues(self):
        return sum(1 for _ in self.get_residues())


class Chain:
    def __init__(self, id, residues):
        self.id = id
        self.residues = residues

    def get_residues(self):
        for residue in self.residues:
            yield residue


class Residue:
    def __init__(self, resnum, resname, atoms, chain, ring, normal_vector):
        self.resnum = resnum
        self.resname = resname
        self.atoms = atoms
        self.chain = chain
        self.ring = ring
        self.normal_vector = normal_vector

    def get_atoms(self):
        for atom in self.atoms:
            yield atom
        
class Atom:
    def __init__(self, atomname, x, y, z, occupancy, residue):
        self.atomname = atomname
        self.x = x
        self.y = y
        self.z = z
        self.occupancy = occupancy
        self.residue = residue


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
    
    def to_contact(self):
        return [self.idchain1.split(":")[1], self.residueatom1.split(":")[0], self.idchain2.split(":")[1], 
                self.residueatom2.split(":")[0], self.type]
        

class Match:
    def __init__(self, avd, contact1, contact2, d3d4):
        self.avd = avd
        self.contact1 = contact1
        self.contact2 = contact2
        self.d3d4 = d3d4 
        
    def __eq__(self, other): # compares a Match object with another Match object
        if isinstance(other, Match): # checks if other is a Match object
            return (self.avd == other.avd and
                    self.contact1 == other.contact1 and
                    self.contact2 == other.contact2 and
                    self.d3d4 == other.d3d4)
        return False