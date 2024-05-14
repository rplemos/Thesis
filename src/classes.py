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
    
    def count_residues(self):
        return self.residues[-1].resnum
        #return sum(1 for _ in self.get_residues())


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
    def __init__(self, id1, chain1, residue_num1, residue_name1, atom1, 
                 id2, chain2, residue_num2, residue_name2, atom2, 
                 distance, type, atom_object1, atom_object2):
        self.id1 = id1
        self.chain1 = chain1
        self.residue_num1 = residue_num1
        self.residue_name1 = residue_name1
        self.atom1 = atom1
        self.id2 = id2
        self.chain2 = chain2
        self.residue_num2 = residue_num2
        self.residue_name2 = residue_name2
        self.atom2 = atom2
        self.distance = distance
        self.type = type
        self.atom_object1 = atom_object1
        self.atom_object2 = atom_object2
        
    def get_values(self):
        all_values = list(self.__dict__.values())
        return all_values[:-2]
    
    def print_values(self):
        all_values = list(self.__dict__.values())
        return [f"{all_values[0]}:{all_values[1]}", f"{all_values[2]}{all_values[3]}:{all_values[4]}",
                f"{all_values[5]}:{all_values[6]}", f"{all_values[7]}{all_values[8]}:{all_values[9]}",
                all_values[10], all_values[11]]
    
    def print_contact(self):
        all_values = list(self.__dict__.values())
        return f"Distance between {all_values[0]}:{all_values[1]}-{all_values[2]}{all_values[3]}:{all_values[4]} and {all_values[5]}:{all_values[6]}-{all_values[7]}{all_values[8]}:{all_values[9]}: {all_values[10]} A. Types: {all_values[11]}"
    
    # def to_contact2(self):
    #     type = str(self.type)[2:-2].replace("'", "")
    #     return [self.idchain1.split(":")[1], self.residueatom1.split(":")[0], self.idchain2.split(":")[1], 
    #             self.residueatom2.split(":")[0], [f"{type} ({self.atom1.atomname}:{self.atom2.atomname})"]]
        
    # def to_contact(self):
    #     return [self.idchain1.split(":")[1], self.residueatom1.split(":")[0], self.idchain2.split(":")[1], 
    #             self.residueatom2.split(":")[0], self.type]
        

class Match:
    def __init__(self, avd, contact1, contact2, d3d4):
        self.avd = avd
        self.contact1 = contact1
        self.contact2 = contact2
        self.d3d4 = d3d4 
        
    # def __eq__(self, other): # compares a Match object with another Match object
    #     if isinstance(other, Match): # checks if other is a Match object
    #         return (self.avd == other.avd and
    #                 self.contact1 == other.contact1 and
    #                 self.contact2 == other.contact2 and
    #                 self.d3d4 == other.d3d4)
    #     return False
