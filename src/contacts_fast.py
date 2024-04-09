import math
from timeit import default_timer as timer
import numpy as np

# 'RES:ATOM':	[	Hydrophobic,	Aromatic,	Positive,	Negative,	Donor,	Acceptor	]
# 'ALA:CA':		[	0|1,			0|1,		0|1,		0|1,		0|1,	0|1			]
contact_types = {'ALA:N':[0, 0, 0, 0, 1, 0],
            'ALA:CA':[0, 0, 0, 0, 0, 0],
            'ALA:C':[0, 0, 0, 0, 0, 0],
            'ALA:O':[0, 0, 0, 0, 0, 1],
            'ALA:CB':[1, 0, 0, 0, 0, 0],
            'ARG:N':[0, 0, 0, 0, 1, 0],
            'ARG:CA':[0, 0, 0, 0, 0, 0],
            'ARG:C':[0, 0, 0, 0, 0, 0],
            'ARG:O':[0, 0, 0, 0, 0, 1],
            'ARG:CB':[1, 0, 0, 0, 0, 0],
            'ARG:CG':[1, 0, 0, 0, 0, 0],
            'ARG:CD':[0, 0, 0, 0, 0, 0],
            'ARG:NE':[0, 0, 1, 0, 1, 0],
            'ARG:CZ':[0, 0, 1, 0, 0, 0],
            'ARG:NH1':[0, 0, 1, 0, 1, 0],
            'ARG:NH2':[0, 0, 1, 0, 1, 0],
            'ASN:N':[0, 0, 0, 0, 1, 0],
            'ASN:CA':[0, 0, 0, 0, 0, 0],
            'ASN:C':[0, 0, 0, 0, 0, 0],
            'ASN:O':[0, 0, 0, 0, 0, 1],
            'ASN:CB':[1, 0, 0, 0, 0, 0],
            'ASN:CG':[0, 0, 0, 0, 0, 0],
            'ASN:OD1':[0, 0, 0, 0, 0, 1],
            'ASN:ND2':[0, 0, 0, 0, 1, 0],
            'ASP:N':[0, 0, 0, 0, 1, 0],
            'ASP:CA':[0, 0, 0, 0, 0, 0],
            'ASP:C':[0, 0, 0, 0, 0, 0],
            'ASP:O':[0, 0, 0, 0, 0, 1],
            'ASP:CB':[1, 0, 0, 0, 0, 0],
            'ASP:CG':[0, 0, 0, 0, 0, 0],
            'ASP:OD1':[0, 0, 0, 1, 0, 1],
            'ASP:OD2':[0, 0, 0, 1, 0, 1],
            'CYS:N':[0, 0, 0, 0, 1, 0],
            'CYS:CA':[0, 0, 0, 0, 0, 0],
            'CYS:C':[0, 0, 0, 0, 0, 0],
            'CYS:O':[0, 0, 0, 0, 0, 1],
            'CYS:CB':[1, 0, 0, 0, 0, 0],
            'CYS:SG':[0, 0, 0, 0, 1, 1],
            'GLN:N':[0, 0, 0, 0, 1, 0],
            'GLN:CA':[0, 0, 0, 0, 0, 0],
            'GLN:C':[0, 0, 0, 0, 0, 0],
            'GLN:O':[0, 0, 0, 0, 0, 1],
            'GLN:CB':[1, 0, 0, 0, 0, 0],
            'GLN:CG':[1, 0, 0, 0, 0, 0],
            'GLN:CD':[0, 0, 0, 0, 0, 0],
            'GLN:OE1':[0, 0, 0, 0, 0, 1],
            'GLN:NE2':[0, 0, 0, 0, 1, 0],
            'GLU:N':[0, 0, 0, 0, 1, 0],
            'GLU:CA':[0, 0, 0, 0, 0, 0],
            'GLU:C':[0, 0, 0, 0, 0, 0],
            'GLU:O':[0, 0, 0, 0, 0, 1],
            'GLU:CB':[1, 0, 0, 0, 0, 0],
            'GLU:CG':[1, 0, 0, 0, 0, 0],
            'GLU:CD':[0, 0, 0, 0, 0, 0],
            'GLU:OE1':[0, 0, 0, 1, 0, 1],
            'GLU:OE2':[0, 0, 0, 1, 0, 1],
            'GLY:N':[0, 0, 0, 0, 1, 0],
            'GLY:CA':[0, 0, 0, 0, 0, 0],
            'GLY:C':[0, 0, 0, 0, 0, 0],
            'GLY:O':[0, 0, 0, 0, 0, 1],
            'HIS:N':[0, 0, 0, 0, 1, 0],
            'HIS:CA':[0, 0, 0, 0, 0, 0],
            'HIS:C':[0, 0, 0, 0, 0, 0],
            'HIS:O':[0, 0, 0, 0, 0, 1],
            'HIS:CB':[1, 0, 0, 0, 0, 0],
            'HIS:CG':[0, 1, 0, 0, 0, 0],
            'HIS:ND1':[0, 1, 1, 0, 1, 1],
            'HIS:CD2':[0, 1, 0, 0, 0, 0],
            'HIS:CE1':[0, 1, 0, 0, 0, 0],
            'HIS:NE2':[0, 1, 1, 0, 1, 1],
            'ILE:N':[0, 0, 0, 0, 1, 0],
            'ILE:CA':[0, 0, 0, 0, 0, 0],
            'ILE:C':[0, 0, 0, 0, 0, 0],
            'ILE:O':[0, 0, 0, 0, 0, 1],
            'ILE:CB':[1, 0, 0, 0, 0, 0],
            'ILE:CG1':[1, 0, 0, 0, 0, 0],
            'ILE:CG2':[1, 0, 0, 0, 0, 0],
            'ILE:CD1':[1, 0, 0, 0, 0, 0],
            'LEU:N':[0, 0, 0, 0, 1, 0],
            'LEU:CA':[0, 0, 0, 0, 0, 0],
            'LEU:C':[0, 0, 0, 0, 0, 0],
            'LEU:O':[0, 0, 0, 0, 0, 1],
            'LEU:CB':[1, 0, 0, 0, 0, 0],
            'LEU:CG':[1, 0, 0, 0, 0, 0],
            'LEU:CD1':[1, 0, 0, 0, 0, 0],
            'LEU:CD2':[1, 0, 0, 0, 0, 0],
            'LYS:N':[0, 0, 0, 0, 1, 0],
            'LYS:CA':[0, 0, 0, 0, 0, 0],
            'LYS:C':[0, 0, 0, 0, 0, 0],
            'LYS:O':[0, 0, 0, 0, 0, 1],
            'LYS:CB':[1, 0, 0, 0, 0, 0],
            'LYS:CG':[1, 0, 0, 0, 0, 0],
            'LYS:CD':[1, 0, 0, 0, 0, 0],
            'LYS:CE':[0, 0, 0, 0, 0, 0],
            'LYS:NZ':[0, 0, 1, 0, 1, 0],
            'MET:N':[0, 0, 0, 0, 1, 0],
            'MET:CA':[0, 0, 0, 0, 0, 0],
            'MET:C':[0, 0, 0, 0, 0, 0],
            'MET:O':[0, 0, 0, 0, 0, 1],
            'MET:CB':[1, 0, 0, 0, 0, 0],
            'MET:CG':[1, 0, 0, 0, 0, 0],
            'MET:SD':[0, 0, 0, 0, 0, 1],
            'MET:CE':[1, 0, 0, 0, 0, 0],
            'PHE:N':[0, 0, 0, 0, 1, 0],
            'PHE:CA':[0, 0, 0, 0, 0, 0],
            'PHE:C':[0, 0, 0, 0, 0, 0],
            'PHE:O':[0, 0, 0, 0, 0, 1],
            'PHE:CB':[1, 0, 0, 0, 0, 0],
            'PHE:CG':[1, 1, 0, 0, 0, 0],
            'PHE:CD1':[1, 1, 0, 0, 0, 0],
            'PHE:CD2':[1, 1, 0, 0, 0, 0],
            'PHE:CE1':[1, 1, 0, 0, 0, 0],
            'PHE:CE2':[1, 1, 0, 0, 0, 0],
            'PHE:CZ':[1, 1, 0, 0, 0, 0],
            'PRO:N':[0, 0, 0, 0, 0, 0],
            'PRO:CA':[0, 0, 0, 0, 0, 0],
            'PRO:C':[0, 0, 0, 0, 0, 0],
            'PRO:O':[0, 0, 0, 0, 0, 1],
            'PRO:CB':[1, 0, 0, 0, 0, 0],
            'PRO:CG':[1, 0, 0, 0, 0, 0],
            'PRO:CD':[0, 0, 0, 0, 0, 0],
            'SER:N':[0, 0, 0, 0, 1, 0],
            'SER:CA':[0, 0, 0, 0, 0, 0],
            'SER:C':[0, 0, 0, 0, 0, 0],
            'SER:O':[0, 0, 0, 0, 0, 1],
            'SER:CB':[0, 0, 0, 0, 0, 0],
            'SER:OG':[0, 0, 0, 0, 1, 1],
            'THR:N':[0, 0, 0, 0, 1, 0],
            'THR:CA':[0, 0, 0, 0, 0, 0],
            'THR:C':[0, 0, 0, 0, 0, 0],
            'THR:O':[0, 0, 0, 0, 0, 1],
            'THR:CB':[0, 0, 0, 0, 0, 0],
            'THR:OG1':[0, 0, 0, 0, 1, 1],
            'THR:CG2':[1, 0, 0, 0, 0, 0],
            'TRP:N':[0, 0, 0, 0, 1, 0],
            'TRP:CA':[0, 0, 0, 0, 0, 0],
            'TRP:C':[0, 0, 0, 0, 0, 0],
            'TRP:O':[0, 0, 0, 0, 0, 1],
            'TRP:CB':[1, 0, 0, 0, 0, 0],
            'TRP:CG':[1, 1, 0, 0, 0, 0],
            'TRP:CD1':[0, 1, 0, 0, 0, 0],
            'TRP:CD2':[1, 1, 0, 0, 0, 0],
            'TRP:NE1':[0, 1, 0, 0, 1, 0],
            'TRP:CE2':[0, 1, 0, 0, 0, 0],
            'TRP:CE3':[1, 1, 0, 0, 0, 0],
            'TRP:CZ2':[1, 1, 0, 0, 0, 0],
            'TRP:CZ3':[1, 1, 0, 0, 0, 0],
            'TRP:CH2':[1, 1, 0, 0, 0, 0],
            'TYR:N':[0, 0, 0, 0, 1, 0],
            'TYR:CA':[0, 0, 0, 0, 0, 0],
            'TYR:C':[0, 0, 0, 0, 0, 0],
            'TYR:O':[0, 0, 0, 0, 0, 1],
            'TYR:CB':[1, 0, 0, 0, 0, 0],
            'TYR:CG':[1, 1, 0, 0, 0, 0],
            'TYR:CD1':[1, 1, 0, 0, 0, 0],
            'TYR:CD2':[1, 1, 0, 0, 0, 0],
            'TYR:CE1':[1, 1, 0, 0, 0, 0],
            'TYR:CE2':[1, 1, 0, 0, 0, 0],
            'TYR:CZ':[0, 1, 0, 0, 0, 0],
            'TYR:OH':[0, 0, 0, 0, 1, 1],
            'VAL:N':[0, 0, 0, 0, 1, 0],
            'VAL:CA':[0, 0, 0, 0, 0, 0],
            'VAL:C':[0, 0, 0, 0, 0, 0],
            'VAL:O':[0, 0, 0, 0, 0, 1],
            'VAL:CB':[1, 0, 0, 0, 0, 0],
            'VAL:CG1':[1, 0, 0, 0, 0, 0],
            'VAL:CG2':[1, 0, 0, 0, 0, 0], 
}


# RULES
# 1 - must be made by different residue atoms
# 2 - aromatic = aromatic + aromatic
# 3 - hydrogenb => aceptor + donor
# 4 - hidrophobic: hidrofobic + hidrofobic
# 5 - Repulsive: positive=>positive ou negative=>negative
# 6 - Atractive: positive=>negative ou negative=>positive
# 7 - salt_bridge: positive=>negative ou negative=>positive

# TODO: 
# - implement a way to reset the residue_pairs set when the chain changes
# - implement alpha-helix skipping to avoid false positives (3 residues minimum?) 
#       this is kinda done  
def fast_contacts(protein1, protein2):
    start = timer()
    categories = {
        'hydrophobic': (2, 4.5),
        'aromatic': (2, 4),
        'hydrogen_bond': (0, 3.9),
        'repulsive': (2, 6),
        'attractive': (2, 6),
        'salt_bridge': (0, 3.9),
        'disulfide_bond': (0, 2.8),
    }

    contact_conditions = {
        'disulfide_bond': lambda name1, name2: name1 == "CYS:SG" and name2 == "CYS:SG",
        'aromatic': lambda name1, name2: contact_types[name1][1] == 1 and contact_types[name2][1] == 1,
        'hydrogen_bond': lambda name1, name2: ((contact_types[name1][4] == 1 and contact_types[name2][5] == 1) or (contact_types[name1][5] == 1 and contact_types[name2][4] == 1)) and (residue2.resnum - residue1.resnum >= 3),       
        'hydrophobic': lambda name1, name2: contact_types[name1][0] == 1 and contact_types[name2][0] == 1,
        'repulsive': lambda name1, name2: (contact_types[name1][2] == 1 and contact_types[name2][2] == 1) or (contact_types[name1][3] == 1 and contact_types[name2][3] == 1),
        'attractive': lambda name1, name2: (contact_types[name1][2] == 1 and contact_types[name2][3] == 1) or (contact_types[name1][3] == 1 and contact_types[name2][2] == 1),
        'salt_bridge': lambda name1, name2: (contact_types[name1][2] == 1 and contact_types[name2][3] == 1) or (contact_types[name1][3] == 1 and contact_types[name2][2] == 1)
    }
    
    residues1 = list(protein1.get_residues())
    residues2 = list(protein2.get_residues())
    contacts = []
    
    # FOR GETTING ALL THE CHAINS ON THE REFERENCE PROTEIN, AND SETTING TO COMPARE ONLY TO THEM
    chains = [chain.id for chain in protein1.chains]
    print(f"Chains to be analyzed: {chains}")
    
    big_ones = ["ARG", "LYS", "GLU"]
    
    #residue_pairs = set()
    for i, residue1 in enumerate(residues1):
        for j, residue2 in enumerate(residues2[i+1:], start=i+1):
            residue1 = residues1[i]
            residue2 = residues2[j] 
            #residue_pair = tuple(sorted([residue1.resnum, residue2.resnum]))
            
            if residue1.chain.id in chains and residue2.chain.id in chains:
                ca1, ca2 = residue1.atoms[1], residue2.atoms[1] # alpha carbons
                distance_ca = math.dist((ca1.x, ca1.y, ca1.z), (ca2.x, ca2.y, ca2.z))
                if distance_ca > 20:
                    #residue_pairs.add(residue_pair)
                    continue # skips the current residue 2
                elif residue1.resname not in big_ones and residue2.resname not in big_ones and distance_ca > 13:
                    continue
                
                if residue1.ring and residue2.ring: # two aromatic residues
                    ring1, ring2 = residue1.atoms[-1], residue2.atoms[-1] # RNG atoms
                    distance = math.dist((ring1.x, ring1.y, ring1.z), (ring2.x, ring2.y, ring2.z))
                    angle = calc_angle(residue1.normal_vector, residue2.normal_vector)
                    if angle > 160:
                        pass
                        #print(f"{distance:.2f} between {residue1.resnum}{residue1.resname} and {residue2.resnum}{residue2.resname}: {angle:.2f}")
                    if distance >= 2 and distance <= 4: # within aromatic stacking limits
                        contacts.append([f"{protein1.id}:{residue1.chain.id}", f"{residue1.resnum}", 
                                            f"{protein2.id}:{residue2.chain.id}", f"{residue2.resnum}", 
                                            distance, 'stacking', ring1, ring2])
                        
                for atom1 in residue1.atoms:
                    for atom2 in residue2.atoms:
                        name1 = f"{atom1.residue.resname}:{atom1.atomname}" # matches the pattern from contacts dictionary
                        name2 = f"{atom2.residue.resname}:{atom2.atomname}"
                        if name1 in contact_types and name2 in contact_types: # excludes the RNG atom and any different other
                            if atom1.atomname != 'CA' and atom2.atomname != 'CA': # no need to calculate again for alpha carbons
                                distance = math.dist((atom1.x, atom1.y, atom1.z), (atom2.x, atom2.y, atom2.z))
                            if distance <= 6: # max distance for contacts
                                for contact_type, distance_range in categories.items():
                                    if distance_range[0] <= distance <= distance_range[1]:
                                        if contact_conditions[contact_type](name1, name2):
                                            if distance_ca > 14: # checking the weird distant ones
                                                #pass
                                                print(distance, distance_ca, residue1.chain.id, residue1.resnum, name1, residue2.chain.id, residue2.resnum, name2)
                                            to_append = [f"{protein1.id}:{residue1.chain.id}", f"{residue1.resnum}{name1}", 
                                                            f"{protein2.id}:{residue2.chain.id}", f"{residue2.resnum}{name2}", 
                                                            distance, contact_type, atom1, atom2]
                                            contacts.append(to_append) 
                                            
                        else: # for control over non-standard atom names
                            pass
                            #print(f"Unknown atom: {name1} or {name2}")      
    
    end = timer()
    print(f"Time elapsed: {end - start}\n")

    return contacts

def avd(contact_list_protein1, contact_list_protein2, cutoff):
    avd_list = []
    for contact1 in contact_list_protein1:
        p1 = contact1[6] # atom object
        p2 = contact1[7]
        for contact2 in contact_list_protein2:
            q1 = contact2[6]
            q2 = contact2[7]

            d1 = math.dist((p1.x, p1.y, p1.z), (q1.x, q1.y, q1.z)) # p1 x q1
            d2 = math.dist((p2.x, p2.y, p2.z), (q2.x, q2.y, q2.z)) # p2 x q2
            d3 = math.dist((p1.x, p1.y, p1.z), (q2.x, q2.y, q2.z)) # p1 x q2
            d4 = math.dist((p2.x, p2.y, p2.z), (q1.x, q1.y, q1.z)) # p2 x q1
            
            avd = min(((d1 + d2)/2),((d3 + d4)/2))
            
            if avd < cutoff:
                avd_list.append([avd, contact1[:6], contact2[:6]])
    
    if len(avd_list) == 0:
        return None, None, None
    
    average_avd = sum(single_avd[0] for single_avd in avd_list) / len(avd_list)
    contact_matches = len(avd_list)

    return avd_list, average_avd, contact_matches

def show_contacts(distances):
    category_counts = {}

    # Iterate over the output and count occurrences for each category
    for entry in distances:
        category = entry[5]
        category_counts[category] = category_counts.get(category, 0) + 1

    sorted_categories = sorted(category_counts.items(), key=lambda x: x[1])

    for category, count in sorted_categories:
        print(f"Number of '{category}' occurrences:", count)
        if count <= 1:
            print(f"All entries for '{category}':")
            for entry in distances:
                if entry[5] == category:
                    print("\t",entry)

    print(f"Total number of contacts: {len(distances)}\n")    

def calc_angle(vector1, vector2):
    dot_product = np.dot(vector1, vector2)
    magnitude_product = np.linalg.norm(vector1) * np.linalg.norm(vector2) # normalizes the dot product
    angle = np.arccos(dot_product / magnitude_product) # angle in radians   
    
    return np.degrees(angle)

# For testing without main:

# #file = ["human_133l.pdb"]
# #file = ["human_133l.pdb", "chick_132l_aligned_rotate.pdb"]
# #file = ["human_133l.pdb", "8uw4.pdb"]
# file = ["8uw4.pdb"]
# parsed_proteins = pdb_parser.parse_pdb(file)
         
# distances = fast_contacts(parsed_proteins[0], parsed_proteins[0])
