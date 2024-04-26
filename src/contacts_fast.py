from math import dist
from timeit import default_timer as timer
from numpy import dot, arccos, degrees
from numpy.linalg import norm
from classes import Contact, Match
import conditions

def fast_contacts(protein1, protein2):
    start = timer()
    
    residues1 = list(protein1.get_residues())
    residues2 = list(protein2.get_residues())
    contacts = []
    
    # FOR GETTING ALL THE CHAINS ON THE REFERENCE PROTEIN, AND SETTING TO COMPARE ONLY TO THEM
    chains = [chain.id for chain in protein1.chains]
    print(f"Chains to be analyzed: {chains}")
    
    big_ones = ["ARG", "LYS", "GLU", "PHE"]
    
    for i, residue1 in enumerate(residues1):
        for j, residue2 in enumerate(residues2[i+1:], start=i+1):
            residue1 = residues1[i]
            residue2 = residues2[j] 
            if residue1.chain.id in chains and residue2.chain.id in chains:
                ca1, ca2 = residue1.atoms[1], residue2.atoms[1] # alpha carbons
                distance_ca = dist((ca1.x, ca1.y, ca1.z), (ca2.x, ca2.y, ca2.z))
                if distance_ca > 20:
                    continue # skips the current residue 2
                elif residue1.resname not in big_ones and residue2.resname not in big_ones and distance_ca > 13:
                    # print(residue1.chain.id, residue1.resname, residue1.resnum, residue1.chain.id, residue2.resname, residue2.resnum, distance_ca)
                    continue
                
                # CHECKING FOR AROMATIC STACKINGS
                if residue1.ring and residue2.ring:
                    ring1, ring2 = residue1.atoms[-1], residue2.atoms[-1] # RNG atoms
                    distance = dist((ring1.x, ring1.y, ring1.z), (ring2.x, ring2.y, ring2.z))
                    angle = calc_angle(residue1.normal_vector, residue2.normal_vector)
                    if distance >= 2 and distance <= 5: # within aromatic stacking limits
                        if (160 <= angle < 180) or (0 <= angle < 20):
                            stack_type = "-parallel"
                            #print(f"Parallel.     \t Distance: {distance:.2f}. Angle ({residue1.chain.id}:{residue1.resnum}{residue1.resname} - {residue2.chain.id}:{residue2.resnum}{residue2.resname}): {angle:.2f}")
                        elif (80 <= angle < 100):
                            stack_type = "-perpendicular"
                            #print(f"Perpendicular.\t Distance: {distance:.2f}. Angle ({residue1.chain.id}:{residue1.resnum}{residue1.resname} - {residue2.chain.id}:{residue2.resnum}{residue2.resname}): {angle:.2f}")
                        else:
                            stack_type = "-other"
                            #print(f"Unknown.      \t Distance: {distance:.2f}. Angle ({residue1.chain.id}:{residue1.resnum}{residue1.resname} - {residue2.chain.id}:{residue2.resnum}{residue2.resname}): {angle:.2f}")
                                                           
                        contact = Contact(f"{protein1.id}:{residue1.chain.id}", f"{residue1.resnum}{residue1.resname}:{ring1.atomname}", 
                                        f"{protein2.id}:{residue2.chain.id}", f"{residue2.resnum}{residue2.resname}:{ring2.atomname}", 
                                        distance, ["stacking"+stack_type], ring1, ring2)
                        
                        contacts.append(contact)
                        
                for atom1 in residue1.atoms:
                    for atom2 in residue2.atoms:
                        name1 = f"{atom1.residue.resname}:{atom1.atomname}" # matches the pattern from contacts dictionary
                        name2 = f"{atom2.residue.resname}:{atom2.atomname}"
                        if name1 in conditions.contact_types and name2 in conditions.contact_types: # excludes the RNG atom and any different other
                            if atom1.atomname != 'CA' and atom2.atomname != 'CA': # no need to calculate again for alpha carbons
                                distance = dist((atom1.x, atom1.y, atom1.z), (atom2.x, atom2.y, atom2.z))
                            if distance <= 6: # max distance for contacts
                                contact_types = []
                                for contact_type, distance_range in conditions.categories.items():
                                    if contact_type == 'hydrogen_bond' and (residue2.resnum - residue1.resnum <= 3): # skips alpha-helix for h-bonds
                                        continue
                                    if distance_range[0] <= distance <= distance_range[1]: # fits the range
                                        if conditions.contact_conditions[contact_type](name1, name2): # fits the type of contact
                                            # if distance_ca > 13 and residue1.resname == 'ARG' and residue2.resname == 'ARG': # checking the weird distant ones
                                            #     print(distance, distance_ca, residue1.chain.id, residue1.resnum, name1, residue2.chain.id, residue2.resnum, name2)

                                            contact_types.append(contact_type)

                                if contact_types:
                                    contact = Contact(f"{protein1.id}:{residue1.chain.id}", f"{residue1.resnum}{name1}", 
                                                    f"{protein2.id}:{residue2.chain.id}", f"{residue2.resnum}{name2}", 
                                                    distance, contact_types, atom1, atom2)
                                    contacts.append(contact)
                                                                                        
                        else: # for control over non-standard atom names
                            pass
                            #print(f"Unknown atom: {name1} or {name2}")      
    end = timer()
    print(f"Time elapsed: {end - start}\n")
    return contacts

def avd(contact_list_protein1, contact_list_protein2, cutoff, match_list_new):
    start = timer()
    
    match_list = []
    
    for contact1 in contact_list_protein1:
        p1_coords = contact1.atom1.x, contact1.atom1.y, contact1.atom1.z
        p2_coords = contact1.atom2.x, contact1.atom2.y, contact1.atom2.z
        
        for contact2 in contact_list_protein2:
            q1_coords = contact2.atom1.x, contact2.atom1.y, contact2.atom1.z
            q2_coords = contact2.atom2.x, contact2.atom2.y, contact2.atom2.z
            
            d1 = dist(p1_coords, q1_coords)  # p1 x q1
            if d1 > cutoff * 5:
                continue                
            d2 = dist(p2_coords, q2_coords)  # p2 x q2
            d3 = dist(p1_coords, q2_coords)  # p1 x q2
            d4 = dist(p2_coords, q1_coords)  # p2 x q1

            avd = min(((d1 + d2)/2),((d3 + d4)/2))
            
            if avd < cutoff and 'hydrophobic' not in contact1.type and 'hydrophobic' not in contact2.type: # means that it's a match
                match = Match(avd, contact1.to_list(), contact2.to_list())
                if match not in match_list_new:
                    match_list.append(match)
    
    if len(match_list) == 0:
        return None, None, None
    
    average_avd = sum(single_match.avd for single_match in match_list) / len(match_list)
    contact_matches = len(match_list)

    end = timer()
    print(f"AVD Time elapsed (old): {end - start}\n")

    return match_list, average_avd, contact_matches

def new_avd(contact_list_protein1, contact_list_protein2, cutoff):
    start = timer()
    
    match_list = []
        
    for contact1 in contact_list_protein1:
        p1_coords = contact1.atom1.x, contact1.atom1.y, contact1.atom1.z
        p2_coords = contact1.atom2.x, contact1.atom2.y, contact1.atom2.z
        
        for contact2 in contact_list_protein2:
            q1_coords = contact2.atom1.x, contact2.atom1.y, contact2.atom1.z
            q2_coords = contact2.atom2.x, contact2.atom2.y, contact2.atom2.z
            
            d1 = dist(p1_coords, q1_coords)  # p1 x q1
            if d1 > (cutoff * 3):
                continue
            d2 = dist(p2_coords, q2_coords)  # p2 x q2

            avd = ((d1 + d2)/2)
            
            if avd < cutoff and 'hydrophobic' not in contact1.type and 'hydrophobic' not in contact2.type: # means that it's a match
                match = Match(avd, contact1.to_list(), contact2.to_list())
                match_list.append(match)
        
    if len(match_list) == 0:
        return None, None, None
    
    average_avd = sum(single_avd.avd for single_avd in match_list) / len(match_list)
    contact_matches = len(match_list)

    end = timer()
    print(f"AVD Time elapsed (new): {end - start}\n")

    return match_list, average_avd, contact_matches

def show_contacts(contacts):
    category_counts = {}

    # Iterate over the output and count occurrences for each category
    for contact in contacts:
        category = tuple(contact.type)
        category_counts[category] = category_counts.get(category, 0) + 1
    
    sorted_categories = sorted(category_counts.items(), key=lambda x: x[1])

    for category, count in sorted_categories:
        print(f"\nNumber of {str(category)[1:-1]} occurrences:", count)
        if count <= 3:
            print(f"All entries for {str(category)[1:-1]}:")
            for entry in contacts:
                if tuple(entry.type) == category:
                    print("\t",entry.to_list())

    print(f"\nTotal number of contacts: {len(contacts)}\n")    

def calc_angle(vector1, vector2):
    dot_product = dot(vector1, vector2)
    magnitude_product = norm(vector1) * norm(vector2) # normalizes the dot product
    angle = arccos(dot_product / magnitude_product) # angle in radians   
    
    return degrees(angle)

# For testing without main:

# #file = ["human_133l.pdb"]
# #file = ["human_133l.pdb", "chick_132l_aligned_rotate.pdb"]
# #file = ["human_133l.pdb", "8uw4.pdb"]
# file = ["8uw4.pdb"]
# parsed_proteins = pdb_parser.parse_pdb(file)
         
# distances = fast_contacts(parsed_proteins[0], parsed_proteins[0])
