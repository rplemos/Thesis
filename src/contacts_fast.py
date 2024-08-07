from math import dist
from timeit import default_timer as timer
from numpy import dot, arccos, degrees
from numpy.linalg import norm
from classes import Contact, Match
import conditions
import distances

def fast_contacts(protein, fast, maximum_distances):
    start = timer()
    
    residues = list(protein.get_residues())
    contacts = []
    
    for i, residue1 in enumerate(residues):
        for j, residue2 in enumerate(residues[i+1:], start=i+1):
            
            if residue1.resnum == residue2.resnum and residue1.chain.id == residue2.chain.id: # ignores same residue
                continue
            
            if len(residue1.atoms) > 1 and len(residue2.atoms) > 1:
                ca1, ca2 = residue1.atoms[1], residue2.atoms[1] # alpha carbons

                distance_ca = dist((ca1.x, ca1.y, ca1.z), (ca2.x, ca2.y, ca2.z))
            
                # pair = tuple(sorted((residue1.resname, residue2.resname)))
                # if distance_ca > (distances.distances[pair] + 1):
                #     continue

                if distance_ca > 20:
                    continue
                
            else:
                continue              
            
            # CHECKING FOR AROMATIC STACKINGS
            if residue1.ring and residue2.ring:
                ring1, ring2 = residue1.atoms[-1], residue2.atoms[-1] # RNG atoms
                distance = dist((ring1.x, ring1.y, ring1.z), (ring2.x, ring2.y, ring2.z))
                angle = calc_angle(residue1.normal_vector, residue2.normal_vector)
                if distance >= 2 and distance <= 5: # within aromatic stacking limits
                    if (160 <= angle < 180) or (0 <= angle < 20):
                        stack_type = "-parallel"
                    elif (80 <= angle < 100):
                        stack_type = "-perpendicular"
                    else:
                        stack_type = "-other"

                    contact = Contact(protein.id, residue1.chain.id, residue1.resnum, residue1.resname, ring1.atomname, 
                                    protein.id, residue2.chain.id, residue2.resnum, residue2.resname, ring2.atomname, 
                                    float(f"{distance:.2f}"), "stacking"+stack_type, ring1, ring2)
                    
                    contacts.append(contact)
                    
            for atom1 in residue1.atoms:
                for atom2 in residue2.atoms:
                    name1 = f"{atom1.residue.resname}:{atom1.atomname}" # matches the pattern from contacts dictionary
                    name2 = f"{atom2.residue.resname}:{atom2.atomname}"
                    
                    if name1 in conditions.contact_types and name2 in conditions.contact_types: # excludes the RNG atom and any different other
                        
                        #if atom1.atomname != 'CA' and atom2.atomname != 'CA': # no need to calculate again for alpha carbons
                            distance = dist((atom1.x, atom1.y, atom1.z), (atom2.x, atom2.y, atom2.z))
                            
                            if distance <= 6: # max distance for contacts
                                contact_types = []
                                for contact_type, distance_range in conditions.categories.items():
                                    
                                    if not fast:
                                        if contact_type == 'hydrogen_bond' or contact_type == 'hydrophobic':
                                            continue                                
                                    
                                    if contact_type == 'hydrogen_bond' and (abs(residue2.resnum - residue1.resnum) <= 3): # skips alpha-helix for h-bonds
                                        continue
                                    
                                    if distance_range[0] <= distance <= distance_range[1]: # fits the range
                                        if conditions.contact_conditions[contact_type](name1, name2): # fits the type of contact
                                            contact_types.append(contact_type)

                                if contact_types:                                                                                                  
                                    contact = Contact(protein.id, residue1.chain.id, residue1.resnum, residue1.resname, atom1.atomname, 
                                                    protein.id, residue2.chain.id, residue2.resnum, residue2.resname, atom2.atomname, 
                                                    float(f"{distance:.2f}"), contact_types[0], atom1, atom2)
                                    
                                    contacts.append(contact)
                                        
                                    ####################
                                    # BLOCK FOR CONSTRUCTING MAXIMUM DISTANCES LIST                   
                                    if (residue1.resname, residue2.resname) not in maximum_distances:
                                        maximum_distances[residue1.resname, residue2.resname] = [float(f"{distance_ca:.2f}"), residue1.resnum, residue2.resnum, protein.id, atom1.atomname, atom2.atomname, residue1.chain.id, residue2.chain.id, contact_types[0]]
                                        #print(f"Creating {residue1.resname, residue2.resname} : {distance_ca}")
                                    elif distance_ca > maximum_distances[residue1.resname, residue2.resname][0]:
                                        #print(f"Changing {residue1.resname, residue2.resname} from {maximum_distances[residue1.resname, residue2.resname][0]} to {distance_ca}")
                                        maximum_distances[residue1.resname, residue2.resname] = [float(f"{distance_ca:.2f}"), residue1.resnum, residue2.resnum, protein.id, atom1.atomname, atom2.atomname, residue1.chain.id, residue2.chain.id, contact_types[0]]
                                    # ###################
                                                                                                            
                    else: # for control over non-standard atom names
                        pass
                        #print(f"Unknown atom: {name1} or {name2}")      
    
    end = timer()
    current_time = end - start

    return contacts, current_time, maximum_distances


def avd(contact_list_protein1, contact_list_protein2, cutoff):
    
    start = timer()
    
    match_list = []
    
    for contact1 in contact_list_protein1:
        atom1 = contact1.atom_object1
        atom2 = contact1.atom_object2
        p1_coords = atom1.x, atom1.y, atom1.z
        p2_coords = atom2.x, atom2.y, atom2.z
        
        for contact2 in contact_list_protein2:
            
            # FOR TESTING
            # if contact1.residueatom1 == "13VAL:CB" and contact1.residueatom2 == "123PHE:CD1":
            #     if contact2.residueatom1 == "94VAL:CG1" and contact2.residueatom2 == "98LEU:CD2":
            #         q1_coords = contact2.atom1.x, contact2.atom1.y, contact2.atom1.z
            #         q2_coords = contact2.atom2.x, contact2.atom2.y, contact2.atom2.z
                    
            #         d1 = dist(p1_coords, q1_coords)  # p1 x q1
            #         if d1 > cutoff * 3:
            #             continue                
            #         d2 = dist(p2_coords, q2_coords)  # p2 x q2
            #         d3 = dist(p1_coords, q2_coords)  # p1 x q2
            #         d4 = dist(p2_coords, q1_coords)  # p2 x q1

            #         avd = min(((d1 + d2)/2),((d3 + d4)/2))
            #         print(avd)
            #         print(contact1.type, contact2.type)
            #         print(contact1.to_list(), contact2.to_list())
            
            atom3 = contact2.atom_object1
            atom4 = contact2.atom_object2
            q1_coords = atom3.x, atom3.y, atom3.z
            q2_coords = atom4.x, atom4.y, atom4.z
            
            d1 = dist(p1_coords, q1_coords)  # p1 x q1
            if d1 > cutoff * 3:
                continue                
            d2 = dist(p2_coords, q2_coords)  # p2 x q2
            d3 = dist(p1_coords, q2_coords)  # p1 x q2
            d4 = dist(p2_coords, q1_coords)  # p2 x q1

            avd = min(((d1 + d2)/2),((d3 + d4)/2))
            
            if avd < cutoff:
                #if ('hydrogen_bond' not in contact1.type or 'hydrogen_bond' not in contact2.type):
                    d3d4 = True if (d1+d2)/2 > (d3+d4)/2 else False
                    match = Match(float(f"{avd:.2f}"), contact1, contact2, d3d4)
                    match_list.append(match)
    
    if len(match_list) == 0:
        return None, None, None
    
    average_avd = sum(single_match.avd for single_match in match_list) / len(match_list)
    contact_matches = len(match_list)

    end = timer()
    print(f"AVD Time elapsed (old): {end - start}\n")

    return match_list, average_avd, contact_matches


def show_contacts(contacts):
    category_counts = {}

    # Iterate over the output and count occurrences for each category
    for contact in contacts:
        category = contact.type
        category_counts[category] = category_counts.get(category, 0) + 1
    
    sorted_categories = sorted(category_counts.items(), key=lambda x: x[1])

    for category, count in sorted_categories:
        print(f"\nNumber of {category} occurrences:", count)
        if count <= 9999:
            print(f"All entries for {category}:")
            for entry in contacts:
                if entry.type == category:
                    print("\t",entry.print_text())

def calc_angle(vector1, vector2):
    dot_product = dot(vector1, vector2)
    magnitude_product = norm(vector1) * norm(vector2) # normalizes the dot product
    angle = arccos(dot_product / magnitude_product) # angle in radians   
    
    return degrees(angle)
