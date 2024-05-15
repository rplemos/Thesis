import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator
import mplcursors
import textwrap

from timeit import default_timer as timer


def plot_matrix(contact_list, chain_residues, protein_length):
            
    start = timer()
    
    #------------------DEFINITIONS-------------------------
    
    legend_handles = []
    scatter_points = []
    matrix = np.empty((protein_length+1, protein_length+1), dtype=object)

    color_map = {
        'hydrophobic': 'brown',
        'attractive': 'red',
        'repulsive': 'blue',
        'hydrogen_bond': 'purple',
        'disulfide_bond': 'yellow',
        'stacking': 'green',
    }

    plt.imshow(np.ones((protein_length, protein_length), dtype=int), cmap='binary')  # Plot the mask
    plt.gca().invert_yaxis()
    plt.grid(True)
    
    # Get the limits of the plot area
    xlim = plt.xlim()
    ylim = plt.ylim()

    # Calculate the position for the text
    text_offset = 0.08
    text_x = xlim[0] - (xlim[1] - xlim[0]) * text_offset
    text_y = ylim[0] - (ylim[1] - ylim[0]) * text_offset
    
    for label, color in color_map.items():
        legend_handles.append(plt.Line2D([0], [0], marker='o', color='w', markersize=10, markerfacecolor=color, label=label))
    
    plt.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(1, 1))
    
    plt.gca().xaxis.set_major_locator(MultipleLocator(round(protein_length/20)))
    plt.gca().yaxis.set_major_locator(MultipleLocator(round(protein_length/20)))    
    
    #------------------PROCESSING-------------------------
    
    for i, (current_chain, x) in enumerate(chain_residues.items()):
        for contact in contact_list:
            
            if contact.chain1 == current_chain:

                contact.residue_num1 += chain_residues[contact.chain1]
                contact.residue_num2 += chain_residues[contact.chain2]
                
                if matrix[contact.residue_num1][contact.residue_num2] is None:
                    matrix[contact.residue_num1][contact.residue_num2] = contact.print_values()
                else:
                    matrix[contact.residue_num1][contact.residue_num2] += contact.print_values()

                label = contact.type
                if 'stacking' in label:
                    color = color_map['stacking']
                elif 'hydrogen_bond' in label or 'hydrophobic' in label:
                    color = color_map[label]
                else:
                    color = color_map.get(label, 'black')
                
                scatter_points.append([contact.residue_num1, contact.residue_num2, color])
                
        if i == len(chain_residues) - 1:
            x_next = protein_length
        else:
            x_next = list(chain_residues.values())[i + 1]
        plt.gca().add_patch(Rectangle((x, x), x_next - x, x_next - x, fill=True, facecolor='blue', edgecolor='black', alpha=0.1))
                
        median_pos = (x + x_next) / 2
        
        plt.text(median_pos, text_x, current_chain, ha='center')
        plt.text(text_y, median_pos, current_chain, va='center')        
                       
    x_values, y_values, colors = zip(*scatter_points)
    
    scatter = plt.scatter(x=x_values, y=y_values, c=colors, s=10, marker='s')
                                
    mplcursors.cursor(scatter, hover=True).connect(
        "add", lambda sel: (
            sel.annotation.set_text("\n".join(textwrap.wrap(str(matrix[int(sel.target[0])][int(sel.target[1])]), width=50)))
            if (int(sel.target[0]), int(sel.target[1])) in zip(x_values, y_values)
            else sel.annotation.set_text(""),
            sel.annotation.set_bbox({"pad": 10, "facecolor": "white", "edgecolor": "black"})
        )
    )
    
    plt.title(contact_list[0].id1)
    plt.show()
    
    end = timer()
    print(f"Contact map graph - Time elapsed: {end - start}\n")
    
    return matrix
    

def contact_matrix(contact_list, chain_residues, protein_length):
        
    matrix = np.empty((protein_length+1, protein_length+1), dtype=object)
    
    for i, (current_chain, x) in enumerate(chain_residues.items()):
        for contact in contact_list:            
            if contact.chain1 == current_chain:
                
                contact.residue_num1 += chain_residues[contact.chain1]
                contact.residue_num2 += chain_residues[contact.chain2]
                
                if matrix[contact.residue_num1][contact.residue_num2] is None:
                    matrix[contact.residue_num1][contact.residue_num2] = contact.print_values()
                else:
                    matrix[contact.residue_num1][contact.residue_num2] += contact.print_values()
    
    return matrix