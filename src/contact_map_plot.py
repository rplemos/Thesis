import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator

import mplcursors
import textwrap

def plot_matrix(contact_list, chain_residues, protein_length):
    
    # Define a mapping of values to colors
    color_map = {
        'hydrophobic': 'brown',
        'attractive': 'red',
        'repulsive': 'blue',
        'hydrogen_bond': 'purple',
        'disulfide_bond': 'yellow',
        'stacking': (0, 1, 0),
    }

    plt.imshow(np.ones((protein_length, protein_length), dtype=int), cmap='binary')  # Plot the mask
    plt.gca().invert_yaxis()
    plt.grid(True)
    
    scatter_points = []
    
    matrix = np.empty((protein_length+1, protein_length+1), dtype=object)    
    
    for chain, length in chain_residues.items():
        print(chain, length, protein_length)
        
        current_chain = chain
        
        for contact in contact_list:
            
            if contact.chain1 == current_chain:
                contact.residue_num1 += chain_residues[contact.chain1]
                contact.residue_num2 += chain_residues[contact.chain2]
                
                label = contact.type
                size = 10
                marker = 's'
                if 'stacking' in label:
                    color = color_map['stacking']
                    zorder = 10
                elif 'hydrogen_bond' in label or 'hydrophobic' in label:
                    color = color_map[label]
                    zorder = 0
                else:
                    color = color_map.get(label, 'black')  # Use color map for other labels
                    zorder = 5
                
                scatter_point = plt.scatter(contact.residue_num1, contact.residue_num2, color=color, marker=marker, s=size, zorder=zorder)
                scatter_points.append(scatter_point)
                                
                if matrix[contact.residue_num1][contact.residue_num2] is None:
                    matrix[contact.residue_num1][contact.residue_num2] = contact.print_values()
                else:
                    matrix[contact.residue_num1][contact.residue_num2] += contact.print_values()

    
    mplcursors.cursor(scatter_points, hover=True).connect(
        "add", lambda sel: (sel.annotation.set_text("\n".join(textwrap.wrap(str(matrix[int(sel.target[0])][int(sel.target[1])]), width=50))),
                            sel.annotation.set_bbox({"pad": 10, "facecolor": "white", "edgecolor": "black"}))
    )

    legend_handles = []
    for label, color in color_map.items():
        legend_handles.append(plt.Line2D([0], [0], marker='o', color='w', markersize=10, markerfacecolor=color, label=label))

    plt.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(1, 1))
    
    # Get the limits of the plot area
    xlim = plt.xlim()
    ylim = plt.ylim()

    # Calculate the position for the text
    text_offset = 0.08  # Adjust this value as needed
    text_x = xlim[0] - (xlim[1] - xlim[0]) * text_offset
    text_y = ylim[0] - (ylim[1] - ylim[0]) * text_offset
    
    plt.gca().xaxis.set_major_locator(MultipleLocator(round(protein_length/20)))
    plt.gca().yaxis.set_major_locator(MultipleLocator(round(protein_length/20)))

    for i, (label, x) in enumerate(chain_residues.items()):
        if i == len(chain_residues) - 1:
            x_next = protein_length
        else:
            x_next = list(chain_residues.values())[i + 1]
        plt.gca().add_patch(Rectangle((x, x), x_next - x, x_next - x, fill=True, facecolor='blue', edgecolor='black', alpha=0.1))
                
        median_pos = (x + x_next) / 2
        
        plt.text(median_pos, text_x, label, ha='center')
        plt.text(text_y, median_pos, label, va='center')
    
    plt.show()   