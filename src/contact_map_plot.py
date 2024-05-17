import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import mplcursors
import textwrap

from timeit import default_timer as timer


def plot_matrix(contact_list, chain_residues, protein_length):
            
    start = timer()
    
    #------------------DEFINITIONS-------------------------
    
    legend_handles = []
    scatter_points = []
    tick_positions = []
    tick_labels = []
    
    color_map = {
        'hydrophobic': 'brown',
        'attractive': 'red',
        'repulsive': 'blue',
        'hydrogen_bond': 'purple',
        'disulfide_bond': 'yellow',
        'stacking': 'green',
        'salt_bridge': 'orange',
    }
    
    matrix = np.empty((protein_length+1, protein_length+1), dtype=object)
    plot_area = np.ones((protein_length+1, protein_length+1), dtype=int)
    plt.imshow(plot_area, cmap='binary')
        
    plt.gca().invert_yaxis()
    plt.grid(True, alpha=0.5)
    
    # Get the limits of the plot area
    xlim = plt.xlim()
    ylim = plt.ylim()

    # Calculate the position for the text
    text_offset = 0.08
    text_x = xlim[0] - (xlim[1] - xlim[0]) * text_offset
    text_y = ylim[0] - (ylim[1] - ylim[0]) * text_offset
    
    for label, color in color_map.items():
        legend_handles.append(plt.Line2D([0], [0], marker='o', color='w', markersize=10, markerfacecolor=color, label=label))
    
    #plt.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(1, 1))
    plt.legend(handles=legend_handles, loc='lower right')

    #------------------PROCESSING-------------------------
    
    for current_chain in chain_residues.keys():
        for contact in contact_list:
            
            if contact.chain1 == current_chain:
                
                res1 = contact.residue_num1 + chain_residues[contact.chain1]
                res2 = contact.residue_num2 + chain_residues[contact.chain2]
                
                size = 10
                
                if matrix[res1][res2] is None:
                    matrix[res1][res2] = contact.print_values()
                else:
                    matrix[res1][res2] += contact.print_values()

                label = contact.type
                if 'stacking' in label:
                    color = color_map['stacking']
                    alpha = 1
                    size = 25
                elif 'hydrogen_bond' in label or 'hydrophobic' in label:
                    alpha = 0.5
                    color = color_map[label]
                elif 'disulfide' in label:
                    alpha = 1
                    size = 50
                else:
                    alpha = 1
                    color = color_map.get(label, 'black')
                    size = 25
                
                scatter_points.append([res1, res2, color, alpha, size])
        
        
    for i, (current_chain, x) in enumerate(chain_residues.items()):
        
        if i == len(chain_residues) - 1:
            x_next = protein_length
        else:
            x_next = list(chain_residues.values())[i + 1]
            
            plt.axvline(x=x_next, color='black', linestyle='-', alpha=0.2)
            plt.axhline(y=x_next, color='black', linestyle='-', alpha=0.2)
            
        plt.gca().add_patch(Rectangle((x, x), x_next - x, x_next - x, fill=True, facecolor='blue', edgecolor='black', alpha=0.1))
                
        median_pos = (x + x_next) / 2
        
        plt.text(median_pos, text_x, current_chain, ha='center')
        plt.text(text_y, median_pos, current_chain, va='center')

        tick_labels.append([i for i in range(0, (x_next - x), 20)])
        tick_positions.append([i for i in range(x, x_next, 20)])  
        
        
    tick_labels = [item for sublist in tick_labels for item in sublist]
    tick_positions = [item for sublist in tick_positions for item in sublist]
        
    plt.xticks(tick_positions, rotation=45, fontsize=7)
    plt.yticks(tick_positions, fontsize=7)
    
    plt.gca().set_xticklabels(tick_labels)
    plt.gca().set_yticklabels(tick_labels)    
                       
    x_values, y_values, colors, alphas, sizes = zip(*scatter_points)
        
    scatter = plt.scatter(x=x_values, y=y_values, c=colors, marker='s', alpha=alphas, s=sizes, zorder=3)
                                
    mplcursors.cursor(scatter, hover=True).connect(
        "add", lambda sel: (
            sel.annotation.set_text("\n".join(textwrap.wrap(str(matrix[int(sel.target[0])][int(sel.target[1])]), width=100)))
            if (int(sel.target[0]), int(sel.target[1])) in zip(x_values, y_values)
            else sel.annotation.set_text(""),
            sel.annotation.set_bbox({"pad": 10, "facecolor": "white", "edgecolor": "black"})
        )
    )
    
    plt.suptitle(f"Contact Matrix for {contact_list[0].id1}", fontsize=20)
    plt.title(f"Protein Size: {protein_length}\nNumber of Contacts: {len(contact_list)}", fontsize=10)
    
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