import numpy as np
import matplotlib.pyplot as plt
import mplcursors

def contact_map(contact_list, protein_length):
    matrix = np.empty((protein_length+1, protein_length+1), dtype=object)
        
    for contact in contact_list:
        parsed = contact.to_contact()
        num_res1 = int(parsed[1][:-3])
        num_res2 = int(parsed[3][:-3])

        matrix[num_res1][num_res2] = parsed
        
    return matrix

def plot_matrix(matrix):

    # Define a mapping of values to colors
    color_map = {
        'hydrophobic': 'red',
        'attractive': 'blue',
        'hydrogen_bond': 'green',  # Add more mappings as needed
    }

    plt.imshow(np.ones_like(matrix, dtype=int), cmap='binary')  # Plot the mask
    plt.gca().invert_yaxis()  # Invert y-axis to match matrix indexing

    # Plot colored dots based on the last entry of each entry
    scatter_points = []
    for i, row in enumerate(matrix):
        if row is not None:
            for j, val in enumerate(row):
                if isinstance(val, list):
                    color = color_map.get(val[-1][0], 'black')  # Extract string from list and get color from mapping
                    scatter_point = plt.scatter(j, i, color=color)  # Plot colored dot
                    scatter_points.append(scatter_point)

    # Add annotations to display the value when hovering
    mplcursors.cursor(scatter_points, hover=True).connect(
        "add", lambda sel: sel.annotation.set_text(str(matrix[int(sel.target[1])][int(sel.target[0])]))
    )

    plt.show()