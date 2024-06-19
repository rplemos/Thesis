import matplotlib.pyplot as plt
import sys

file = sys.argv[1]

with open(file) as f:
    proteins = {}
    for line in f:
        line = line.strip()
        if line.startswith("Detecting"):
            id = line[-6:]
            proteins[id] = []
        elif "True" in line:
            size = (line.split(":")[1][1:])
            proteins[id].append(size)
        elif line.startswith("Number of contacts"):
            contacts = int(line.split(":")[1][1:])
            proteins[id].append(contacts)
        elif line.startswith("Contact"):
            time = float(line.split(":")[1][1:])
            proteins[id].append(f"{time:.4f}")  

ids = list(proteins.keys())
sizes = [float(entry[0]) for entry in proteins.values()]
contacts = [float(entry[1]) for entry in proteins.values()]
times = [float(entry[2]) for entry in proteins.values()]

# Create a scatter plot to visualize the relationship between contacts and times
plt.figure(figsize=(12, 8))

# Plot the data
plt.scatter(sizes, times, color='blue', edgecolor='k', alpha=0.7)

# Add labels and title
plt.xlabel('Protein Size', fontsize=14)
plt.ylabel('Time (s)', fontsize=14)
plt.title('Time vs Protein Size', fontsize=16)

# Annotate each point with the protein ID
for contact, time, id in zip(sizes, times, ids):
    plt.annotate(id, (contact, time), textcoords="offset points", xytext=(0, 10), ha='center', fontsize=9)

# Show the plot
plt.tight_layout()
plt.show()
