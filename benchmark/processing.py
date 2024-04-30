import matplotlib.pyplot as plt

file = './output.txt'

with open(file) as f:
    proteins = []
    for line in f:
        line = line.strip()
        if line.startswith("Detecting"):
            id = line[-6:]
            protein = [id]
        elif line.startswith("Protein"):
            size = int(line.split(":")[1])
            protein.append(size)
        # elif line.startswith("Contacts"):
        #     contacts = int(line.split(":")[1])
        #     proteins[id].append(contacts)
        #     protein.append(contacts)
        elif line.startswith("Time"):
            time = float(line.split(":")[1])
            protein.append(f"{time:.4f}")
        if len(protein) == 3:
            proteins.append(protein)

# Extracting data
ids = [entry[0] for entry in proteins]
sizes = [entry[1] for entry in proteins]
times = [float(entry[2]) for entry in proteins]

# Creating bar plot
plt.figure(figsize=(10, 6))
plt.bar(sizes, times)
plt.xlabel('Size')
plt.ylabel('Time')
plt.title('Time vs Size')
plt.xticks(rotation=45)
plt.show()