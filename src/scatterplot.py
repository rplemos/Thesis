import matplotlib.pyplot as plt
import sys
import numpy as np
from sklearn.metrics import r2_score 

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
        elif line.startswith("Contact"):
            time = float(line.split(":")[1][1:])
            proteins[id].append(f"{time:.4f}")  
        # elif line.startswith("Number of contacts"):
        #     contacts = int(line.split(":")[1][1:])
        #     proteins[id].append(contacts)

ids = list(proteins.keys())
sizes = [float(entry[0]) for entry in proteins.values()]
times = [float(entry[1]) for entry in proteins.values()]
#contacts = [float(entry[1]) for entry in proteins.values()]

sizes = np.array(sizes)
times = np.array(times)

model_linear = np.poly1d(np.polyfit(sizes, times, 1)) 
model_quadratic = np.poly1d(np.polyfit(sizes, times, 2))
  
polyline = np.linspace(min(sizes), max(sizes), 100) 

linear_polyline = model_linear(polyline)
quadratic_polyline = model_quadratic(polyline) 

r2_linear = r2_score(times, model_linear(sizes))
print(f"R2 linear: {r2_linear}")

r2_quadratic = r2_score(times, model_quadratic(sizes))
print(f"R2 quadratic: {r2_quadratic}")

plt.figure(figsize=(12, 8))
plt.scatter(sizes, times, color='blue', edgecolor='k', alpha=0.7)
plt.plot(polyline, linear_polyline, label=f'Linear Fit\n$R^2 = {r2_linear:.2f}$') 
plt.plot(polyline, quadratic_polyline, label=f'Quadratic Fit\n$R^2 = {r2_quadratic:.2f}$')
plt.xlabel('Protein Size', fontsize=14)
plt.ylabel('Time (s)', fontsize=14)
plt.title('Time vs Protein Size', fontsize=16)

# Annotate each point with the protein ID
# for contact, time, id in zip(sizes, times, ids):
#     plt.annotate(id, (contact, time), textcoords="offset points", xytext=(0, 10), ha='center', fontsize=9)

plt.tight_layout()
plt.legend()
plt.show()

# models
print(model_linear)
print(model_quadratic) 