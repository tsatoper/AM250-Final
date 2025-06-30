import csv
import numpy as np
import matplotlib.pyplot as plt

csv_file = 'runtime_results.csv'
numproc = []
N = []
timesteps = []
time_elapsed = []

with open(csv_file, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        numproc.append(int(row[0]))
        N.append(int(row[1]))
        timesteps.append(int(row[2]))
        time_elapsed.append(float(row[3]))

# Convert lists to numpy arrays for easier indexing
numproc = np.array(numproc)
N = np.array(N)
time_elapsed = np.array(time_elapsed)

plt.figure()

unique_numproc = np.unique(numproc)

for p in unique_numproc:
    indices = np.where(numproc == p)
    size = np.sqrt(numproc[indices]) * N[indices]
    times = time_elapsed[indices]
    print(f'{indices=}, {size=}')
    plt.plot(size, times, marker='o', label=f'numproc = {p}')

plt.xlabel('grid length')
plt.ylabel('Time elapsed (s)')
plt.title('Time elapsed vs grid length')
plt.legend()
plt.grid(True)
plt.savefig('time_vs_size_by_numproc.png')
print("Plot saved to time_vs_size_by_numproc.png")
