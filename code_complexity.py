import subprocess
import time
import matplotlib.pyplot as plt
import numpy as np

# Values of N to test
N_values = [10, 100, 200, 400, 800, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]

# List to store execution times
times = []

for N in N_values:
    # Generate the filename based on N
    filename = f"./input_data/ellipse_N_{N:05d}.gal"

    # Command to run the code
    cmd = ["./galsim", str(N), filename, "200", "1e-5", "8"]

    # List to store the execution times for each run
    run_times = []

    for _ in range(5):
        # Measure the execution time
        start_time = time.time()
        subprocess.run(cmd)
        end_time = time.time()

        # Store the execution time
        run_times.append((end_time - start_time))

    # Store the minimum execution time
    times.append(min(run_times))

# Plot the execution times
plt.figure(figsize=(10, 6))
plt.plot(N_values, times, 'o-', label='Execution time')
plt.plot(N_values, np.log(N_values), 'o-', label='log N')
plt.xlabel('N')
plt.ylabel('Execution time / log N')
plt.title('Execution time and log N vs N')
plt.grid(True)
plt.legend()

# Save the figure
plt.savefig('code_complexity.png')