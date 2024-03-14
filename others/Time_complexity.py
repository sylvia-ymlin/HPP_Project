import subprocess
import time
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np


# Change the executable inputs to run with different N
command_template = "./a.out {} {} 200 1e-5 0"

# Number of threads
N_array = [10, 20, 40, 80, 100, 200, 400, 800, 1000, 2000]
file_names = ["./input_data/ellipse_N_00010.gal",
             "./input_data/ellipse_N_00020.gal",
             "./input_data/ellipse_N_00040.gal",
             "./input_data/ellipse_N_00080.gal",
             "./input_data/ellipse_N_00100.gal",
             "./input_data/ellipse_N_00200.gal",
             "./input_data/ellipse_N_00400.gal",
             "./input_data/ellipse_N_00800.gal",
             "./input_data/ellipse_N_01000.gal",
             "./input_data/ellipse_N_02000.gal"]

execution_times = []
# Run the executable for each number of threads
for N, name in zip(N_array, file_names):
    thread_times = []
    
    for _ in range(5):  # Run the command several times
        command = command_template.format(N, name)
        start_time = time.time()
        subprocess.run(command, shell=True)  # Run the command
        end_time = time.time()
        execution_time = end_time - start_time
        thread_times.append(execution_time)
        time.sleep(0.1)
        
    execution_times.append(min(thread_times)/N)  # Only store the minimum execution time

# Define the 1/x function
def func(x, a, b):
    return a / x + b

# # Calculate ideal times
# ideal_times = [execution_times[0]/n for n in num_threads]

# # Fit the 1/x function to the ideal times
# popt, pcov = curve_fit(func, num_threads, ideal_times)

# Generate a range of x values for the fitted function
x = np.linspace(N_array[0], N_array[-1], 500)

# Plot the results
plt.scatter(N_array, execution_times)
# plt.plot(x, func(x, *popt), '--', label='Ideal 1/threads behavior', color='r')
# plt.xlabel('Number of Threads', fontsize=18)
# plt.ylabel('Execution Time (s)', fontsize=18)
# I don't think title is necessary, we can describe the plot in the report
# plt.tick_params(axis='both', which='major', labelsize=16)

# Change the name here if you modify N, else it would overwrite the previous plot
plt.savefig('time_complexity.png', dpi=300, bbox_inches='tight')

plt.show()