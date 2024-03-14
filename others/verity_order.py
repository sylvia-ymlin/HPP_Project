import numpy as np
import matplotlib.pyplot as plt

# Load error values from file
with open('error_values_order2.txt', 'r') as file:
    error_values = [float(line.strip()) for line in file]

# Calculate differences between consecutive error values
x = [error_values[i] - error_values[i-1] for i in range(1, len(error_values))]

# Generate time values from 0 to 0.019 with a step size of 0.001
dt_squared = np.arange(0, 0.019, 0.001)

# Plot error vs dt^2 on a linear scale
plt.figure(figsize=(10, 6))
plt.plot(dt_squared, x, marker='o', linestyle='-')  # Using linear scale
plt.xlabel('T')
plt.ylabel('Error')
plt.title('Convergence Study for Velocity Verlet Method')
plt.grid(True)



# Save the plot to an image file
plt.savefig('convergence_plot_order2.png')

# Show the plot
plt.show()
