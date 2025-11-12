import numpy as np
import matplotlib.pyplot as plt
import sys

# Define the file name
file_name = 'output.txt'
directory = 'plots/'
output_image_file = 'semilogy_plot.png'
output_path = directory + output_image_file

try:
    # Use numpy.loadtxt to read the data efficiently.
    # It automatically handles the variable spacing and exponential notation.
    # 'data' will be a 2D numpy array where each row is a line from the file.
    data = np.loadtxt(file_name)

    # ----------------------------------------------
    # 1. Separate the columns
    # The first column (index 0) is the X-axis (data[0] in Python/NumPy slicing)
    # The second column (index 1) is the first Y-data (data[1])
    # The third column (index 2) is the second Y-data (data[2])
    X = data[:, 0]  # All rows, first column
    Y1 = data[:, 1] # All rows, second column
    Y2 = data[:, 2] # All rows, third column
    
    # Filter out any zero values, as log(0) is undefined and will cause errors
    # in the semilogy plot. We use a mask for Y1 and Y2.
    # This ensures that only positive, non-zero data points are plotted on the log scale.
    mask1 = Y1 > 0
    mask2 = Y2 > 0

    # Apply the masks to X and Y data for plotting
    X1_plot, Y1_plot = X[mask1], Y1[mask1]
    X2_plot, Y2_plot = X[mask2], Y2[mask2]

    # ----------------------------------------------
    # 2. Create the Semilogy Plot
    plt.figure(figsize=(10, 6))

    # Plot the second column (Y1) data using semilogy (log scale on Y-axis)
    # We use markers for visibility since there are few data points.
    plt.semilogy(X1_plot, Y1_plot, 'o-', label='QR-version', color='darkblue')

    # Plot the third column (Y2) data using semilogy
    plt.semilogy(X2_plot, Y2_plot, 's--', label='RQ-version', color='red')

    # ----------------------------------------------
    # 3. Add labels, title, and styling
    plt.title(f'Numerical accuracy of Both methods', fontsize=16)
    plt.xlabel('Dimension', fontsize=12)
    plt.ylabel('2-Norm of eigenvalue errorvector', fontsize=12)

    # Add a legend to differentiate the lines
    plt.legend(frameon=True, shadow=True, fontsize=10)

    # Add grid lines for better readability on the log scale
    plt.grid(True, which="both", ls="--", linewidth=0.5)

    # Use tight layout to prevent labels from being cut off
    plt.tight_layout()

    # Save the plot to a file instead of displaying it interactively,
    # resolving the FigureCanvasAgg UserWarning.
    plt.savefig(output_path)
    print(f"Plot successfully saved to {output_path}")


except FileNotFoundError:
    print(f"Error: The file '{file_name}' was not found.")
    print("Please make sure 'output.txt' is in the same directory as the script.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")