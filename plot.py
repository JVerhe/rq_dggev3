import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# --- Configuration ---
INPUT_DIR = 'results'
OUTPUT_DIR = 'plots'
# ---------------------

def select_data_file():
    """
    Lists .txt files in the INPUT_DIR, prompts the user to select one, 
    and returns the full path of the selected file.
    """
    try:
        # Ensure the input directory exists
        if not os.path.isdir(INPUT_DIR):
            print(f"Error: Input directory '{INPUT_DIR}' not found.")
            # Exit the script gracefully if the expected directory is missing
            sys.exit(1)

        # 1. Look inside the directory /results and filter for .txt files
        all_files = os.listdir(INPUT_DIR)
        txt_files = sorted([f for f in all_files if f.endswith('.txt')])

        if not txt_files:
            print(f"No .txt files found in the '{INPUT_DIR}' directory. Exiting.")
            sys.exit(1)

        # Write the names of all the files inside and ask for selection
        print(f"\nFiles found in '{INPUT_DIR}/':")
        for i, filename in enumerate(txt_files):
            print(f"  [{i + 1}] {filename}")

        # 2. Ask the user to select a txt file via a number input
        while True:
            try:
                selection = input("\nEnter the number of the file you want to plot: ")
                # Convert input to 0-based index
                index = int(selection) - 1
                
                if 0 <= index < len(txt_files):
                    selected_filename = txt_files[index]
                    full_file_path = os.path.join(INPUT_DIR, selected_filename)
                    print(f"Selected file: {selected_filename}")
                    return full_file_path
                else:
                    print(f"Invalid selection. Please enter a number between 1 and {len(txt_files)}.")
            except ValueError:
                print("Invalid input. Please enter a number.")

    except Exception as e:
        print(f"An error occurred during file selection: {e}")
        sys.exit(1)


def generate_plot(file_path):
    """
    Reads the data from the given file_path, generates the semilogy plot, 
    and saves it to the OUTPUT_DIR.
    """
    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Get just the filename for titles and output path
    input_filename = os.path.basename(file_path)
    base_name = os.path.splitext(input_filename)[0]
    output_image_file = f'{base_name}_semilogy.png'
    output_path = os.path.join(OUTPUT_DIR, output_image_file)

    try:
        # Use numpy.loadtxt to read the data efficiently.
        data = np.loadtxt(file_path)

        # ----------------------------------------------
        # 1. Separate the columns
        X = data[:, 0]  # All rows, first column
        Y1 = data[:, 1] # All rows, second column
        Y2 = data[:, 2] # All rows, third column
        
        # Filter out zero values (log(0) is undefined)
        mask1 = Y1 > 0
        mask2 = Y2 > 0

        # Apply the masks to X and Y data for plotting
        X1_plot, Y1_plot = X[mask1], Y1[mask1]
        X2_plot, Y2_plot = X[mask2], Y2[mask2]

        # ----------------------------------------------
        # 2. Create the Semilogy Plot
        plt.figure(figsize=(10, 6))

        # Plot the second column (Y1) data
        plt.semilogy(X1_plot, Y1_plot, 'o-', label='QR-version', color='darkblue')

        # Plot the third column (Y2) data
        plt.semilogy(X2_plot, Y2_plot, 's--', label='RQ-version', color='red')

        # ----------------------------------------------
        # 3. Add labels, title, and styling
        plt.title(f'Numerical Accuracy of Both Methods ({input_filename})', fontsize=16)
        plt.xlabel('Dimension', fontsize=12)
        plt.ylabel('2-Norm of eigenvalue errorvector', fontsize=12)

        # Add a legend and grid
        plt.legend(frameon=True, shadow=True, fontsize=10)
        plt.grid(True, which="both", ls="--", linewidth=0.5)

        plt.tight_layout()

        # Save the plot to the /plots directory
        plt.savefig(output_path)
        print(f"\nPlot successfully saved to {output_path}")

    except FileNotFoundError:
        # This error should ideally be caught in select_data_file, but included for robustness
        print(f"Error: The file '{input_filename}' was not found.")
    except Exception as e:
        print(f"An unexpected error occurred during plotting: {e}")


if __name__ == "__main__":
    # Main execution flow
    
    # First, check if the required output directory exists.
    if not os.path.exists(OUTPUT_DIR):
        print(f"Creating output directory: {OUTPUT_DIR}/")
        os.makedirs(OUTPUT_DIR)
        
    # Get the selected file path interactively
    selected_path = select_data_file()
    
    # Generate and save the plot
    generate_plot(selected_path)