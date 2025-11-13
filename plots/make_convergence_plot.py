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
    Reads the data from the given file_path, generates the accuracy (semilogy) 
    and execution time (linear) plots side-by-side, and saves them to the OUTPUT_DIR.
    """
    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Get just the filename for titles and output path
    input_filename = os.path.basename(file_path)
    base_name = os.path.splitext(input_filename)[0]
    # Update output file name to reflect two plots
    output_image_file = f'{base_name}_accuracy_and_time.png' 
    output_path = os.path.join(OUTPUT_DIR, output_image_file)

    try:
        # Use numpy.loadtxt to read the data efficiently.
        data = np.loadtxt(file_path)

        # ----------------------------------------------
        # 1. Separate the columns based on the original request and the new requirement
        X = data[:, 0]    # Column 1 (Index 0): Dimension
        
        # Accuracy data (Columns 2 and 3)
        Y_acc_1 = data[:, 1] # Index 1: Accuracy 1 (QR-version)
        Y_acc_2 = data[:, 2] # Index 2: Accuracy 2 (RQ-version)
        
        # Execution time data (Columns 4 and 5)
        Y_time_1 = data[:, 3] # Index 3: Time 1 (QR-version)
        Y_time_2 = data[:, 4] # Index 4: Time 2 (RQ-version)

        # Filter zero values for the accuracy plot (only log-scaled data needs filtering)
        mask_acc_1 = Y_acc_1 > 0
        mask_acc_2 = Y_acc_2 > 0

        # Apply masks
        X1_acc_plot, Y1_acc_plot = X[mask_acc_1], Y_acc_1[mask_acc_1]
        X2_acc_plot, Y2_acc_plot = X[mask_acc_2], Y_acc_2[mask_acc_2]
        # Time data doesn't require filtering unless times can be zero/negative. 
        # Assuming non-zero positive times for simplicity.

        # ----------------------------------------------
        # 2. Create Figure with Two Subplots
        
        # Create a figure and a set of subplots (1 row, 2 columns)
        fig, (ax_acc, ax_time) = plt.subplots(1, 2, figsize=(16, 6)) 
        
        # Set overall title for the figure
        fig.suptitle(f'Comparative Analysis: Accuracy and Execution Time ({input_filename})', fontsize=16)

        # --- Subplot 1: Accuracy (Semilogy Plot) ---
        
        # Plot accuracy data using semilogy
        ax_acc.semilogy(X1_acc_plot, Y1_acc_plot, 'o-', label='QR-version', color='darkblue')
        ax_acc.semilogy(X2_acc_plot, Y2_acc_plot, 's--', label='RQ-version', color='red')

        # Add labels, title, and styling for the accuracy plot
        ax_acc.set_title('Numerical Accuracy', fontsize=14)
        ax_acc.set_xlabel('Dimension', fontsize=12)
        ax_acc.set_ylabel('2-Norm of eigenvalue error vector (Log Scale)', fontsize=12)
        ax_acc.legend(frameon=True, shadow=True, fontsize=10)
        ax_acc.grid(True, which="both", ls="--", linewidth=0.5)

        # --- Subplot 2: Execution Time (Linear Plot) ---
        
        # Plot execution time data (standard linear plot is appropriate for time)
        ax_time.plot(X, Y_time_1, 'o-', label='QR-version', color='darkblue')
        ax_time.plot(X, Y_time_2, 's--', label='RQ-version', color='red')
        
        # Add labels, title, and styling for the time plot
        ax_time.set_title('Execution Time', fontsize=14)
        ax_time.set_xlabel('Dimension', fontsize=12)
        ax_time.set_ylabel('Execution Time (ms)', fontsize=12)
        ax_time.legend(frameon=True, shadow=True, fontsize=10)
        ax_time.grid(True, which="both", ls="--", linewidth=0.5)
        
        # Use tight layout to prevent labels from overlapping
        plt.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust rect to leave space for suptitle

        # Save the plot to the /plots directory
        plt.savefig(output_path)
        print(f"\nPlot successfully saved to {output_path}")

    except FileNotFoundError:
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