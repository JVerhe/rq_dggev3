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
        if not os.path.isdir(INPUT_DIR):
            print(f"Error: Input directory '{INPUT_DIR}' not found.")
            sys.exit(1)

        all_files = os.listdir(INPUT_DIR)
        txt_files = sorted([f for f in all_files if f.endswith('.txt')])

        if not txt_files:
            print(f"No .txt files found in the '{INPUT_DIR}' directory. Exiting.")
            sys.exit(1)

        print(f"\nFiles found in '{INPUT_DIR}/':")
        for i, filename in enumerate(txt_files):
            print(f"  [{i + 1}] {filename}")

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
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    input_filename = os.path.basename(file_path)
    base_name = os.path.splitext(input_filename)[0]
    output_image_file = f'{base_name}_plot.png' 
    output_path = os.path.join(OUTPUT_DIR, output_image_file)

    try:
        data = np.loadtxt(file_path)

        X = data[:, 0]    # Column 1 (Index 0): Dimension
        
        # Accuracy data (Columns 2 and 3)
        Y_acc_1 = data[:, 1] # Index 1: Accuracy 1 (QR-version)
        Y_acc_2 = data[:, 2] # Index 2: Accuracy 2 (RQ-version)
        
        # Execution time data (Columns 4 and 5)
        Y_time_1 = data[:, 3] # Index 3: Time 1 (QR-version)
        Y_time_2 = data[:, 4] # Index 4: Time 2 (RQ-version)

        mask_acc_1 = Y_acc_1 > 0
        mask_acc_2 = Y_acc_2 > 0

        X1_acc_plot, Y1_acc_plot = X[mask_acc_1], Y_acc_1[mask_acc_1]
        X2_acc_plot, Y2_acc_plot = X[mask_acc_2], Y_acc_2[mask_acc_2]

        fig, (ax_acc, ax_time) = plt.subplots(1, 2, figsize=(16, 6)) 
        
        fig.suptitle(f'Accuracy and Execution Time ({input_filename})', fontsize=16)
        
        ax_acc.semilogy(X1_acc_plot, Y1_acc_plot, 'o-', label='QR-version', color='darkblue')
        ax_acc.semilogy(X2_acc_plot, Y2_acc_plot, 's--', label='RQ-version', color='red')

        ax_acc.set_title('Numerical Accuracy', fontsize=14)
        ax_acc.set_xlabel('Dimension', fontsize=12)
        if base_name.contains("RMSD"):
            ax_acc.set_ylabel('Root Mean Square Deviation', fontsize=12)
        else:
            ax_acc.set_ylabel('Max Relative Error', fontsize=12)
        ax_acc.legend(frameon=True, shadow=True, fontsize=10)
        ax_acc.grid(True, which="both", ls="--", linewidth=0.5)

        ax_time.plot(X, Y_time_1, 'o-', label='QR-version', color='darkblue')
        ax_time.plot(X, Y_time_2, 's--', label='RQ-version', color='red')
        
        ax_time.set_title('Execution Time', fontsize=14)
        ax_time.set_xlabel('Dimension', fontsize=12)
        ax_time.set_ylabel('Execution Time (ms)', fontsize=12)
        ax_time.legend(frameon=True, shadow=True, fontsize=10)
        ax_time.grid(True, which="both", ls="--", linewidth=0.5)
        
        plt.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust rect to leave space for suptitle

        plt.savefig(output_path)
        print(f"\nPlot successfully saved to {output_path}")

    except FileNotFoundError:
        print(f"Error: The file '{input_filename}' was not found.")
    except Exception as e:
        print(f"An unexpected error occurred during plotting: {e}")


if __name__ == "__main__":
    
    if not os.path.exists(OUTPUT_DIR):
        print(f"Creating output directory: {OUTPUT_DIR}/")
        os.makedirs(OUTPUT_DIR)
        
    selected_path = select_data_file()
    
    generate_plot(selected_path)