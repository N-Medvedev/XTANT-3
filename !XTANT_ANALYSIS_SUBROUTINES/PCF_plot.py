import os
import glob
import re
import shlex
import numpy as np
import matplotlib.pyplot as plt

def plot_single_rdf(file_path):
    """
    Parses and plots the RDF data from a single file with custom colors, 
    line styles, aspect ratio, and robust negative time extraction.
    """
    # Load the numeric data
    data = np.loadtxt(file_path)
    r = data[:, 0]  # First column: Pair separation distance
    num_curves = data.shape[1] - 1
    
    # Read the header to extract labels dynamically
    with open(file_path, 'r') as f:
        f.readline()  # Skip the first line
        second_line = f.readline().strip().lstrip('#').strip()
    
    try:
        tokens = shlex.split(second_line)
        if len(tokens) == data.shape[1]:
            labels = tokens[1:]
        else:
            labels = second_line.split()[-num_curves:]
    except Exception:
        labels = second_line.split()[-num_curves:]
        
    while len(labels) < num_curves:
        labels.append(f"Column {len(labels) + 1}")

    # Set up the figure with an 8:6 aspect ratio
    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
    
    # Define your line styles and colors here
    line_styles = ['-', '--', '-.']
    # You can use standard names ('red', 'blue', 'green') or hex codes
    colors = ['royalblue', 'forestgreen', 'orange']
    
    # Plot each column against the first one
    for i in range(num_curves):
        style = line_styles[i % len(line_styles)]
        color = colors[i % len(colors)]
        
        ax.plot(
            r, 
            data[:, i + 1], 
            label=labels[i], 
            linestyle=style, 
            color=color, 
            linewidth=2.0
        )
        
    # Configure axes with LaTeX formatting
    ax.set_xlabel(r'Radial distance ($\AA$)', fontsize=12)
    ax.set_ylabel(r'$g(r)$', fontsize=12)
    
    # Extract the time step (handles negative numbers)
    file_name = os.path.basename(file_path)
    time_match = re.search(r'PCF_(-?\d+)_fs', file_name)
    time_str = time_match.group(1) if time_match else "Unknown"
    
    # Set the title format
    ax.set_title(f'Radial Distribution Function ({time_str} fs)', fontsize=13, fontweight='bold', pad=15)
    
    # Styling configurations
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.9, fontsize=11)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    
    # Optimize layout and save the file
    plt.tight_layout()
    output_plot_path = file_path.replace('.txt', '.png')
    plt.savefig(output_plot_path, bbox_inches='tight')
    print(f"Successfully generated and saved plot: {output_plot_path}")
    plt.close()

if __name__ == "__main__":
    file_pattern = 'PCF_*_fs.txt'
    files_to_process = sorted(glob.glob(file_pattern))
    
    if not files_to_process:
        print(f"No files matching '{file_pattern}' were found in the current directory.")
    else:
        print(f"Found {len(files_to_process)} file(s) to process.")
        for file_path in files_to_process:
            plot_single_rdf(file_path)