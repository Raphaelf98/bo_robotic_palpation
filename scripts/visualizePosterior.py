import pandas as pd
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
import numpy as np
import matplotlib
from matplotlib.colors import Normalize
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
def plot_heatmap_from_csv(filename,file_samples):
    # Read the CSV file into a pandas DataFrame
    data = pd.read_csv(filename, header=None)
    samples = pd.read_csv(file_samples)
    samples_x = samples['x']
    samples_y = samples['y']
    samples_x = 100*samples_x
    samples_y = 100*samples_y
        # Create a figure and an axes
    fig, ax = plt.subplots()


    # Plot the heatmap using pcolormesh
    # Plot the heatmap
    cax = ax.imshow(data, cmap='jet', interpolation='nearest')

    ax.scatter(samples_x, samples_y, s=50, c='black', edgecolor='black', linewidth=1, marker='+', label='Samples')

    fig.colorbar(cax)
    ax.set_xlim(0,data.shape[1])
    ax.set_ylim(0, data.shape[0])
    xticks = np.linspace(0, data.shape[1] - 1, 6)
    yticks = np.linspace(0, data.shape[0] - 1, 6)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_title("Posterior")
    plt.show()
    return fig

def plot_3d_from_csv(filename):
    # Read the CSV file into a pandas DataFrame
    data = pd.read_csv(filename, header=None)
    data = -data # Optional based on how you want to view the plot

    # Create meshgrid
    x = np.arange(data.shape[1])
    y = np.arange(data.shape[0])
    X, Y = np.meshgrid(x, y)
    Z = data.values

    # Create a figure and a 3D Axes
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface
    surf = ax.plot_surface(X, Y, Z, cmap='jet')

    # Add a color bar which maps values to colors.
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7]) 
    fig.colorbar(surf, cax=cbar_ax)
       # Set Z-axis limits
    ax.set_zlim(0, 2)

    # Remove ticks
    #ax.set_xticks([])
    #ax.set_yticks([])
    #ax.set_zticks([])
    # Set labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Stiffness')

    plt.show()

if __name__ == "__main__":
    # Replace 'input.csv' with the name of your CSV file
     # Replace 'filename.csv' with the name of your CSV file
    parser = argparse.ArgumentParser(description="Process a string.")
    parser.add_argument('input_string', type=str, help='A string input by the user')

    args = parser.parse_args()

    print(f"Running Visualization on Experiment : {args.input_string}")

    current_directory = Path(__file__).resolve()
    parent_directory = current_directory.parent
    parent_directory = parent_directory.parent
    print(f"Current Directory: {current_directory}")
    print(f"Parent Directory: {parent_directory}")
    data_dir = str(parent_directory)+"/data/"+args.input_string+"/log/"
    fig = plot_heatmap_from_csv(data_dir+"posterior.txt", data_dir+"gp_samples.csv")
    fig.savefig(f'/home/raphael/Desktop/pacs-project/Report/Experiments/{args.input_string}_posterior_plot.pgf')

    plot_3d_from_csv(data_dir+"posterior.txt")