import pandas as pd
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
import numpy as np

def plot_heatmap_from_csv(filename):
    # Read the CSV file into a pandas DataFrame
    data = pd.read_csv(filename, header=None)
    data = data[::-1]
    print(data.shape)
       # Create a figure and an axes
    fig, ax = plt.subplots()

    # Plot the heatmap
    cax = ax.imshow(data, cmap='viridis', interpolation='nearest')
    fig.colorbar(cax)
    ax.set_xlim(0,data.shape[1])
    ax.set_ylim(0, data.shape[0])
    xticks = np.linspace(0, data.shape[1] - 1, 6)
    yticks = np.linspace(0, data.shape[0] - 1, 6)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_title("Ground Truth")
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
    plot_heatmap_from_csv(data_dir+"groundTruthHeatMap.csv")