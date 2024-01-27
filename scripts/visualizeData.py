import pandas as pd
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
def plot_points_from_csv(filename):
    # Read the CSV file into a pandas DataFrame
    data = pd.read_csv(filename)
    data = data[::-1]

    # Extract x and y data
    x = data['x'].values
    y = data['y'].values

    # Create a scatter plot
    plt.scatter(x, y, marker='o', color='blue')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Scatter plot of points from CSV')
    plt.grid(True)
    
    plt.show()

if __name__ == "__main__":
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
    plot_points_from_csv(data_dir+"scattered_data.csv")
    #plot_points_from_csv("../config/clusterCenters.csv")