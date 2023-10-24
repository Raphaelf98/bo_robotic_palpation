import pandas as pd
import matplotlib.pyplot as plt

def plot_heatmap_from_csv(filename):
    # Read the CSV file into a pandas DataFrame
    data = pd.read_csv(filename, header=None)

    # Plot the heatmap
    plt.imshow(data, cmap='hot', interpolation='nearest')
    plt.colorbar()
    plt.title("Heatmap from CSV")
    plt.show()

if __name__ == "__main__":
    # Replace 'input.csv' with the name of your CSV file
    plot_heatmap_from_csv("../config/normalized_data.csv")