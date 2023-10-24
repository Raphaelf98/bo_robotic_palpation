import pandas as pd
import matplotlib.pyplot as plt

def plot_points_from_csv(filename):
    # Read the CSV file into a pandas DataFrame
    data = pd.read_csv(filename)

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
    #plot_points_from_csv("../config/data.csv")
    plot_points_from_csv("../config/clusterCenters.csv")