import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path

def plot_contours(ax, gorundTruthPath, contourPointsPath, contourPath, centroidsPath):
    # Read the CSV file
    datagt = pd.read_csv(gorundTruthPath)
    contourPointsPath1= contourPointsPath +"1.csv"
    contourPointsPath2= contourPointsPath +"2.csv"
    contourPath1 = contourPath+"1.csv"
    contourPath2 = contourPath+"2.csv"
    contour_points_1 = pd.read_csv(contourPointsPath1)
    data1 = pd.read_csv(contourPath1)
    centroids = pd.read_csv(centroidsPath, header=None)

    # Assuming the columns in the CSV are named 'X' and 'Y'
    cx1, cy1 = centroids.iloc[0]  

    x1 = datagt['X']
    y1 = datagt['Y']
    gt_points_x_1 =[]
    gt_points_y_1 =[]
    gt_points_x_2 =[]
    gt_points_y_2 =[]
    if(len(x1) == 2000):
        gt_points_x_1 = x1[:1000]
        gt_points_y_1 = y1[:1000]
        gt_points_x_2 = x1[1000:]
        gt_points_y_2 = y1[1000:]
    else:
        gt_points_x_1 = x1
        gt_points_y_1 = y1

    points_x_1 = contour_points_1['X']
    points_y_1 = contour_points_1['Y']

    

    # Plot the points
    ax.plot(points_x_1,points_y_1, "bo", label="Contour points")
    ax.plot(gt_points_x_1, gt_points_y_1,"g-",  label="Ground truth")
    ax.scatter(cx1, cy1,  label="Centroid")
    if(len(x1) == 2000):
        cx2, cy2 = centroids.iloc[1]  
        contour_points_2 = pd.read_csv(contourPointsPath2)
        data2 = pd.read_csv(contourPath2)
        points_x_2 = contour_points_2['X']
        points_y_2 = contour_points_2['Y']
        ax.plot(gt_points_x_2, gt_points_y_2,"g-")
        ax.plot(points_x_2,points_y_2, "bo")
        ax.scatter(cx2, cy2,  label="Centroid")
        xp2 = data2['X']
        yp2 = data2['Y']
        ax.plot(xp2, yp2, "r-")


    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])

    xp1 = data1['X']
    yp1 = data1['Y']

    # Plot the points
    ax.plot(xp1, yp1,"r-", label="Approximated Contour")

    # Display the plot
    ax.set_title('Approximated Tumor Contour')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    ax.legend(loc='lower left', frameon=False, title='Legend')
    
def add_colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar

def plot_heatmap_from_csv(ax, filename, title):
    # Read the CSV file into a pandas DataFrame
    data = pd.read_csv(filename, header=None)
    #data = data[::-1]
    # Plot the heatmap
    plot = ax.imshow(data, cmap='viridis', interpolation='nearest')
    add_colorbar(plot)
    ax.set_xlim(0,data.shape[1])
    ax.set_ylim(0, data.shape[0])
    xticks = np.linspace(0, data.shape[1] - 1, 6)
    yticks = np.linspace(0, data.shape[0] - 1, 6)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_title(title)
    

def plot_scatteredpoints_from_csv(ax, filename):
    # Read the CSV file into a pandas DataFrame
    data = pd.read_csv(filename)
    

    # Extract x and y data
    x = data['x'].values
    y = data['y'].values
    #y = y[::-1]
    # Create a scatter plot
    ax.scatter(x, y, marker='o',s=1, color='blue')
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Sampled Posterior Distribution')
    ax.grid(True)

def main():
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
    groundTruthPath = data_dir+"groundTruth.csv"
    contourPointsPath = data_dir+"contour_points_"
    contourPath=data_dir+"contour_"
    scatterPath = data_dir+"scattered_data.csv"
    heatmapPath = data_dir+"normalized_data.csv"
    groundTruthHeatMapPath = data_dir +"groundTruthHeatMap.csv"
    centroidPath = data_dir + "ms_centroids.csv"
    print(groundTruthHeatMapPath)
    print(groundTruthPath)
    print(contourPointsPath)
    print(contourPath)
    print(scatterPath)
    print(heatmapPath)
    fig, axs = plt.subplots(2,2, figsize=(10, 15))  # Adjust the layout as needed

    # Plotting functions
    plot_heatmap_from_csv(axs[0,0], groundTruthHeatMapPath, "Tumor Model")
    plot_contours(axs[1,1], groundTruthPath, contourPointsPath, contourPath,centroidPath)
    plot_heatmap_from_csv(axs[0,1], heatmapPath, "Posterior Distribution")
    plot_scatteredpoints_from_csv(axs[1,0], scatterPath)
    
    # Show the plots
    #plt.tight_layout()
    plt.show()
if __name__=="__main__":
    main()