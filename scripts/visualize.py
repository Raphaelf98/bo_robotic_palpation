import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
from pathlib import Path
def read_contour_data(file_path):
    contours = []
    current_contour = {}

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("Contour"):
                if current_contour:
                    contours.append(current_contour)
                current_contour = {"Contour": line}
            elif line.startswith("Specificity"):
                _, value = line.split(": ")
                current_contour["Specificity"] = float(value)
            elif line.startswith("Sensitivity"):
                _, value = line.split(": ")
                current_contour["Sensitivity"] = float(value)

        # Add the last contour
        if current_contour:
            contours.append(current_contour)

    return contours
def plot_contours(ax, gorundTruthPath, contourPointsPath, contourPath, centroidsPath, metricsPath):
    # Read the CSV file
    
    contourPointsPath1= contourPointsPath +"1.csv"
    contourPointsPath2= contourPointsPath +"2.csv"
    contourPath1 = contourPath+"1.csv"
    contourPath2 = contourPath+"2.csv"
    contour_points_1 = pd.read_csv(contourPointsPath1)
    data1 = pd.read_csv(contourPath1)
    centroids = pd.read_csv(centroidsPath, header=None)
    # get metrics
    datagt = None
    contour_data = None
     # Check if metricsPath exists
    if os.path.exists(metricsPath):
        contour_data = read_contour_data(metricsPath)
        datagt = pd.read_csv(gorundTruthPath)
    else:
        # Set default values or handle the absence of the file
        contour_data = [{} for _ in range(len(centroids))]  # Assuming we want one empty dict for each centroid
        datagt = {'X': [], 'Y': []}

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
    ax.plot(points_x_1,points_y_1, "bo", label="Contour Points")
    ax.plot(gt_points_x_1, gt_points_y_1,"g-",  label="Ground Truth")
    ax.scatter(cx1, cy1,  label="Centroid")

    # Plot lines connecting each point in contour_points_1 to the centroid
    for x, y in zip(contour_points_1['X'], contour_points_1['Y']):
        ax.plot([cx1, x], [cy1, y], 'k--', linewidth=0.5)  # 'k-' sets the color to black

    # Check if you have another set of points and centroid to plot
    if len(centroids) > 1:
        # Similar code for the second set of contour points and centroid
        cx2, cy2 = centroids.iloc[1]
        contour_points_2 = pd.read_csv(contourPointsPath2)
        for x, y in zip(contour_points_2['X'], contour_points_2['Y']):
            ax.plot([cx2, x], [cy2, y], 'k--', linewidth=0.5)  # 'k-' sets the color to black

        

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
    for i, contour in enumerate(contour_data, start=1):
        spec = contour.get('Specificity', 'N/A')
        sens = contour.get('Sensitivity', 'N/A')

        # Choose x and y for text placement. Adjust these as necessary for your plot.
        x_text = 0.1  # For example, 5% from the left
        y_text = 1 - 0.04 * i  # For example, starting from 95% from the bottom and going up every iteration

        plt.text(x_text, y_text, f"Centroid [{round(centroids.iloc[i-1][0],2)}, {round(centroids.iloc[i-1][1],2)}] - Specificity: {spec}, Sensitivity: {sens}", 
                 transform=plt.gca().transAxes, fontsize=7, verticalalignment='top')
    ax.legend(loc='lower left', frameon=False)
    
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
    plot = ax.imshow(data, cmap='jet', interpolation='nearest')
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
    result_dir = str(parent_directory)+"/data/"+args.input_string+"/results/"
    groundTruthPath = data_dir+"groundTruth.csv"
    contourPointsPath = data_dir+"contour_points_"
    contourPath=data_dir+"contour_"
    scatterPath = data_dir+"scattered_data.csv"
    heatmapPath = data_dir+"posterior.txt"
    groundTruthHeatMapPath = data_dir +"groundTruthHeatMap.csv"
    centroidPath = data_dir + "ms_centroids.csv"
    metricsPath = result_dir+"metrics.txt"
    print(groundTruthHeatMapPath)
    print(groundTruthPath)
    print(contourPointsPath)
    print(contourPath)
    print(scatterPath)
    print(heatmapPath)
    print(metricsPath)
    fig, axs = plt.subplots(2,2, figsize=(10, 15))  # Adjust the layout as needed

    # Plotting functions
    plot_heatmap_from_csv(axs[0,0], groundTruthHeatMapPath, "Tumor Model")
    plot_contours(axs[1,1], groundTruthPath, contourPointsPath, contourPath,centroidPath, metricsPath)
    plot_heatmap_from_csv(axs[0,1], heatmapPath, "Posterior Distribution")
    plot_scatteredpoints_from_csv(axs[1,0], scatterPath)
    
    # Show the plots
    #plt.tight_layout()
    plt.show()
if __name__=="__main__":
    main()