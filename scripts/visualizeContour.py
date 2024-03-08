import pandas as pd
import matplotlib.pyplot as plt

import numpy as np
import argparse
import os
from pathlib import Path
import matplotlib
# matplotlib.use("pgf")
# matplotlib.rcParams.update({
#     "pgf.texsystem": "pdflatex",
#     'font.family': 'serif',
#     'text.usetex': True,
#     'pgf.rcfonts': False,
# })
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


# Read the CSV file
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


datagt = pd.read_csv(data_dir+'groundTruth.csv')
centroidsPath = data_dir + "ms_centroids.csv"

centroids = pd.read_csv(centroidsPath, header=None)
contour_points_1 = pd.read_csv(data_dir+'contour_points_1.csv')
data1 = pd.read_csv(data_dir+'contour_1.csv')
contour_points_2 =[]
data2 = []
result_file   = result_dir+"metrics.txt"
contour_data = None

 # Check if metricsPath exists
if os.path.exists(result_file):
    contour_data = read_contour_data(result_file)
else:
    # Set default values or handle the absence of the file
    contour_data = [{} for _ in range(len(centroids))]  # Assuming we want one empty dict for each centroid
cx1, cy1 = centroids.iloc[0]
# Assuming the columns in the CSV are named 'X' and 'Y'
x1 = datagt['X']
y1 = datagt['Y']
gt_points_x_1 =[]
gt_points_y_1 =[]
gt_points_x_2 =[]
gt_points_y_2 =[]
plt.scatter(cx1, cy1,  label="Centroid 1")
# Plot lines connecting each point in contour_points_1 to the centroid
for x, y in zip(contour_points_1['X'], contour_points_1['Y']):
    plt.plot([cx1, x], [cy1, y], 'k--', linewidth=0.5)  # 'k-' sets the color to black

    
if(len(x1) == 2000):
    contour_points_2 = pd.read_csv(data_dir+'contour_points_2.csv')
    # Similar code for the second set of contour points and centroid
    cx2, cy2 = centroids.iloc[1]
    for x, y in zip(contour_points_2['X'], contour_points_2['Y']):
        plt.plot([cx2, x], [cy2, y], 'k--', linewidth=0.5)  # 'k-' sets the color to black
    data2 = pd.read_csv(data_dir+'contour_2.csv')
    gt_points_x_1 = x1[:1000]
    gt_points_y_1 = y1[:1000]
    gt_points_x_2 = x1[1000:]
    gt_points_y_2 = y1[1000:]
    cx2, cy2 = centroids.iloc[1]  
    plt.scatter(cx2, cy2,  label="Centroid 2")
else:
    gt_points_x_1 = x1
    gt_points_y_1 = y1

points_x_1 = contour_points_1['X']
points_y_1 = contour_points_1['Y']


# Plot the points
plt.plot(points_x_1,points_y_1, "bo", label="Contour Points")
plt.plot(gt_points_x_1, gt_points_y_1,"g-",  label="Ground Truth")
#plot specificity and sensitivity


if(len(x1) == 2000):
    points_x_2 = contour_points_2['X']
    points_y_2 = contour_points_2['Y']
    plt.plot(gt_points_x_2, gt_points_y_2,"g-")
    plt.plot(points_x_2,points_y_2, "bo")

    xp2 = data2['X']
    yp2 = data2['Y']
    plt.plot(xp2, yp2, "r-")


plt.xlim(0, 1)
plt.ylim(0, 1)

xp1 = data1['X']
yp1 = data1['Y']

# Plot the points
plt.plot(xp1, yp1,"r-", label="Parametric Cubic Spline Curve")

# Display the plot
plt.title('')
plt.xlabel('X')
plt.ylabel('Y')

for i, contour in enumerate(contour_data, start=1):
    spec = contour.get('Specificity', 'N/A')
    sens = contour.get('Sensitivity', 'N/A')
    
    # Choose x and y for text placement. Adjust these as necessary for your plot.
    x_text = 0.0  # For example, 5% from the left
    y_text = 1 - 0.04 * i  # For example, starting from 95% from the bottom and going up every iteration

    plt.text(x_text, y_text, f"Contour {i} - Specificity: {spec}, Sensitivity: {sens}", 
             transform=plt.gca().transAxes, fontsize=7, verticalalignment='top')


plt.legend(loc='lower left', frameon=False)
# plt.savefig(f'/home/raphael/Desktop/pacs-project/Report/Experiments/{args.input_string}_contour_plot.pgf')
plt.show()
