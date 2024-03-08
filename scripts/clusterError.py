
# import os
# import pandas as pd
# import math
# import matplotlib.pyplot as plt
# import numpy as np
# from scipy.optimize import curve_fit
# from pathlib import Path
# def func(x, a, b, c):
#     return a * np.exp(-b * x) + c
# def read_configuration(file_path):
#     """Reads the configuration file to extract n_iterations value."""
#     n_iterations = None
#     with open(file_path, 'r') as file:
#         for line in file:
#             if line.startswith('n_iterations='):
#                 n_iterations = int(line.strip().split('=')[1])
#                 break
#     return n_iterations

# def read_coordinates(csv_file_path):
#     """Reads the CSV file to extract x and y coordinates."""
#     custom_headers = ['x', 'y']
#     data = pd.read_csv(csv_file_path, header=None, names=custom_headers)
#     return data['x'].values[0], data['y'].values[0]

# def compute_euclidean_distances(coordinates_x, coordinates_y, gt_x=0.5, gt_y=0.5):
#     """Computes Euclidean distances of coordinates from a ground truth point."""
#     distances = [math.sqrt((gt_x - x) ** 2 + (gt_y - y) ** 2) for x, y in zip(coordinates_x, coordinates_y)]
#     return distances

# def fit_curve(iters, distances):
#     """Fits a curve to the given iterations and distances using a predefined function."""
#     def func(x, a, b, c):
#         return a * np.exp(-b * x) + c
    
#     params, covariance = curve_fit(func, iters, distances, bounds=(0, [3., 1., 0.5]))
#     return params

# def plot_data_and_model(iters, distances, params, label, color):
#     """Plots the original data points and the fitted model curve."""
#     x_line = np.linspace(min(iters), max(iters), 100)
#     y_line = func(x_line, *params)
#     plt.plot(iters, distances, 'o', color=color)
#     plt.plot(x_line, y_line, color=color, label=label)

# def process_directory(parent_directory, shape, gt_x, gt_y):
#     """Processes directories of a given shape to extract iterations, coordinates, and fit a model."""
#     coordinates_x, coordinates_y, iters = [], [], []
#     for item in os.listdir(parent_directory):
#         item_path = os.path.join(parent_directory, item)
#         if os.path.isdir(item_path) and item.startswith(shape):
#             conf_file_path = os.path.join(item_path, "parameters", "bo_parameters.txt")
#             csv_file_path = os.path.join(item_path, "log", "ms_centroids.csv")
            
#             n_iterations = read_configuration(conf_file_path)
#             if n_iterations is not None:
#                 iters.append(n_iterations)
            
#             if os.path.exists(csv_file_path):
#                 x, y = read_coordinates(csv_file_path)
#                 coordinates_x.append(x)
#                 coordinates_y.append(y)
    
#     distances = compute_euclidean_distances(coordinates_x, coordinates_y, gt_x, gt_y)
#     params = fit_curve(iters, distances)
#     return iters, distances, params

# # Main code
# current_directory = Path(__file__).resolve()
# parent_directory = current_directory.parent
# parent_directory = parent_directory.parent
# print(f"Current Directory: {current_directory}")
# print(f"Parent Directory: {parent_directory}")
# data_dir = str(parent_directory)+"/data/"

# shapes = ['Circle', 'Triangle', 'Rectangle']
# colors = ['red', 'blue', 'green']
# labels = ['Circle Model', 'Triangle Model', 'Rectangle Model']
# # Define ground truth coordinates for each shape
# ground_truths = {
#     'Circle': (0.5, 0.5),
#     'Triangle': (0.5, 0.5),  # Modify these values as per your requirements
#     'Rectangle': (0.5, 0.5)  # Modify these values as per your requirements
# }
# plt.figure(figsize=(10, 6))
# for shape, color, label in zip(shapes, colors, labels):
#     gt_x, gt_y = ground_truths[shape]
#     iters, distances, params = process_directory(data_dir, shape,gt_x,gt_y)
#     plot_data_and_model(iters, distances, params, label, color)

# plt.xlabel('Iterations N')
# plt.ylabel('Euclidean Error d')
# plt.legend(loc='upper right', frameon=False)
# plt.show()
# Import necessary libraries
import os
import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from pathlib import Path
import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
# Define the model function for curve fitting
def func(x, a, b, c):
    """Exponential decay function used for curve fitting."""
    return a * np.exp(-b * x) + c

# Function to read the configuration file and extract the number of iterations
def read_configuration(file_path):
    """Reads the configuration file to find and return the 'n_iterations' value."""
    n_iterations = None
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('n_iterations='):
                n_iterations = int(line.strip().split('=')[1])
                break
    return n_iterations

# Function to read x and y coordinates from a CSV file
def read_coordinates(csv_file_path):
    """Extracts x and y coordinates from the specified CSV file."""
    custom_headers = ['x', 'y']
    data = pd.read_csv(csv_file_path, header=None, names=custom_headers)
    return data['x'].values[0], data['y'].values[0]

# Function to compute Euclidean distances from each coordinate to a specified ground truth point
def compute_euclidean_distances(coordinates_x, coordinates_y, gt_x, gt_y):
    """Calculates the Euclidean distances between given coordinates and a ground truth point."""
    distances = [math.sqrt((gt_x - x) ** 2 + (gt_y - y) ** 2) for x, y in zip(coordinates_x, coordinates_y)]
    return distances

# Function to fit a curve using the specified function and bounds
def fit_curve(iters, distances):
    """Fits an exponential decay curve to the provided data."""
    params, covariance = curve_fit(func, iters, distances, bounds=(0, [3., 1., 0.5]))
    return params

# Function to plot the original data and the fitted model
def plot_data_and_model(iters, distances, params, label, color):
    """Creates a plot of the original data points and the curve fitted to them."""
    x_line = np.linspace(min(iters), max(iters), 100)
    y_line = func(x_line, *params)
    plt.plot(iters, distances, 'o', color=color)
    plt.plot(x_line, y_line, color=color, label=label)

# Function to process directories and aggregate data for modeling and plotting
def process_directory(parent_directory, shape, gt_x, gt_y):
    """Processes directories of a given shape to extract data, compute distances, and fit a curve."""
    coordinates_x, coordinates_y, iters = [], [], []
    for item in os.listdir(parent_directory):
        item_path = os.path.join(parent_directory, item)
        if os.path.isdir(item_path) and item.startswith(shape):
            conf_file_path = os.path.join(item_path, "parameters", "bo_parameters.txt")
            csv_file_path = os.path.join(item_path, "log", "ms_centroids.csv")
            
            n_iterations = read_configuration(conf_file_path)
            if n_iterations is not None:
                iters.append(n_iterations)
            
            if os.path.exists(csv_file_path):
                x, y = read_coordinates(csv_file_path)
                coordinates_x.append(x)
                coordinates_y.append(y)
    
    distances = compute_euclidean_distances(coordinates_x, coordinates_y, gt_x, gt_y)
    params = fit_curve(iters, distances)
    return iters, distances, params

# Main script execution starts here
# Resolves the current and parent directory paths for data file access
current_directory = Path(__file__).resolve()
parent_directory = current_directory.parent.parent
data_dir = str(parent_directory) + "/data/"

# Defines shapes, colors, and labels for the experiments
shapes = ['Circle', 'Triangle', 'Rectangle']
colors = ['red', 'blue', 'green']
labels = ['Circle Model', 'Triangle Model', 'Rectangle Model']

# Defines ground truth coordinates for each shape
ground_truths = {
    'Circle': (0.3, 0.7),
    'Triangle': (0.3, 0.6),
    'Rectangle': (0.4, 0.4)
}

# Sets up the plot
plt.figure(figsize=(10, 6))

# Processes each shape directory, fits curves to the data, and plots the results
for shape, color, label in zip(shapes, colors, labels):
    gt_x, gt_y = ground_truths[shape]
    iters, distances, params = process_directory(data_dir, shape, gt_x, gt_y)
    plot_data_and_model(iters, distances, params, label, color)

# Configures plot aesthetics and displays the plot
plt.xlabel('Iterations N')
plt.ylabel('Euclidean Error d')
plt.legend(loc='upper right', frameon=False)
plt.show()
plt.savefig(f'/home/raphael/Desktop/pacs-project/Report/Experiments/euc_error_iter_plot.pgf')
