import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the CSV file
datagt = pd.read_csv('/home/raphael/robolab/build/groundTruth.csv')

contour_points_1 = pd.read_csv('/home/raphael/robolab/displaygp/config/contour_points_1.csv')
data1 = pd.read_csv('/home/raphael/robolab/displaygp/config/contour_1.csv')
contour_points_2 = pd.read_csv('/home/raphael/robolab/displaygp/config/contour_points_2.csv')
data2 = pd.read_csv('/home/raphael/robolab/displaygp/config/contour_2.csv')

# Assuming the columns in the CSV are named 'X' and 'Y'
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

points_x_2 = contour_points_2['X']
points_y_2 = contour_points_2['Y']

# Plot the points
plt.plot(points_x_1,points_y_1, "bo", label="Contour points")
plt.plot(gt_points_x_1, gt_points_y_1,"g-",  label="Ground truth")

if(len(x1) == 2000):
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
plt.plot(xp1, yp1,"r-", label="Approximated Contour")

# Display the plot
plt.title('Points from contour.csv')
plt.xlabel('X')
plt.ylabel('Y')

plt.legend(loc='lower left', frameon=False, title='Legend')
plt.show()