import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def fx(t):
    return 0.5 + 0.15*np.sin(t)
def fy(t):
    return 0.5 + 0.15*np.cos(t)


# Read the CSV file
data1 = pd.read_csv('/home/raphael/robolab/displaygp/config/contour_1.csv')
contour_points = pd.read_csv('/home/raphael/robolab/displaygp/config/contour_points_1.csv')
#data2 = pd.read_csv('/home/raphael/robolab/displaygp/config/contour_2.csv')

# Assuming the columns in the CSV are named 'X' and 'Y'
x1 = data1['X']
y1 = data1['Y']

points_x = contour_points['X']
points_y = contour_points['Y']

# Plot the points
plt.plot(points_x,points_y, "ro")
plt.plot(x1, y1)
plt.xlim(0, 1)
plt.ylim(0, 1)
#x2 = data2['X']
#y2 = data2['Y']
x = []
y = []
t = np.linspace(0, 2*np.pi, 1000)

plt.plot(fx(t), fy(t))


# Plot the points
#plt.plot(x2, y2)


# Display the plot
plt.title('Points from contour.csv')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()