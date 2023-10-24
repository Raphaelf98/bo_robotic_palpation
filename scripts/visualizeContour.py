import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
data1 = pd.read_csv('/home/raphael/robolab/displaygp/config/contour_1.csv')
data2 = pd.read_csv('/home/raphael/robolab/displaygp/config/contour_2.csv')

# Assuming the columns in the CSV are named 'X' and 'Y'
x1 = data1['X']
y1 = data1['Y']

# Plot the points
plt.plot(x1, y1)
plt.xlim(0, 1)
plt.ylim(0, 1)
x2 = data2['X']
y2 = data2['Y']

# Plot the points
plt.plot(x2, y2)


# Display the plot
plt.title('Points from contour.csv')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()