import pandas as pd
import matplotlib.pyplot as plt

def plot_csv(filename):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(filename)

    # Check that there are exactly 2 columns
    if len(df.columns) != 2:
        print("The CSV file should have exactly two columns.")
        return

    # Plot the data
    plt.plot(df[df.columns[0]], df[df.columns[1]])
    plt.xlabel(df.columns[0])
    plt.ylabel(df.columns[1])
    plt.title("Plot of CSV Data")
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    filename = "integral_vals.csv"  # replace with your CSV filename
    plot_csv(filename)