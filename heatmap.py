import numpy as np
import matplotlib.pyplot as plt
import csv

DCOL = 2  # data col

def hmap(fname, studytype=None, legend=None):
    # Read the CSV file
    data = []
    if legend is None:
        legend = 'Response'
    if studytype is  None:
        pass
    with open(fname, 'r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            data.append(row)

    # Convert string values to float where needed
    for i in range(len(data)):
        data[i][0] = float(data[i][0])  # col 0
        data[i][1] = float(data[i][1])  # col 1
        data[i][DCOL] = float(data[i][DCOL])  # col 3

    # Get the unique values for x and y axes
    x_values = [row[0] for row in data]
    y_values = [row[1] for row in data]
    x_unique = sorted(list(set(x_values)))
    y_unique = sorted(list(set(y_values)))

    # Create a 2D array for the heatmap data
    heatmap_data = np.full((len(y_unique), len(x_unique)), np.nan)

    # Fill the heatmap data
    for row in data:
        x_val = row[0]
        y_val = row[1]
        color_val = row[DCOL]  # col 3 for color

        x_idx = x_unique.index(x_val)
        y_idx = y_unique.index(y_val)

        heatmap_data[y_idx, x_idx] = color_val

    # Create the heatmap
    plt.figure(figsize=(8, 6))
    plt.imshow(heatmap_data, cmap='viridis', aspect='auto')

    # Add colorbar
    plt.colorbar(label=legend)

    # Set axis labels and ticks
    plt.xticks(range(len(x_unique)), x_unique)
    plt.yticks(range(len(y_unique)), y_unique)
    plt.xlabel('Omega_n')
    plt.ylabel('Zeta')

    # Add title
    # plt.title(legend+' vs. dynamic parameters.')
    plt.title('Simulation: '+legend)

    plt.tight_layout()
    plt.show()

if __name__=='__main__':
    hmap('data.csv' )

