import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd

def read_cell_matrix(filename):
    # return np.loadtxt(filename, delimiter=" ", dtype=np.float64)
    df = pd.read_csv(filename, delimiter="\s+", header=None)
    return df.to_numpy()

def display_cells(cell_matrix, canvas=None):
    if canvas is None:
        fig, ax = plt.subplots()
    else:
        fig, ax = canvas

    for cell in cell_matrix:
        x, y, width, height = cell[:4]
        # Create a rectangle patch
        rectangle = Rectangle((x, y), width, height, linewidth=1, edgecolor='r', facecolor='none')
        # Add the rectangle to the plot
        ax.add_patch(rectangle)
        for i in range(4, len(cell), 4):
            if cell[i] < 0:
                break
            ax.scatter(cell[i], cell[i+1], color="blue", marker="o")
            print(cell[i], cell[i+1])

    return fig, ax

filename_image = "data/circle_30_60.png"
filename_cells = "matrix7_cells.txt"
filename_image = "data/test_8.png"
filename_cells = "matrix8_cells.txt"

cm = read_cell_matrix(filename_cells)
img = plt.imread(filename_image)

fig, ax = plt.subplots()
ax.imshow(img)

display_cells(cm, (fig, ax))

# ax.set_xlim(-1, 1025)
# ax.set_ylim(-1, 1025)
plt.show()
