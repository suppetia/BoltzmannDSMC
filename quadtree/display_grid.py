import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import pandas as pd
import h5py

def read_cell_matrix(filename):
    # return np.loadtxt(filename, delimiter=" ", dtype=np.float64)
    if filename.lower().split(".")[1] == "txt":
        data = pd.read_csv(filename, delimiter="\s+", header=None).to_numpy()
    elif filename.lower().split(".")[1] == "h5":
        # df = pd.read_hdf(filename, key="gridCells")
        hf = h5py.File(filename, "r")
        data = hf.get("gridCells")[...].T # transpose since fortran matrices are stored columnwise
    else:
        raise ValueError("invalid file format")
    return data

def display_cells(cell_matrix, canvas=None):
    if canvas is None:
        fig, ax = plt.subplots()
    else:
        fig, ax = canvas

    patches = []
    points = []

    for cell in cell_matrix:
        x, y, width, height = cell[:4]
        # Create a rectangle patch
        rectangle = Rectangle((x, y), width, height)
        # Add the rectangle to the plot
        patches.append(rectangle)
        # # ax.add_patch(rectangle)
        # for i in range(5, len(cell), 4):
        #     if cell[i] < 0:
        #         break
        #     points.append([cell[i], cell[i+1]])
        #     print(cell[i], cell[i+1])

    points_x = cell_matrix[:, 5::4].reshape(-1)
    np.place(points_x, points_x < 0, np.nan)
    points_y = cell_matrix[:, 6::4].reshape(-1)
    np.place(points_y, points_y < 0, np.nan)

    # points = np.array(points)
    # ax.scatter(points[:,0], points[:,1], color="blue", marker="o")
    ax.scatter(points_x, points_y, color="blue", marker="o")
    ax.add_collection(PatchCollection(patches, linewidth=1, edgecolor='r', facecolor='none'))
    return fig, ax

filename_image = "data/circle_30_60.png"
filename_cells = "matrix7_cells.txt"
filename_image = "data/test_8.png"
filename_cells = "matrix8_cells.txt"
filename_cells = "matrix8_cells.h5"

cm = read_cell_matrix(filename_cells)
img = plt.imread(filename_image)

fig, ax = plt.subplots()
ax.imshow(img)

display_cells(cm, (fig, ax))

# ax.set_xlim(-1, 1025)
# ax.set_ylim(-1, 1025)
plt.show()
