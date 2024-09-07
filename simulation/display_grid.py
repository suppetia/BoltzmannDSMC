import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.animation import FuncAnimation
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

def get_grid_and_particles(cell_matrix):
    patches = []

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

    points_x = cell_matrix[:, 5::5].reshape(-1)
    np.place(points_x, points_x < 0, np.nan)
    points_y = cell_matrix[:, 6::5].reshape(-1)
    np.place(points_y, points_y < 0, np.nan)

    return [points_x, points_y], patches

def update_plot(num, filebasename, artists):

    data = read_cell_matrix(f"{filebasename}_{num*10:05d}.h5")
    points, grid = get_grid_and_particles(data)

    artists[0].set_offsets(np.c_[points[0], points[1]])
    artists[1].set_paths(grid)




filename_image = "data/test_8.png"
filebasename = "data/matrix8"

num_images = 250



# img = plt.imread(filename_image)

fig, ax = plt.subplots()
# ax.imshow(img)

pts = ax.scatter([],[], color="blue", marker=".")
print(type(pts))
patch = ax.add_collection(PatchCollection([], linewidth=1, edgecolor='r', facecolor='none'))
print(type(patch))


ani = FuncAnimation(fig, update_plot, fargs=(filebasename, [pts, patch]),
                    frames=num_images+1, interval=100, repeat=True)


# ax.set_xlim(-1, 1025)
# ax.set_ylim(-1, 1025)
ax.set_xlim(-.1,1.1)
ax.set_ylim(-.1,1.1)
plt.show()
