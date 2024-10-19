import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.animation import FuncAnimation
import pandas as pd
import h5py

from quadtree import QuadTreeNode


def read_cell_matrix(filename, datasetID):
    hf = h5py.File(filename, "r")
    data = hf.get(f"{int(datasetID):05d}")[...].T # transpose since fortran matrices are stored columnwise
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

    # vx = cell_matrix[:, 5+2::5].reshape(-1)[~np.isnan(points_x)][0]
    # vy = cell_matrix[:, 5+3::5].reshape(-1)[~np.isnan(points_y)][0]
    # print(vx, vy, np.sqrt(vx**2+vy**2))
    # print("cell", np.argmax(cell_matrix[:, 5::5].reshape(-1)))

    return [points_x, points_y], patches

def update_plot(num, img_freq, img_offset, filebasename, artists):

    # print((num*img_freq)+img_offset, end=" ")
    # data = read_cell_matrix(f"{filebasename}_{(num*img_freq)+img_offset:05d}.h5")
    data = read_cell_matrix(f"{filebasename}.h5", (num*img_freq)+img_offset)
    points, grid = get_grid_and_particles(data)


    artists[0].set_offsets(np.c_[points[0], points[1]])
    artists[1].set_paths(grid)




# filename_image = "data/test_8.png"
# filebasename = "data/matrix8"
# filename_image = "../data/test_11.png"
filebasename = "../data/test_11"
# filebasename = "../data/empty1x1"


img_freq = 10*4
img_offset = 0#527
num_images = 200//4
# num_images = 1000 - img_offset

root = QuadTreeNode.from_csv(filebasename+".tree")

lines = []
for n in root.leafs():
    if n.contours.size != 0:
        for c in n.contours:
            lines.append([c[:2], c[2:4]])

print(lines)




# img = plt.imread(filename_image)

fig, ax = plt.subplots()
# ax.imshow(img)
contours = ax.add_collection(LineCollection(lines, colors="orange", linewidths=3))

pts = ax.scatter([],[], color="blue", marker=".")
print(type(pts))
patch = ax.add_collection(PatchCollection([], linewidth=1, edgecolor='r', facecolor='none'))
print(type(patch))


ani = FuncAnimation(fig, update_plot, fargs=(img_freq, img_offset, filebasename, [pts, patch, contours]),
                    frames=num_images+1, interval=1, repeat=True)

#
# # ax.set_xlim(-1, 1025)
# # ax.set_ylim(-1, 1025)
# # ax.set_xlim(-1, 4001)
# # ax.set_ylim(-1, 3001)
# # ax.set_xlim(-.1,1.1)
# # ax.set_ylim(-.1,1.1)
ax.set_xlim(-1, root.width+1)
ax.set_ylim(-1, root.height+1)
plt.show()
