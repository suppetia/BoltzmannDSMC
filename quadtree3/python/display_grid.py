import os
# os.environ["QT_QPA_PLATFORM"] = "wayland"
# os.environ["XDG_SESSION_TYPE"] = "xcb"

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.animation import FuncAnimation
import matplotlib.colors as mcolors
import pandas as pd
import h5py
from scipy.optimize import curve_fit

from quadtree import QuadTreeNode

cmap = plt.get_cmap('inferno')
max_brightness = 1e-30


def boltzmann_total(v, a, b):
    return a * v**2 * np.exp(-b*v**2)
def boltzmann(v,a,b):
    return a * np.exp(-b*v**2)


def read_cell_matrix(filename, datasetID):
    hf = h5py.File(filename, "r")
    data = np.array(hf.get(f"{int(datasetID):05d}_tree")) # transpose since fortran matrices are stored columnwise
    return data

def read_particle_matrix(filename, datasetID):
    hf = h5py.File(filename, "r")
    data = hf.get(f"{int(datasetID):05d}_particles")
    if data is None:
        return np.array([[]])
    return np.array(data) # transpose since fortran matrices are stored columnwise

def get_grid_and_particles(cell_matrix, particle_matrix):
    patches = []
    colors = []

    global max_brightness
    # max_brightness = max_brightness[1:] + [np.max(cell_matrix[:, 5])]#
    # print(cell_matrix.shape)
    # print(np.max(cell_matrix[:, 4], axis=0), np.argmax(cell_matrix[:,4], axis=0))
    # print(cell_matrix[:, :2])
    max_brightness = max(max_brightness, np.max(cell_matrix[:, 4], axis=0))
    # print(max_brightness)
    # cnorm = mcolors.Normalize(vmin=0, vmax=np.mean(max_brightness), clip=True)
    cnorm = mcolors.Normalize(vmin=1e-3*max_brightness, vmax=.1*max_brightness, clip=True)
    # cnorm = mcolors.LogNorm(vmin=1e-3*max_brightness, vmax=.1*max_brightness, clip=True)

    for cell in cell_matrix:
        # print(cell[4])
        x, y, width, height = cell[:4]
        # Create a rectangle patch
        rectangle = Rectangle((x, y), width, height)
        # Add the rectangle to the plot
        patches.append(rectangle)

        colors.append(cmap(cnorm(cell[4])))

        # # ax.add_patch(rectangle)
        # for i in range(5, len(cell), 4):
        #     if cell[i] < 0:
        #         break
        #     points.append([cell[i], cell[i+1]])
        #     print(cell[i], cell[i+1])

    points_x = particle_matrix[0,:]
    # np.place(points_x, points_x < 0, np.nan)
    vel_x = particle_matrix[2,:]
    # np.place(vel_x, vel_x < 0, np.nan)
    points_y = particle_matrix[1,:]
    # np.place(points_y, points_y < 0, np.nan)
    vel_y = particle_matrix[3,:]
    # np.place(vel_y, vel_y < 0, np.nan)
    vel_z = particle_matrix[4,:]
    # np.place(vel_z, vel_z < 0, np.nan)

    return [points_x, points_y, vel_x, vel_y, vel_z], (patches, colors)

def update_plot(num, img_freq, img_offset, filebasename, artists, display_params, hist_params):

    # print(num)

    disp_grid, disp_particles, disp_density, disp_vel_hist = display_params

    rect, ax_hist, num_bins = hist_params

    # print((num*img_freq)+img_offset, end=" ")
    # data = read_cell_matrix(f"{filebasename}_{(num*img_freq)+img_offset:05d}.h5")
    cell_mat = read_cell_matrix(f"{filebasename}.h5", (num*img_freq)+img_offset)
    # print(cell_mat)
    particle_mat = read_particle_matrix(f"{filebasename}.h5", (num*img_freq)+img_offset)
    # print(particle_mat)
    points, (grid, colors) = get_grid_and_particles(cell_mat, particle_mat)

    if disp_vel_hist:
        # find points in a given rectangle
        mask_x = (points[0] >= rect[0]) & (points[0] < rect[0] + rect[2])
        mask_y = (points[1] >= rect[1]) & (points[1] < rect[1] + rect[3])
        mask = mask_x & mask_y

        vel = np.sqrt(points[2][mask]**2 + points[3][mask]**2 + points[4][mask]**2)
        if vel.size > 0:
            mean_vel = np.mean(vel)

            for dim in ["x", "y"]:
                if dim == "x":
                    idx = 2
                elif dim == "y":
                    idx = 3
                ax_hist[idx-2].cla()

                v = points[idx][mask]
                v_mean = np.mean(v)
                counts, bins, _ = ax_hist[idx-2].hist(v-v_mean, bins=num_bins, density=True, alpha=.5)
                # counts, bins, _ = ax_hist.hist(vel, bins=num_bins, density=True)
                bin_centers = 0.5 * (bins[1:] + bins[:-1])

                try:
                    # Fit the histogram data to the Boltzmann distribution
                    popt, pcov = curve_fit(boltzmann, bin_centers, counts, p0=[1/(counts.max()*bin_centers[len(bin_centers)//2]), 1e-13])  # Initial guesses for a and b
                    x_fit = np.linspace(bin_centers.min(), bin_centers.max(), 100)
                    y_fit = boltzmann(x_fit, *popt)
                    ax_hist[idx-2].plot(x_fit, y_fit, label=fr"$p(v_{dim}) = a\cdot \exp(-b\cdot (v_{dim}- "+r"\overline{v_"+dim+"})^2)$"+f"\na={popt[0]:.2e}, b={popt[1]:.2e}")
                except RuntimeError:
                    pass
                ax_hist[idx-2].legend()
            plt.draw()

        # mark the selected rectangle
        grid = np.insert(grid, 0, Rectangle(rect[:2], rect[2], rect[3]), axis=0)
        colors = np.insert(colors, 0, (0,0,0,1), axis=0)


    if disp_particles:
        artists[0].set_offsets(np.c_[points[0], points[1]])
    artists[1].set_paths(grid)
    if disp_grid:
        artists[1].set_edgecolor("r")
    if disp_density:
        artists[1].set_facecolor(colors)




# filename_image = "data/test_8.png"
# filebasename = "data/matrix8"
# filename_image = "../data/test_11.png"
filebasename = "../data/test_19"
# filebasename = "../data/empty1x1"

display_grid = True#True
display_particles = True
display_density = False
display_vel_hist = True
bin_count = 50

save_animation = False

img_freq = 1
img_offset = 0#527
num_images = 200#//4
# num_images = 1000 - img_offset

rect = np.array([1000,1400,300,300])

alpha = .5

root = QuadTreeNode.from_csv(filebasename+".tree")
# root = QuadTreeNode(0, 0,0,1,1)

lines = []
for n in root.leafs():
    if n.contours.size != 0:
        for c in n.contours:
            lines.append([c[:2], c[2:4]])

print(lines)




# img = plt.imread(filename_image)
figsize_max = 16
if root.width > root.height:
    fig_x = figsize_max
    fig_y = figsize_max * root.height/root.width
else:
    fig_y = figsize_max
    fig_x = figsize_max * root.width/root.height

fig, ax = plt.subplots(figsize=(fig_x, fig_y))
# ax.imshow(img)
contours = ax.add_collection(LineCollection(lines, colors="orange", linewidths=3))

pts = ax.scatter([],[], color="blue", marker=".")
print(type(pts))
patch = ax.add_collection(PatchCollection([], linewidth=1, edgecolor='none', facecolor='none', alpha=alpha))
print(type(patch))

if display_vel_hist:
    ax_hist = [ax.inset_axes([.8, .7, .2, .2]), ax.inset_axes([.8, .4, .2, .2])]
else:
    ax_hist = None
#
# cells = read_cell_matrix(filebasename+".h5", 0)
# print(cells)
# particles = read_particle_matrix(filebasename+".h5", 0)
# print(particles)


ani = FuncAnimation(fig, update_plot, fargs=(img_freq, img_offset, filebasename,
                                             [pts, patch, contours],
                                             [display_grid, display_particles, display_density, display_vel_hist],
                                             [rect, ax_hist, bin_count]),
                    frames=num_images+1, interval=100, repeat=True)

#
# # ax.set_xlim(-1, 1025)
# # ax.set_ylim(-1, 1025)
# # ax.set_xlim(-1, 4001)
# # ax.set_ylim(-1, 3001)
# # ax.set_xlim(-.1,1.1)
# # ax.set_ylim(-.1,1.1)
ax.axis("equal")
ax.set_ylim(-1, root.height+1)
ax.set_xlim(-1, root.width+1)
if save_animation:
    progress_callback = lambda i, n: print(f'Saving frame {i}/{n}\n' if i%10 == 0 else "", end="")
    ani.save(f"{filebasename}.mp4", fps=24, dpi=200, progress_callback=progress_callback)
else:
    plt.show()
