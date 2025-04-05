import os
from functools import partial

# os.environ["QT_QPA_PLATFORM"] = "wayland"
# os.environ["XDG_SESSION_TYPE"] = "xcb"

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.animation import FuncAnimation
import matplotlib.colors as mcolors
import matplotlib as mpl
import pandas as pd
import h5py
from scipy.optimize import curve_fit

from quadtree import QuadTreeNode

cmap = plt.get_cmap('inferno')
max_brightness = 1e-30

def get_stat_from_name(name):
    # format of the name: "<alias (str)> <species number (int)>"
    stat_names = {
        "numParticles": 0,
        "N": 0, # alias for numParticles
        "n": 1,
        "rho": 2,
        "cx_0": 3,
        "cy_0": 4,
        "cz_0": 5,
        "p": 6,
        "T":7
    }
    return stat_names.get(name.split(" ")[0], "error"), int(name.split(" ")[1])


def boltzmann_total(v, a, b):
    return a * v**2 * np.exp(-b*v**2)
def boltzmann(v,a,b):
    return a * np.exp(-b*v**2)


def read_cell_matrix(filename, datasetID):
    hf = h5py.File(filename, "r")
    if "00000_stats" in hf: # old format
        data = np.array(hf.get(f"{int(datasetID):05d}_tree"))
    else:
        data = np.array(hf.get(f"{int(datasetID):07d}_tree"))
    return data

def read_particle_matrix(filename, datasetID):
    hf = h5py.File(filename, "r")
    if "00000_stats" in hf: # old format
        data = hf.get(f"{int(datasetID):05d}_particles")
    else:
        data = hf.get(f"{int(datasetID):07d}_particles")
    if data is None:
        return np.array([[]])
    return np.array(data) # transpose since fortran matrices are stored columnwise

def read_stats_matrix(filename, datasetID):
    hf = h5py.File(filename, "r")
    if "00000_stats" in hf: # old format
        data = np.array(hf.get(f"{int(datasetID):05d}_stats"))
    else:
        data = np.array(hf.get(f"{int(datasetID):07d}_stats"))
    # format: cols, rows, stat_types, particle_types+overall
    return data

# def get_stream_lines(cell_matrix):
#     stream_lines = np.empty((cell_matrix.shape[0], 4), dtype=np.float64)
#     stream_lines[:, 1] = cell_matrix[:, 0] + cell_matrix[:, 2]/2
#     stream_lines[:, 2] = cell_matrix[:, 1] + cell_matrix[:, 3]/2
#     stream_lines[:, 3] = cell_matrix[:, 9]
#     stream_lines[:, 4] = cell_matrix[:, 10]
#     return stream_lines

def get_average_velocities(stats_mat, simulation_area, particle_species):
    w, h = simulation_area

    num_x = stats_mat.shape[0]
    num_y = stats_mat.shape[1]
    w /= stats_mat.shape[0]
    h /= stats_mat.shape[1]

    x,y = np.meshgrid(np.arange(num_x)*w+w/2, np.arange(num_y)*h+h/2)
    u = np.zeros_like(x)
    v = np.zeros_like(x)

    for c in range(stats_mat.shape[0]):
        for r in range(stats_mat.shape[1]):
            u[r,c] = stats_mat[c,r,get_stat_from_name(f"cx_0 {particle_species}")[0],particle_species]
            v[r,c] = stats_mat[c,r,get_stat_from_name(f"cy_0 {particle_species}")[0],particle_species]

    return x,y,u,v


def get_stream_lines(particle_matrix, averaging_rect, simulation_area):

    # def particle_in_rect(particle, rect_size, i,j):
    #     return (
    #             i*rect_size[0] <= particle[0] < (i+1)*rect_size[0]
    #             and j*rect_size[1] <= particle[1] < (j+1)*rect_size[1]
    #     )
    #
    # p_in_rect = np.vectorize(particle_in_rect, otypes=[bool])

    w,h = averaging_rect

    num_x = int(np.ceil(simulation_area[0]/w))
    num_y = int(np.ceil(simulation_area[1]/h))

    x,y = np.meshgrid(np.linspace(0, simulation_area[0], num=num_x),
                      np.linspace(0, simulation_area[1], num=num_y))

    u = np.zeros_like(x)
    v = np.zeros_like(x)

    for i in range(num_x):
        for j in range(num_y):
            if particle_matrix.size != 0:
                mask = ((i * averaging_rect[0] <= particle_matrix[0,:])
                        & (particle_matrix[0,:] <  (i + 1) * averaging_rect[0]))
                mask &= ((j * averaging_rect[1] <= particle_matrix[1,:])
                        & (particle_matrix[1,:] < (j + 1) * averaging_rect[1]))
                particles = particle_matrix[:, mask]

                u[j,i] = np.mean(particles[2,:])
                v[j,i] = np.mean(particles[3,:])

    #
    # x,y,u,v = [],[],[],[]
    #
    # for i in range(num_x):
    #     for j in range(num_y):
    #         # mask = p_in_rect(particle_matrix, averaging_rect, i,j)
    #         if particle_matrix.size != 0:
    #             mask = ((i * averaging_rect[0] <= particle_matrix[0,:])
    #                     & (particle_matrix[0,:] <  (i + 1) * averaging_rect[0]))
    #             mask &= ((j * averaging_rect[1] <= particle_matrix[1,:])
    #                     & (particle_matrix[1,:] < (j + 1) * averaging_rect[1]))
    #             particles = particle_matrix[:, mask]
    #
    #             mean_vx = np.mean(particles[2,:])
    #             mean_vy = np.mean(particles[3,:])
    #
    #         else:
    #             mean_vx = mean_vy = 0
    #
    #         x.append((i+.5)*w)
    #         y.append((j+.5)*h)
    #         u.append(mean_vx)
    #         v.append(mean_vy)

    return x,y,u,v



def get_grid_and_particles(cell_matrix, particle_matrix, cnorm):
    patches = []
    colors = []

    # global max_brightness
    # max_brightness = max_brightness[1:] + [np.max(cell_matrix[:, 5])]#
    # print(cell_matrix.shape)
    # print(np.max(cell_matrix[:, 4], axis=0), np.argmax(cell_matrix[:,4], axis=0))
    # print(cell_matrix[:, :2])
    # max_brightness = max(max_brightness, np.max(cell_matrix[:, 4], axis=0))
    # print(max_brightness)
    # cnorm = mcolors.Normalize(vmin=0, vmax=np.mean(max_brightness), clip=True)
    # cnorm = mcolors.Normalize(vmin=1e-4*max_brightness, vmax=.01*max_brightness, clip=True)
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

    if particle_matrix.size == 0:
        return [], (patches, colors)
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

def update_plot(num, fig, ax, img_freq, img_offset, filebasename, artists, display_params, hist_params, **kwargs):

    print(num)

    disp_grid, disp_particles, disp_density, disp_vel_hist, disp_streamvectors, disp_streamplot, display_stats = display_params

    rect, ax_hist, num_bins = hist_params

    # print((num*img_freq)+img_offset, end=" ")
    # data = read_cell_matrix(f"{filebasename}_{(num*img_freq)+img_offset:05d}.h5")
    cell_mat = read_cell_matrix(f"{filebasename}.h5", (num*img_freq)+img_offset)
    print("num leafs:", cell_mat.shape[0])
    # print(cell_mat)
    particle_mat = read_particle_matrix(f"{filebasename}.h5", (num*img_freq)+img_offset)
    # print(particle_mat)
    print("num particles:", particle_mat.shape[1])
    points, (grid, colors) = get_grid_and_particles(cell_mat, particle_mat, cnorm=kwargs["cnorm"])

    stats_mat = read_stats_matrix(f"{filebasename}.h5", (num*img_freq)+img_offset)
    print(stats_mat.shape)

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
        if points:
            artists[0].set_offsets(np.c_[points[0], points[1]])
    artists[1].set_paths(grid)
    if disp_grid:
        artists[1].set_edgecolor("r")
    if disp_density:
        artists[1].set_facecolor(colors)

    if disp_streamvectors or disp_streamplot:
        # x,y,u,v = get_stream_lines(particle_mat, kwargs["averaging_area"], kwargs["simulation_area"])
        x,y,u,v = get_average_velocities(stats_mat, kwargs["simulation_area"], get_stat_from_name(kwargs["stat_displayed"])[1])
        if disp_streamvectors:
            artists[3].set_UVC(u,v)
        if disp_streamplot:
            if artists[4] is not None:
                artists[4].lines.remove()
                # artists[4].arrows.remove()
                artists[4].arrows.set_paths([])
                for art in ax.get_children():
                    if not isinstance(art, FancyArrowPatch):
                        continue
                    art.remove()
                    # for a in artists[4].arrows:
                #     print(a)
                # try:
                #     artists[4].arrows.remove()
                # except:
                #     ...
            new_streamplot = ax.streamplot(x,y,u,v, broken_streamlines=False, color="black")
            artists[4] = new_streamplot

    if display_stats:
        stat_id, species_id = get_stat_from_name(kwargs["stat_displayed"])
        w,h = kwargs["simulation_area"]
        w /= stats_mat.shape[0]
        h /= stats_mat.shape[1]
        statistics_cells = []
        statistics_data = []
        # cnorm = mcolors.Normalize(vmin=0, vmax=np.max(stats_mat[:,:,stat_id,species_id]))
        cnorm = kwargs.get("cnorm")
        print("max", np.max(stats_mat[:, :, stat_id, species_id]))
        print("min", np.min(stats_mat[:, :, stat_id, species_id]))
        for c in range(stats_mat.shape[0]):
            for r in range(stats_mat.shape[1]):
                statistics_cells.append(Rectangle((c*w, r*h), w,h))
                statistics_data.append(cmap(cnorm(stats_mat[c,r,stat_id,species_id])))

        artists[5].set_paths(statistics_cells)
        artists[5].set_facecolor(statistics_data)
        # artists[5].set_edgecolor("g")





# filename_image = "data/test_8.png"
# filebasename = "data/matrix8"
# filename_image = "../data/test_11.png"
# filebasename = "../data/test_19"
filebasename = "../data/soundwaves/soundwaves2"
filebasename = "../data/vertical_line/vl2"
filebasename = "../data/vertical_line/vl3_2"
# filebasename = "../data/vertical_line/triangle1"
filebasename = "../data/temperature_exchange/T2"
# filebasename = "../data/cylinder/c1_b"
# filebasename = "../data/cylinder/c4"
filebasename = "../data/shear_flow/s2"
# filebasename = "../data/apollo_cm/cm8"
# filebasename = "../data/debug/d1"
# filebasename = "../data/wing1"
# filebasename = "../data/empty1x1"

display_grid = False
display_particles = True
display_density = False
display_vel_hist = False
bin_count = 50

display_stream = False
display_streamplot = False
# divs_x = 10
# divs_y = 10
arrow_scale = 50000

display_stats = True
stat_displayed = "T 0"

cnorm = mcolors.Normalize(vmin=0, vmax=6000, clip=True)
# cnorm = mcolors.Normalize(vmin=5, vmax=10, clip=False)
# cnorm = mcolors.Normalize(vmin=0, vmax=100, clip=True)
cnorm = mcolors.Normalize(vmin=40, vmax=120, clip=True)
# cnorm = mcolors.LogNorm(vmin=1e20, vmax=1e22, clip=True)
# cnorm = mcolors.Normalize(vmin=1e13, vmax=3e13, clip=True)
# cnorm = mcolors.LogNorm(vmin=1e18, vmax=1e20, clip=True)

save_animation = False

img_freq = 500
img_offset = 0#527
num_images = 100000//img_freq#//4
# num_images = 1000 - img_offset
img_freq = 100
img_offset = 0#527
num_images = 10000//img_freq#//4
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
    fig_y = (figsize_max - 2) * root.height/root.width
else:
    fig_y = figsize_max - 2
    fig_x = (figsize_max - 2) * root.width/root.height + 2

fig, ax = plt.subplots(figsize=(fig_x, fig_y))
# ax.imshow(img)
contours = ax.add_collection(LineCollection(lines, colors="orange", linewidths=3))

pts = ax.scatter([],[], color="blue", marker=".")
print(type(pts))
patch = ax.add_collection(PatchCollection([], linewidth=1, edgecolor='none', facecolor='none', alpha=alpha))
print(type(patch))
patch2 = ax.add_collection(PatchCollection([], linewidth=1, edgecolor='none', facecolor='none', alpha=alpha))


colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=cnorm), ax=ax, orientation="vertical", label=stat_displayed)

if display_vel_hist:
    ax_hist = [ax.inset_axes([.8, .7, .2, .2]), ax.inset_axes([.8, .4, .2, .2])]
else:
    ax_hist = None
#
# cells = read_cell_matrix(filebasename+".h5", 0)
# print(cells)
# particles = read_particle_matrix(filebasename+".h5", 0)
# print(particles)

if display_stream:
    # x,y,u,v = get_stream_lines(np.array([[]]),
    #                            (root.width/divs_x, root.height/divs_y),
    #                            (root.width, root.height))
    x,y,u,v = get_average_velocities(read_stats_matrix(filebasename+".h5", 0),
                                     simulation_area=(root.width, root.height),
                                     particle_species=0)
    stream = ax.quiver(x,y,u,v, scale=arrow_scale)
    print(type(stream))
else:
    stream = None


# ani = FuncAnimation(fig, update_plot, fargs=(img_freq, img_offset, filebasename,
#                                              [pts, patch, contours, stream],
#                                              [display_grid, display_particles, display_density, display_vel_hist, display_stream],
#                                              [rect, ax_hist, bin_count]),
#                     frames=num_images+1, interval=100, repeat=True)

ani = FuncAnimation(fig, partial(update_plot,
                                 fig=fig, ax=ax,
                                 img_freq=img_freq, img_offset=img_offset, filebasename=filebasename,
                                 artists=[pts, patch, contours, stream, None, patch2],
                                 display_params=[display_grid, display_particles, display_density, display_vel_hist, display_stream, display_streamplot, display_stats],
                                 hist_params=[rect, ax_hist, bin_count],
                                 simulation_area=(root.width, root.height),
                                 # averaging_area=(root.width/divs_x, root.height/divs_y),
                                 stat_displayed=stat_displayed,
                                 cnorm=cnorm
                                 ),
                    frames=num_images+1, interval=100, repeat=True)

#
# # ax.set_xlim(-1, 1025)
# # ax.set_ylim(-1, 1025)
# # ax.set_xlim(-1, 4001)
# # ax.set_ylim(-1, 3001)
# # ax.set_xlim(-.1,1.1)
# # ax.set_ylim(-.1,1.1)
# ax.axis("equal")
ax.set_ylim(-1, root.height+1)
ax.set_xlim(-1, root.width+1)
if save_animation:
    progress_callback = lambda i, n: print(f'Saving frame {i}/{n}\n' if i%10 == 0 else "", end="")
    ani.save(f"{filebasename}.mp4", fps=5, dpi=200, progress_callback=progress_callback)
else:
    plt.show()
