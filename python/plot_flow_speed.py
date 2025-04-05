import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import PatchCollection, LineCollection
import h5py
from PIL import Image
import os

from collections.abc import Iterable, Generator, Sequence

from matplotlib.pyplot import tight_layout

from quadtree import QuadTreeNode


def read_stats_matrix(filename, datasetID):
    hf = h5py.File(filename, "r")
    if "00000_stats" in hf: # old format
        data = np.array(hf.get(f"{int(datasetID):05d}_stats"))
    else:
        data = np.array(hf.get(f"{int(datasetID):07d}_stats"))
    # format: cols, rows, stat_types, particle_types+overall
    return data



def average_flow_speed(filename, step_range, simulation_area, species_id):
    w, h = simulation_area

    step_list = list(step_range)

    stats_mat = read_stats_matrix(filename, step_list[0])

    num_x = stats_mat.shape[0]
    num_y = stats_mat.shape[1]
    w /= stats_mat.shape[0]
    h /= stats_mat.shape[1]

    x, y = np.meshgrid(np.arange(num_x) * w + w / 2, np.arange(num_y) * h + h / 2)
    z = np.zeros_like(x)


    for step in step_list:
        stats_mat = read_stats_matrix(filename, step)
        for c in range(num_x):
            for r in range(num_y):
                z[r, c] += np.linalg.norm(stats_mat[c, r, 3:6, species_id])

    z /= len(step_list)

    return x, y, z

def average_velocities(filename, step_range, simulation_area, species_id):
    w, h = simulation_area

    step_list = list(step_range)

    stats_mat = read_stats_matrix(filename, step_list[0])

    num_x = stats_mat.shape[0]
    num_y = stats_mat.shape[1]
    w /= stats_mat.shape[0]
    h /= stats_mat.shape[1]

    x, y = np.meshgrid(np.arange(num_x) * w + w / 2, np.arange(num_y) * h + h / 2)
    u = np.zeros_like(x)
    v = np.zeros_like(x)
    z = np.zeros_like(x)

    for step in step_list:
        stats_mat = read_stats_matrix(filename, step)
        for c in range(stats_mat.shape[0]):
            for r in range(stats_mat.shape[1]):
                u[r,c] += stats_mat[c,r,3,species_id]
                v[r,c] += stats_mat[c,r,4,species_id]
                z[r,c] += np.sqrt(u[r,c]**2+v[r,c]**2)

    u /= len(step_list)
    v /= len(step_list)
    z /= len(step_list)

    return x, y, u,v,z



if __name__ == "__main__":

    filebasename = "../data/temperature_exchange/T1"
    filebasename = "../data/vertical_line/vl3_2"
    filebasename = "../data/cylinder/c1_b"
    filebasename = "../data/cylinder/c3"
    # filebasename = "../data/apollo_cm/cm7"
    filebasename = "../data/shear_flow/s2"
    EXPORT_DIR = "../../../../Thesis/src/graphics/plots/"+filebasename.split("/")[-2]+"/"

    species_id = 0
    temperature = 200 # in Kelvin
    molar_mass = 39.948 * 1e-3 # of Argon
    # sim_step = 1000
    sim_step = range(9500, 10000 + 1, 100)
    # sim_step = range(99000, 99900+1, 100)
    # sim_step = range(99000, 99000+1, 200)
    # sim_step = range(49000, 50000+1, 500)
    # sim_step = range(1900, 2000+1, 100)
    levels = None
    # levels = np.logspace(3, 4, 9)

    show_streamlines = True

    # cnorm = mcolors.LogNorm(vmin=1e20, vmax=1e21, clip=True)
    cnorm = mcolors.Normalize(vmin=0, vmax=3500, clip=True)
    cnorm = mcolors.Normalize(vmin=0, vmax=80, clip=True)
    # cnorm = mcolors.LogNorm(vmin=1, vmax=1e4, clip=True)

    root = QuadTreeNode.from_csv(filebasename + ".tree")
    # root = QuadTreeNode(0, 0,0,1,1)

    lines = []
    for n in root.leafs():
        if n.contours.size != 0:
            for c in n.contours:
                lines.append([c[:2], c[2:4]])

    x, y, z = average_flow_speed(filebasename + ".h5", sim_step, (root.width, root.height), species_id)
    print(z.shape)

    if show_streamlines:
        x_,y_,u,v,z_ = average_velocities(filebasename+".h5", sim_step, (root.width, root.height), species_id)

    try:
        img = Image.open(filebasename+".png")
        img = img.resize(np.array(z.shape)[::-1])
        img_gray = np.sum(np.asarray(img),axis=2)/3
        mask = np.flipud(img_gray <= 150)

        z = np.ma.array(z, mask=mask)
    except FileNotFoundError:
        print("no image found")

    fig, ax = plt.subplots(tight_layout=True)
    contours = ax.add_collection(LineCollection(lines, colors="black", linewidths=1))

    # ax.imshow(img_gray, extent=(root.x, root.x+root.width, root.y, root.y+root.height), cmap="gray")
    if levels is not None:
        CS = ax.contourf(x,y,z, levels=levels, norm=cnorm, cmap="rainbow")
    else:
        CS = ax.contourf(x,y,z, norm=cnorm, cmap="rainbow")
    # colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=cnorm), ax=ax, orientation="horizontal", label=f"Mach number")
    # Use make_axes_locatable to position the colorbar closer

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)  # Adjust "pad" to move it closer

    colorbar = fig.colorbar(CS, ax=ax, orientation="horizontal", label=f"flow speed (m/s)", cax=cax)
    colorbar.set_label(label="flow speed (m/s)", size=12)
    colorbar.ax.tick_params(labelsize=12)

    if show_streamlines:
        # ax.streamplot(x_, y_, u, v, density=.35, broken_streamlines=False, color="white", linewidth=.7, arrowsize=.7, arrowstyle='-|>')
        ax.streamplot(x_, y_, u, v, density=.2, broken_streamlines=False, color="white", linewidth=.7, arrowsize=.7, arrowstyle='-|>')

    # ax.set_xlim(0,0.85)
    ax.set_aspect("equal")
    ax.set_ylim(0.5,0.85)
    # ax.set_ylim(0.5,.9)
    # ax.set_ylim(0,10)
    ax.set_ylim(0,.01)
    ax.set_xticks([])
    ax.set_yticks([])

    os.makedirs(EXPORT_DIR, exist_ok=True)
    fig.savefig(EXPORT_DIR+filebasename.split("/")[-1]+f"_flowSpeed_species{species_id}.eps", dpi=300, bbox_inches="tight")

    plt.show()




