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

def read_stats_matrix(filename, datasetID):
    hf = h5py.File(filename, "r")
    if "00000_stats" in hf: # old format
        data = np.array(hf.get(f"{int(datasetID):05d}_stats"))
    else:
        data = np.array(hf.get(f"{int(datasetID):07d}_stats"))
    # format: cols, rows, stat_types, particle_types+overall
    return data

def get_average_stat(stats_mat, simulation_area, stat):
    w, h = simulation_area

    num_x = stats_mat.shape[0]
    num_y = stats_mat.shape[1]
    w /= stats_mat.shape[0]
    h /= stats_mat.shape[1]

    x,y = np.meshgrid(np.arange(num_x)*w+w/2, np.arange(num_y)*h+h/2)
    z = np.zeros_like(x)

    stat_id, species_id = get_stat_from_name(stat)

    for c in range(stats_mat.shape[0]):
        for r in range(stats_mat.shape[1]):
            z[r,c] = stats_mat[c,r,stat_id,species_id]

    return x,y,z

def average_stat(filename, step_range, simulation_area, stat):
    w, h = simulation_area

    step_list = list(step_range)

    stats_mat = read_stats_matrix(filename, step_list[0])

    num_x = stats_mat.shape[0]
    num_y = stats_mat.shape[1]
    w /= stats_mat.shape[0]
    h /= stats_mat.shape[1]

    x, y = np.meshgrid(np.arange(num_x) * w + w / 2, np.arange(num_y) * h + h / 2)
    z = np.zeros_like(x)

    stat_id, species_id = get_stat_from_name(stat)


    for step in step_list:
        stats_mat = read_stats_matrix(filename, step)
        for c in range(num_x):
            for r in range(num_y):
                z[r, c] += stats_mat[c, r, stat_id, species_id]

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
                u[r,c] += stats_mat[c,r,get_stat_from_name(f"cx_0 {species_id}")[0],species_id]
                v[r,c] += stats_mat[c,r,get_stat_from_name(f"cy_0 {species_id}")[0],species_id]
                z[r,c] += np.sqrt(u[r,c]**2+v[r,c]**2)

    u /= len(step_list)
    v /= len(step_list)
    z /= len(step_list)

    return x, y, u,v,z



if __name__ == "__main__":

    filebasename = "../data/temperature_exchange/T2"
    filebasename = "../data/shear_flow/s2"
    # filebasename = "../data/vertical_line/vl3_2"
    # filebasename = "../data/cylinder/c1_b"
    # filebasename = "../data/cylinder/c3"
    # filebasename = "../data/apollo_cm/cm7"
    EXPORT_DIR = "../../../../Thesis/src/graphics/plots/"+filebasename.split("/")[-2]+"/"

    stat_displayed = "T 0"
    display_name = r"$T_\text{tr}$ (K)"
    # stat_displayed = "p 0"
    # display_name = r"$p$ (Pa)"
    stat_displayed = "n 0"
    display_name = r"$n$ (m$^{-3}$)"
    levels = None
    # levels = np.linspace(0, 7000, 15)
    levels = np.logspace(20, 22.5, 11)
    # sim_step = 1000

    show_streamlines = False
    sim_step = range(98000, 99800+1, 200)
    sim_step = range(49000, 50000+1, 500)
    sim_step = range(49000, 50000+1, 500)
    sim_step = range(9500, 10000+1, 500)
    # sim_step = range(98000, 99500+1, 500)
    # sim_step = range(20000, 20000+1, 500)
    # sim_step = range(30000, 31000+1, 500)
    # sim_step = range(19000, 20000+1, 100)
    # sim_step = range(80000, 99000+1, 500)
    # sim_step = range(13000, 14600+1, 200)

    if levels is not None:
        cnorm = mcolors.LogNorm(vmin=1e20, vmax=5e22, clip=True)
        # cnorm = mcolors.Normalize(vmin=0, vmax=7000, clip=False)
        # cnorm = mcolors.Normalize(vmin=0, vmax=11000, clip=False)
    else:
        # cnorm = mcolors.Normalize(vmin=0, vmax=6000, clip=True)
        cnorm = mcolors.Normalize(vmin=0, vmax=11000, clip=True)
        cnorm = mcolors.Normalize(vmin=0, vmax=6500, clip=True)
        # cnorm = mcolors.Normalize(vmin=0, vmax=0.02, clip=True)
        cnorm = mcolors.Normalize(vmin=0, vmax=10, clip=True)
        # cnorm = mcolors.Normalize(vmin=40, vmax=120, clip=True)
        cnorm = mcolors.LogNorm(vmin=1e20, vmax=5e22, clip=True)

    root = QuadTreeNode.from_csv(filebasename + ".tree")

    if isinstance(sim_step, (Iterable, Generator, Sequence)):
    # if type(sim_step) in [type(list), type(range)] or isinstance(sim_step, Iterable):
        x,y,z = average_stat(filebasename+".h5", sim_step, (root.width, root.height), stat_displayed)
    else:
        stats = read_stats_matrix(filebasename+".h5", sim_step)
        x,y,z = get_average_stat(stats, (root.width, root.height), stat_displayed)

    if show_streamlines:
        x_,y_,u,v,z_ = average_velocities(filebasename+".h5", sim_step, (root.width, root.height), get_stat_from_name(stat_displayed)[1])

    root = QuadTreeNode.from_csv(filebasename + ".tree")
    # root = QuadTreeNode(0, 0,0,1,1)

    lines = []
    for n in root.leafs():
        if n.contours.size != 0:
            for c in n.contours:
                lines.append([c[:2], c[2:4]])
    try:
        img = Image.open(filebasename+".png")
        img = img.resize(np.array(z.shape)[::-1])
        img_gray = np.sum(np.asarray(img),axis=2)/3
        mask = np.flipud(img_gray <= 100)
        print(np.max(mask))

        # z = np.ma.array(z, mask=mask)
    except FileNotFoundError:
        print("no image found")

    fig, ax = plt.subplots(tight_layout=True)
    contours = ax.add_collection(LineCollection(lines, colors="black", linewidths=1))

    if levels is not None:
        CS = ax.contourf(x,y,z, levels=levels, norm=cnorm, cmap="rainbow")
    else:
        CS = ax.contourf(x,y,z, norm=cnorm, cmap="rainbow")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)  # Adjust "pad" to move it closer
    # colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=cnorm), ax=ax, orientation="horizontal", label=stat_displayed)
    colorbar = fig.colorbar(CS, ax=ax, orientation="horizontal", label=display_name, cax=cax)
    colorbar.set_label(label=display_name, size=12)
    colorbar.ax.tick_params(labelsize=12)

    if show_streamlines:
        ax.streamplot(x_, y_, u, v, density=.5, broken_streamlines=False, color="white")

    ax.set_aspect("equal")
    ax.set_ylim(.5,.85)
    # ax.set_ylim(1,1.9)
    # ax.set_ylim(0,2)
    # ax.set_ylim(0,10)
    # ax.set_ylim(3.5,7.5)
    # ax.set_xlim(1.5,5.5)
    ax.set_ylim(0,0.01)
    ax.set_xticks([])
    ax.set_yticks([])

    os.makedirs(EXPORT_DIR, exist_ok=True)
    fig.savefig(EXPORT_DIR+filebasename.split("/")[-1]+f"_{stat_displayed[0]}.eps", dpi=300, bbox_inches="tight")

    plt.show()




