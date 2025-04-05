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
    data = np.array(hf.get(f"{int(datasetID):05d}_stats"))
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

    filebasename = "../data/temperature_exchange/T1"
    filebasename = "../data/vertical_line/vl3_2"
    filebasename = "../data/cylinder/c1_b"
    filebasename = "../data/cylinder/c3"
    EXPORT_DIR = "../../../../Thesis/src/graphics/plots/"+filebasename.split("/")[-2]+"/"

    # stat_displayed = "T 1"
    # display_name = r"$T_\text{tr}$ (K)"
    species_id = 1
    display_name = r"streamlines"
    # sim_step = 1000
    sim_step = range(300, 320+1, 10)
    root = QuadTreeNode.from_csv(filebasename + ".tree")

    x,y,u,v,z = average_velocities(filebasename+".h5", sim_step, (root.width, root.height), species_id)

    root = QuadTreeNode.from_csv(filebasename + ".tree")
    # root = QuadTreeNode(0, 0,0,1,1)

    lines = []
    for n in root.leafs():
        if n.contours.size != 0:
            for c in n.contours:
                lines.append([c[:2], c[2:4]])


    img = Image.open(filebasename+".png")
    img = img.resize(np.array(z.shape)[::-1])
    img_gray = np.sum(np.asarray(img),axis=2)/3
    mask = np.flipud(img_gray <= 0)

    z = np.ma.array(z, mask=mask)
    u = np.ma.array(u, mask=mask)
    v = np.ma.array(v, mask=mask)

    fig, ax = plt.subplots()
    contours = ax.add_collection(LineCollection(lines, colors="black", linewidths=1))

    ax.streamplot(x, y, u, v, density=.5, broken_streamlines=False, color="navy")


    ax.set_aspect("equal")
    ax.set_ylim(0.5,.9)
    ax.set_xticks([])
    ax.set_yticks([])

    os.makedirs(EXPORT_DIR, exist_ok=True)
    fig.savefig(EXPORT_DIR+filebasename.split("/")[-1]+f"_streamplot.eps", dpi=300, bbox_inches="tight")

    plt.show()




