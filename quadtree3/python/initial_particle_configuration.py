import os.path

from scipy.constants import Boltzmann as k_B

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from quadtree import QuadTreeNode
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.path import Path
import matplotlib.patches as patches

rng = np.random.default_rng()

def spawn_particles(density, area, velocity_distribution, particle_distribution):


    if particle_distribution.lower() in ["u", "uniform"]:
        number = density*area[2]*area[3]
        particles = np.empty((number, 5), dtype=float)
        particles[:, 0] = area[0] + rng.random((number,)) * area[2]
        particles[:, 1] = area[1] + rng.random((number,)) * area[3]
    elif particle_distribution.lower() in ["e", "equidistant"]:
        dist = np.sqrt(1/density)
        num_x = int(area[2]/dist)
        num_y = int(area[3]/dist)
        number = num_x*num_y
        particles = np.empty((number, 5), dtype=float)
        print(num_x, num_y)
        print(density*area[2]*area[3])
        print(density*area[2]*area[3] - number)
        for y in range(num_y):
            particles[y*num_x:(y+1)*num_x, 0] = area[0] + np.linspace(0, area[2], num_x)
            particles[y*num_x:(y+1)*num_x, 1] = area[1] + y * area[3]/num_y

    else:
        raise ValueError(f"particle_distribution '{particle_distribution}' not implemented. Use 'uniform' or 'equidistant'.")

    particles[:, 2] = rng.normal(velocity_distribution[0], velocity_distribution[3], size=number)
    particles[:, 3] = rng.normal(velocity_distribution[1], velocity_distribution[4], size=number)
    particles[:, 4] = rng.normal(velocity_distribution[2], velocity_distribution[5], size=number)

    return particles

# def polygon_from_curves(lines):
#     lines = lines[:]
#     polygons = []
#     while lines: # is not empty
#         line = lines.pop(0)
#         poly = [line[0]]
#         last_point = line[1]
#
#         while lines:
#             for i in range(len(lines)):
#                 if (lines[i][0] == last_point).all():
#                     line = lines.pop(i)
#                     poly.append(line[0])
#                     last_point = line[1]
#                     break
#             if (last_point == poly[0]).all():
#                 break
#         polygons.append(poly)
#
#     return polygons


def check_particles_in_polygon(particles, polygon):

    is_in_poly = np.empty(particles.shape[0], dtype=bool)
    #https://stackoverflow.com/questions/36399381/wqhats-the-fastest-way-of-checking-if-a-point-is-inside-a-polygon-in-python
    #https://matplotlib.org/stable/api/path_api.html

    poly = Path([l[0] for l in polygon]+[polygon[-1,1]], closed=True)
    for i, p in enumerate(particles):
        is_in_poly[i] = poly.contains_point(p[:2])

    return is_in_poly, poly


if __name__ == "__main__":
    # filebasename = "../data/test_19"
    # filebasename = "../data/soundwaves/soundwaves2"
    # filebasename = "../data/apollo_cm/cm2"
    filebasename = "../data/vertical_line/vl2"


    root = QuadTreeNode.from_csv(filebasename + ".tree")
    if os.path.exists(filebasename+".npy"):
        polygon = np.load(filebasename + ".npy")
    else:
        polygon = None

    particle_distribution = "e" # uniform or equidistant
    T = 300 # temperature
    particle_mass = 4.6518e-27 # mass of N2
    # num_particles = 1000
    particle_density = 1e1
    # spawn_area = [100,100,root.width*.5, root.height*.8]
    spawn_area = [0,0,root.width, root.height]

    sigma = np.sqrt(k_B * T/particle_mass)
    print(sigma)
    velocity_distribution = [1e3,0,0, # v_x, v_y, v_z
                             sigma,sigma,sigma # s_x, s_y, s_z
                             ]


    lines = []
    for n in root.leafs():
        if n.contours.size != 0:
            for c in n.contours:
                lines.append([c[:2], c[2:4]])

    figsize_max = 16
    if root.width > root.height:
        fig_x = figsize_max
        fig_y = figsize_max * root.height / root.width
    else:
        fig_y = figsize_max
        fig_x = figsize_max * root.width / root.height

    fig, ax = plt.subplots(figsize=(fig_x, fig_y))
    contours = ax.add_collection(LineCollection(lines, colors="orange", linewidths=3))

    particles = spawn_particles(particle_density, spawn_area, velocity_distribution, particle_distribution)

    if len(lines) > 0:
        if polygon is not None:
            mask_in_poly, poly = check_particles_in_polygon(particles, polygon)

            ax.add_patch(patches.PathPatch(poly, facecolor="red"))

            # ax.scatter(particles[~mask_in_poly, 0], particles[~mask_in_poly, 1], color="red", marker=".")

            particles = particles[~mask_in_poly]
    ax.scatter(particles[:, 0], particles[:, 1], color="blue", marker=".")

    nrows, ncols = particles.shape

    output_filename = filebasename + ".particles"
    df = pd.DataFrame(index=[1]*nrows, data=particles) # store the particle types in the index, for now just use 1
    with open(output_filename, "w") as f:
        f.write(f"{nrows}\n")
    df.to_csv(output_filename, sep=" ", header=False, mode="a")
    # np.savetxt(output_filename, particles, header=f"{nrows}")


    ax.set_ylim(-1, root.height+1)
    ax.set_xlim(-1, root.width+1)
    ax.axis("equal")

    plt.show()