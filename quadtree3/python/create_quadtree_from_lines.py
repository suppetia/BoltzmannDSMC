import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection, LineCollection

from quadtree import (
    QuadTreeNode,
    createSubTree,
    line_intersection,
    snap_point_to_grid,
    getRectangles
)

def getPolyGridApproximation(contour_lines, tree_root:QuadTreeNode, first_iteration=True):
    # in the first iteration find all actual crossings
    # in the second iteration find cells where no structure was found yet and check for crossings of an edge
    # if a structure just touches the edge of a cell, introduce a very small shift for the line

    a = None
    b = None
    for node in tree_root.leafs():
        if not first_iteration:
            if node.contours:
                continue
        for line in contour_lines:
            left_border = [[node.x, node.y], [node.x, node.y+node.height]]
            lower_border = [[node.x, node.y], [node.x+node.width, node.y]]
            right_border = [[node.x+node.width, node.y], [node.x+node.width, node.y+node.height]]
            upper_border = [[node.x, node.y+node.height], [node.x+node.width, node.y+node.height]]

            # normal vector
            n = np.array([(line[1,1] - line[0,1]), -(line[1,0] - line[0,0])])
            n /= np.sqrt(n[0] ** 2 + n[1] ** 2)
            # print(n)


            for i, b1 in enumerate([lower_border, right_border, upper_border, left_border]):
                a = line_intersection(line, b1)
                if a is None:
                    continue
                b = None
                for b2 in [lower_border, right_border, upper_border, left_border][i+1:]:
                    b = line_intersection(line, b2)
                    if b is None:
                        continue
                    if np.isclose(a,b).all(): # if the intersection happens at a corner
                        if first_iteration:
                            b = None
                            continue

                        # slope of the line
                        m = (line[1, 1] - line[0, 1])/ (line[1,0] - line[0,0])
                        if np.allclose(a, [node.x, node.y]): # lower left
                            a[0] += 1e-10
                            b[1] -= m * 1e-10 # m is in this case always negative
                        elif np.allclose(a, [node.x+node.width, node.y]): # lower right
                            a[0] -= 1e-10
                            b[1] += m * 1e-10 # m is in this case always positive
                        elif np.allclose(a, [node.x, node.y+node.height]): #upper left
                            a[0] += 1e-10
                            b[1] -= m * 1e-10 # m is in this case always positive
                        elif np.allclose(a, [node.x+node.width, node.y+node.height]): # upper right
                            a[0] -= 1e-10
                            b[1] += m * 1e-10 # m is in this case always negative
                    # print(a, b)
                    break
                break
            else:
                continue
            if (a is None) or (b is None):
                # if no intersection was found
                continue

            # check the orientation of the line -> standard orientation is when a <= b
            standard_orientation = all(
                [(line[1,0]-line[0,0] < 0) == (b[0] - a[0] < 0),
                 (line[1,1]-line[0,1] < 0) == (b[1] - a[1] < 0)]
            )
            # print(standard_orientation)
            # print(a,b)
            # print(b[0]-a[0], np.sqrt((b[0]-a[0])**2+(b[1]-a[1])**2))
            if not standard_orientation:
                a,b = b,a
            if b[1] - a[1] >= 0:
                gamma = np.arccos((b[0]-a[0])/np.sqrt((b[0]-a[0])**2+(b[1]-a[1])**2)) #+ np.pi
            else:
                gamma = np.arccos(-(b[0]-a[0])/np.sqrt((b[0]-a[0])**2+(b[1]-a[1])**2)) + np.pi

            node.contour_angles.append(gamma)

            # # normal vector
            # n = np.array([(b[1]-a[1]), -(b[0]-a[0])])
            # n /= np.sqrt(n[0]**2+n[1]**2)
            # print(n)
            node.contours.append([a[0], a[1], b[0], b[1], n[0], n[1]])

def snap_lines_to_grid_without_redundancy_check(lines, x,y,dx,dy):
    # snap lines to grid with grid spacing dx,dy and an offset x,y
    new_lines = []
    for line in lines:
        print([snap_point_to_grid(line[0], x, y, dx, dy), snap_point_to_grid(line[1], x, y, dx, dy)])
        new_lines.append([snap_point_to_grid(line[0],x,y,dx,dy), snap_point_to_grid(line[1],x,y,dx,dy)])

    return np.array(new_lines)


if __name__ == "__main__":

    filename = "../data/vertical_line/vl2"

    width=30
    height=15

    quadtree_depth = 3

    # vertical line
    x_offset = 1e-9
    lines = np.array([
        [[15 + x_offset, 0], [15+x_offset, 5]],
        [[15 + x_offset, 5], [15, 5 + x_offset]],
        [[15, 5 + x_offset], [15-x_offset, 5]],
        [[15 - x_offset, 5], [15-x_offset, 0]],
        [[width, 0], [0,0]]
    ])

    show_degrees=False#True
    show_normals=True
    show_cells=True

    ### end of configuation area


    root = QuadTreeNode(0, 0,0, width, height)

    createSubTree(root, quadtree_depth)

    snapped_lines = snap_lines_to_grid_without_redundancy_check(lines, 0,0, width*2**(-quadtree_depth), height*2**(-quadtree_depth))

    getPolyGridApproximation(snapped_lines, root)
    getPolyGridApproximation(snapped_lines, root, first_iteration=False)

    root.mergeSubtreeCells()

    # draw the quadtree (representation as the simulation area)
    fig, ax = plt.subplots(figsize=(20,20))
    rects = getRectangles(root)

    col_rect = ax.add_collection(PatchCollection(rects, linewidth=.1, edgecolor="r", facecolor='none'))


    lines = []
    normals = []
    for n in root.leafs():
        # print(n.contours)
        if np.array(n.contours).size > 0:
            if show_cells:
                ax.fill_between([n.x, n.x+n.width], [n.y]*2, [n.y+n.height]*2, color="orange", alpha=0.3)
            for i,c in enumerate(n.contours):
                # print(c)
                # print([c[:2], c[2:4]])
                c = np.array(c)
                lines.append([c[:2], c[2:4]])
                if show_normals:
                    normals.append([c[:2], c[:2]+c[4:]])
                if show_degrees:
                    # ax.text((c[0]+c[2])/2, (c[1]+c[3])/2, f"{np.arcsin(c[4])*180/np.pi:.1f}", color="pink", fontsize=20)
                    ax.text((c[0]+c[2])/2, (c[1]+c[3])/2, f"{n.contour_angles[i]*180/np.pi:.1f}", color="pink", fontsize=20)

    approx_lines = ax.add_collection(LineCollection(lines, colors="red", linewidths=.75))
    if normals:
        print(normals)
        normal_lines = ax.add_collection(LineCollection(normals[:], colors="green", linewidths=.5))

    ax.set_xlim(root.x - .1, root.width + .1)
    ax.set_ylim(root.y - .1, root.height + .1)

    ax.axis("equal")

    plt.show()
    #
    # # save the contours
    # np.save(filename+".npy", snapped_lines)

    root.to_csv(filename+".tree")
    # print(len(root.leafs()))

    # new_root = QuadTreeNode.from_csv(".".join(filename.split(".")[:-1])+".tree")
