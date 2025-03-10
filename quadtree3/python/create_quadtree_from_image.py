import matplotlib
import numpy as np
import cv2
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection, LineCollection

from quadtree import (
    QuadTreeNode,
    snap_lines_to_grid,
    snap_point_to_grid,
    line_intersection,
    getRectangles
)


def createSubTreeFromImage(root:QuadTreeNode, depth:int, img:np.ndarray):

    if depth == root.level:
        return
    offset = 15
    if np.all(img[max(0,int(root.y)-offset):min(img.shape[0],int(root.y)+int(np.ceil(root.height))+offset), max(0, int(root.x)-offset):min(img.shape[1],int(root.x)+int(np.ceil(root.width))+offset)] == 0):
        return
    elif np.any(img[max(0,int(root.y)-offset):min(img.shape[0],int(root.y)+int(np.ceil(root.height))+offset), max(0, int(root.x)-offset):min(img.shape[1],int(root.x)+int(np.ceil(root.width))+offset)] == 0):
        root.split()
        for i in range(len(root.children)):
            createSubTreeFromImage(root.children[i], depth, img)
    return


def getPolyApproximation(contour):
    epsilon = 0.005 * cv2.arcLength(contour, True)  # Tolerance factor
    approx_contour = cv2.approxPolyDP(contour, epsilon, True)  # Approximate contour
    # pprint(approx_contour)

    lines = np.array([[approx_contour[i,0,:], approx_contour[i+1,0,:]] for i in range(len(approx_contour)-1)] + [[approx_contour[-1,0,:], approx_contour[0,0,:]]])

    return lines


def simplify_lines(lines, x, y, dx, dy):
    lines_ = lines.copy()
    new_lines = []
    for i,l in enumerate(lines_):
        l1 = None
        l2 = None
        if l[0,1] == l[1,1]:
            l1 = snap_point_to_grid(l[0,:], x,y,dx,dy, snap_y=False)
            l2 = snap_point_to_grid(l[1,:], x,y,dx,dy, snap_y=False)
        if l[0,0] == l[1,0]:
            l1 = snap_point_to_grid(l[0,:], x,y,dx,dy, snap_x=False)
            l2 = snap_point_to_grid(l[1,:], x,y,dx,dy, snap_x=False)

        if l1 is None and l2 is None:
            new_lines.append(l)
        else:
            if len(new_lines) > 0:
                new_lines[-1][1,:] = l1
            else:
                lines_[-1][1,:] = l1
            if not all(l1 == l2):
                new_lines.append(np.array([l1, l2]))
            if not i == len(lines_)-1: # if not the last point update the next start
                lines_[i+1][0,:] = l2
            else: # in case of the last line update the beginning of the first line
                new_lines[0][0,:] = l2
    return new_lines



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
        # plt.figure()
        # plt.plot(line[:, 0], line[:, 1])
        # plt.scatter(line[0, 0], line[0, 1], label="line_a")
        # plt.scatter(line[1, 0], line[1, 1], label="line_b")
        # # plt.scatter(line[:,0], line[1,:])
        # # plt.scatter(np.array([a,b])
        #
        # rects = getRectangles(tree_root)
        #
        # plt.gca().add_collection(PatchCollection(rects, linewidth=.1, edgecolor="r", facecolor='none'))
        # plt.legend()
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

            # print(np.isclose(a,b).all())
            # print(a, b)
            #
            # plt.figure()
            # plt.plot(line[:, 0], line[:, 1])
            # plt.scatter(line[0, 0], line[0, 1], label="line_a")
            # plt.scatter(line[1, 0], line[1, 1], label="line_b")
            # plt.scatter(a[0], a[1], label="a")
            # plt.scatter(b[0], b[1], label="b")
            # # plt.scatter(line[:,0], line[1,:])
            # # plt.scatter(np.array([a,b])
            #
            # rects = getRectangles(tree_root)
            #
            # plt.gca().add_collection(PatchCollection(rects, linewidth=.1, edgecolor="r", facecolor='none'))
            # plt.legend()
            # plt.show()

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
            # node.contours.append([a[0], a[1], b[0], b[1], np.sin(gamma), np.cos(gamma)])

            # # normal vector
            # n = np.array([(b[1]-a[1]), -(b[0]-a[0])])
            # n /= np.sqrt(n[0]**2+n[1]**2)
            # print(n)
            node.contours.append([a[0], a[1], b[0], b[1], n[0], n[1]])

            # rotMat = np.array([[np.cos(np.pi/2), np.sin(np.pi/2)], [-np.sin(np.pi/2), np.cos(np.pi/2)]])
            #
            # n2 = np.array(rotMat @ (b-a))
            # node.contours.append([a[0], a[1], b[0], b[1], n2[0], n2[1]])




            # plt.figure()
            # plt.plot(line[:, 0], line[:, 1])
            # plt.plot([a[0], a[0]+n2[0]], [a[1], a[1]+n2[1]])
            # plt.plot([a[0], a[0]+n[0]], [a[1], a[1]+n[1]])
            # plt.scatter(line[0, 0], line[0, 1], label="line_a")
            # plt.scatter(line[1, 0], line[1, 1], label="line_b")
            # plt.scatter(a[0], a[1], label="a")
            # plt.scatter(b[0], b[1], label="b")
            # # plt.scatter(line[:,0], line[1,:])
            # # plt.scatter(np.array([a,b])
            # plt.axis("equal")
            #
            # plt.legend()
            # plt.show()
    #
    # if repeat:
    #     getPolyGridApproximation(residual_lines, tree_root, repeat=False)


if __name__ == "__main__":

    # filename = "../data/empty1x1.png"
    # filename = "../data/test_19.png"
    # filename = "../data/soundwaves/soundwaves2.png"
    # filename = "../data/apollo_cm/cm3.png"
    filename = "../data/vertical_line/vl1.png"
    # filename = "../data/wing1.png"
    img = cv2.imread(filename)
    img = cv2.flip(img, 0)
    closed_area=False
    depth = 10
    show_degrees=False#True
    show_normals=True
    show_cells=True
    # rescale_x = 1.
    # rescale_y = 1.
    rescale_x = False
    rescale_y = False


    imgray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    imgray = np.where(imgray > 30, 255, 0).astype(np.uint8)
    print(np.min(imgray))

    ret, thresh = cv2.threshold(imgray, 100, 255,0)

    contours, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
    img_contours = cv2.drawContours(img, contours, 0, (0,255,0), 1)
    print(len(contours))


    x,y = 0.,0.
    width, height = img.shape[1],img.shape[0]


    root = QuadTreeNode(0, x,y, width, height)

    # root.split()
    # root.children[0].split()
    # root.children[1].contours.append([1,2,3,4,5])
    # print(root.layerwise_subtree())
    # root.mergeSubtreeCells()
    # print(root.layerwise_subtree())

    # createSubTree(root, depth)

    poly_approx_simplified = [[[], []]]
    fig, ax = plt.subplots(figsize=(20,20))

    if not np.all(imgray == 255): # if the image is empty
        createSubTreeFromImage(root, depth, imgray)




        if closed_area:
            contours_ = contours
        else:
            contours_ = contours[1:]
        for c in contours_:
            poly = getPolyApproximation(c)
            col_lines = ax.add_collection(LineCollection(poly, colors="b", linewidths=1))

            poly_approx = snap_lines_to_grid(poly, x,y, width*2**(-depth), height*2**(-depth))
            # pprint(poly_approx)
            col_approx_lines = ax.add_collection(LineCollection(poly_approx, colors="darkred", linewidths=1))
            poly_approx_simplified = simplify_lines(poly_approx, x,y, width*2**(-depth), height*2**(-depth))
            col_approx_simple_lines = ax.add_collection(LineCollection(poly_approx_simplified, colors="green", linewidths=10))

            # pprint(poly_approx_simplified)
            getPolyGridApproximation(poly_approx_simplified, root)
            getPolyGridApproximation(poly_approx_simplified, root, first_iteration=False)

        root.mergeSubtreeCells()
    root.rescale(rescale_x, rescale_y)

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
                lines.append([c[:2], c[2:4]])
                if show_normals:
                    normals.append([c[:2], c[:2]+c[4:]])
                if show_degrees:
                    # ax.text((c[0]+c[2])/2, (c[1]+c[3])/2, f"{np.arcsin(c[4])*180/np.pi:.1f}", color="pink", fontsize=20)
                    ax.text((c[0]+c[2])/2, (c[1]+c[3])/2, f"{n.contour_angles[i]*180/np.pi:.1f}", color="pink", fontsize=20)

    # ax.text(lines[0][0][0], lines[0][0][1], "a")
    # ax.text(lines[0][1][0], lines[0][1][1], "b")
    # print(lines)
    approx_lines = ax.add_collection(LineCollection(lines, colors="red", linewidths=.75))
    if normals:
        normal_lines = ax.add_collection(LineCollection(normals[:], colors="green", linewidths=.5))

    ax.set_xlim(root.x-.1, root.width+.1)
    ax.set_ylim(root.y-.1, root.height+.1)

    if not rescale_x:
        rescale_x = imgray.shape[1]
    if not rescale_y:
        rescale_y = imgray.shape[0]
    ax.imshow(imgray, cmap="gray", extent=(x,rescale_x,y,rescale_y), origin="lower", aspect="auto")

    ax.axis("equal")
    # ax.imshow(img_contours)
    # ax.imshow(imgray)
    #
    #
    # cv2.imshow("contours", img_contours)

    plt.show()

    # save the contours
    print(np.array(poly_approx_simplified).shape)
    np.save(".".join(filename.split(".")[:-1])+".npy", poly_approx_simplified)

    root.to_csv(".".join(filename.split(".")[:-1])+".tree")
    # print(len(root.leafs()))

    new_root = QuadTreeNode.from_csv(".".join(filename.split(".")[:-1])+".tree")
    # print(new_root.layerwise_subtree())
    # pprint([n.contours for n in new_root.leafs()])

