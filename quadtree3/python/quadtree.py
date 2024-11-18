import matplotlib
import numpy as np
import cv2
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection, LineCollection

from typing import List

from pprint import pprint

class QuadTreeNode:

    def __init__(self, cellID, x, y, width, height, parent=None, contours=np.array([])):
        self.cellID = cellID
        self.x = x
        self.y = y
        self.width = width
        self.height = height

        self.contours = contours
        self.contour_angles = []

        self.parent = parent
        self.children: List["QuadTreeNode"] = None

    @property
    def level(self):
        return int(self.cellID >> 58)

    def split(self):
        self.children = []
        newHeight = self.height / 2
        newWidth = self.width / 2
        for i in range(4):
            cellID = (self.cellID & (2**58-1)) | ((self.level + 1) << 58) | (i << 2*self.level)
            # print(bin(cellID))
            x = self.x + (i&1) * newWidth
            y = self.y + (i&2)//2 * newHeight
            # print(x,y)
            self.children.append(QuadTreeNode(cellID, x, y, newWidth, newHeight, parent=self, contours=[]))

    def leafs(self) -> List["QuadTreeNode"]:
        if self.children is None:
            return [self]
        leafs = []
        for c in self.children:
            leafs.extend(c.leafs())
        return leafs

    @property
    def is_mergeble(self):
        return self.children is None and not self.contours

    def mergeSubtreeCells(self):
        is_mergeble = True
        if self.children is None:
            return
        for c in self.children:
            c.mergeSubtreeCells()
            if not c.is_mergeble:
                is_mergeble = False
        if is_mergeble:
            self.children = None

    def __lt__(self, other: "QuadTreeNode"):
        if self.x < other.x:
            return True
        elif self.x == other.x:
            return self.y < other.y
        return False

    def __repr__(self):
        return f"QuadTreeNode(x={self.x:.2f}, y={self.y}, width={self.width}, height={self.height}, level={self.level})"

    def layerwise_subtree(self):
        nodes = [self]
        i = 0
        while i < len(nodes):
            n = nodes[i]
            if n.children is not None:
                nodes.extend(n.children)
            i += 1
        return nodes
    def to_csv(self, filename:str):
        nodes = self.layerwise_subtree()
        with open(filename, "w") as f:
            f.write(f"{len(nodes)} {self.width:.5f} {self.height:.5f}\n")
            f.writelines([f"{n.cellID:019d} {1 if n.children is None else 0} {len(n.contours)}\n" + (f"{' '.join([' '.join([f'{c_:.5f}' for c_ in c]) for c in n.contours])}\n" if n.contours else "") for n in nodes])
            # f.writelines([f"{n.cellID:020d} {1 if n.children is None else 0} {len(n.contours)}{' ' if n.contours else ''}{' '.join([' '.join([f'{c_:0=q7.5f}' for c_ in c]) for c in n.contours])}\n" for n in nodes])
        print(f"stored QuadTree in {filename}")

    @classmethod
    def from_csv(cls, filename):

        with open(filename, "r") as f:
            line1 = f.readline().split(" ")
            width = float(line1[1])
            height = float(line1[2])
            line2 = f.readline().split(" ")
            root = cls(int(line2[0]), 0, 0, width, height, contours=[])
            root.children = []
            for i in range(int(line1[0])-1):
                line = f.readline()
                l = line.split(" ")
                cellID = int(l[0])
                isLeaf = bool(int(l[1]))
                nContours = int(l[2])
                # print(cellID, isLeaf, nContours)

                n = root
                for level in range(0,(cellID >> 58) - 1):
                    n = n.children[(cellID >> 2*level) & 3]

                n.children.append(QuadTreeNode(0,0,0,0,0, parent=n))
                n = n.children[-1]
                if not isLeaf:
                    n.children = []
                else:
                    if nContours > 0:
                        contours = f.readline().split(" ")
                        n.contours = np.array(contours, dtype=np.float64).reshape((nContours, 6))
                    else:
                        n.contours = np.array([])

        return root

    def rescale(self, rescale_x=False, rescale_y=False):
        width = self.width
        height = self.height
        if rescale_x:
            scaling_x = rescale_x/width
        else:
            scaling_x = 1
        if rescale_y:
            scaling_y = rescale_y/height
        else:
            scaling_y = 1
        nodes = self.layerwise_subtree()
        for n in nodes:
            n.x *= scaling_x
            n.y *= scaling_y
            n.width *= scaling_x
            n.height *= scaling_y

            for i in range(len(n.contours)):
                n.contours[i] = np.array(n.contours[i], dtype=np.float64)
                n.contours[i][0::2] *= scaling_x
                n.contours[i][1::2] *= scaling_y



def find_cell(root:QuadTreeNode, x,y) -> QuadTreeNode:
    if root.children is None:
        return root
    i = 0
    if x > root.x + root.width/2:
        i += 1
    if y > root.y + root.height/2:
        i += 2
    return find_cell(root.children[i], x,y)

def createSubTree(root, depth):
    if depth == root.level:
        return

    root.split()
    for i in range(len(root.children)):
        createSubTree(root.children[i], depth)

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



def getRectangles(root: QuadTreeNode):

    if root.children is None:
        return [Rectangle((root.x, root.y), root.width, root.height)]
    else:
        rects = []
        for i in range(len(root.children)):
            rects.extend(getRectangles(root.children[i]))

        return rects

def getPolyApproximation(contour):
    epsilon = 0.005 * cv2.arcLength(contour, True)  # Tolerance factor
    approx_contour = cv2.approxPolyDP(contour, epsilon, True)  # Approximate contour
    pprint(approx_contour)

    lines = np.array([[approx_contour[i,0,:], approx_contour[i+1,0,:]] for i in range(len(approx_contour)-1)] + [[approx_contour[-1,0,:], approx_contour[0,0,:]]])

    return lines

def line_intersection(line1, line2):
    # based on https://paulbourke.net/geometry/pointlineplane/, https://www.jeffreythompson.org/collision-detection/line-line.php
    (x1,y1), (x2,y2) = line1
    (x3,y3), (x4,y4) = line2

    # denominator
    d = (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1)
    # d = (x1-x1)*(y3-y4) - (y1-y2) * (x3-x4)

    if d == 0: # the lines are parallel
        # even if the lines are identical, return None since then it intersects with the perpendicular borders
        return None

    # calculate the intersection point
    uA = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / d
    uB = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / d

    if not (0 <= uA <= 1) or not (0 <= uB <= 1):
        return None

    return np.array([x1 + (uA * (x2-x1)), y1 + (uA * (y2-y1))])

def snap_point_to_grid(point, x,y,dx,dy, snap_x=True, snap_y=True):
    px,py = point

    px_ = px - x
    py_ = py - y

    if snap_x:
        if px_ - np.floor(px_/dx) * dx < np.ceil(px_/dx) * dx - px_:
            dist_x = px_ - np.floor(px_/dx) * dx
            new_x = px - dist_x
        else:
            dist_x = np.ceil(px_/dx) * dx - px_
            new_x = px + dist_x
    else:
        new_x = px
        dist_x = np.inf

    if snap_y:
        if py_ - np.floor(py_/dy) * dy < np.ceil(py_/dy) * dy - py_:
            dist_y = py_ - np.floor(py_/dy) * dy
            new_y = py - dist_y
        else:
            dist_y = np.ceil(py_/dy) * dy - py_
            new_y = py + dist_y
    else:
        new_y = py
        dist_y = np.inf

    if dist_x < dist_y:
        return np.array([new_x, py])
    else:
        return np.array([px, new_y])

    # dist_py = min(py_ - np.floor(py_/dy) * dy, py_ - np.ceil(py_/dx) * dy )

def snap_lines_to_grid(lines, x,y,dx,dy):
    new_lines = []
    for line in lines:
        new_lines.append([snap_point_to_grid(line[0],x,y,dx,dy), snap_point_to_grid(line[1],x,y,dx,dy)])

    new_lines = np.array(new_lines)
    pprint(new_lines)
    # remove redundant lines (a->b and b->a)
    new_lines_ = []
    skip_next = False
    for i in range(len(new_lines)-1):
        if skip_next:
            skip_next = False
            continue
        if not np.allclose(new_lines[i], new_lines[i+1, ::-1]):
            new_lines_.append(new_lines[i])
        else:
            skip_next = True

    if len(new_lines_) > 0:
        if np.allclose(new_lines_[-1], new_lines[-1,::-1]):
            new_lines_ = new_lines_[:-1]
        else:
            new_lines_.append(new_lines[-1])


    return np.array(new_lines_)

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
                    print(a, b)
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
    filename = "../data/test_18.png"
    img = cv2.imread(filename)
    closed_area=False
    depth = 9
    show_degrees=True
    show_normals=True
    show_cells=True
    # rescale_x = 1.
    # rescale_y = 1.
    rescale_x = False
    rescale_y = False


    imgray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    imgray = np.where(imgray > 10, 255, 0).astype(np.uint8)
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
    createSubTreeFromImage(root, depth, imgray)



    fig, ax = plt.subplots(figsize=(20,20))

    if closed_area:
        contours_ = contours
    else:
        contours_ = contours[1:]
    for c in contours_:
        poly = getPolyApproximation(c)
        col_lines = ax.add_collection(LineCollection(poly, colors="b", linewidths=1))

        poly_approx = snap_lines_to_grid(poly, x,y, width*2**(-depth), height*2**(-depth))
        pprint(poly_approx)
        col_approx_lines = ax.add_collection(LineCollection(poly_approx, colors="darkred", linewidths=1))
        poly_approx_simplified = simplify_lines(poly_approx, x,y, width*2**(-depth), height*2**(-depth))
        col_approx_simple_lines = ax.add_collection(LineCollection(poly_approx_simplified, colors="green", linewidths=1))

        pprint(poly_approx_simplified)
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
        if n.contours:
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

    root.to_csv(".".join(filename.split(".")[:-1])+".tree")

    new_root = QuadTreeNode.from_csv(".".join(filename.split(".")[:-1])+".tree")
    print(new_root.layerwise_subtree())
    pprint([n.contours for n in new_root.leafs()])

