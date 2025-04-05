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
            f.writelines([f"{n.cellID:019d} {1 if n.children is None else 0} {len(n.contours)}\n" + (f"{' '.join([' '.join([f'{c_:.5f}' for c_ in c]) for c in n.contours])}\n" if np.array(n.contours).size > 0 else "") for n in nodes])
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
                n.contours[i][0:3:2] *= scaling_x
                n.contours[i][1:4:2] *= scaling_y



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

def getRectangles(root: QuadTreeNode):

    if root.children is None:
        return [Rectangle((root.x, root.y), root.width, root.height)]
    else:
        rects = []
        for i in range(len(root.children)):
            rects.extend(getRectangles(root.children[i]))

        return rects

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
    # pprint(new_lines)
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
