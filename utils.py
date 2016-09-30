#!/usr/bin/python
import matplotlib.path as mplPath
import numpy as np


def inside_polygon(points, path):
    """Check if the points are inside a polygon 

    Keyword arguments:
    points -- numpy array containing the coordinates of the points
    path -- matplotlib.path describing the vertex of the polygon
    """
    return path.contains_points(points)

def rotate2D(xy, yaw=0):
    print xy
    return np.array([xy[0]*np.cos(yaw) - xy[1]*np.sin(yaw),       # X coordinate rotation of yaw
                     xy[0]*np.sin(yaw) - xy[1]*np.cos(yaw)])      # Y coordinate rotation of yaw


def evaluate_sss_path(nav_status, width, length):
    """Build a matplotlib.path describing the rectangular fov of a 
    side-scan sonar.

    Keyword arguments:
    nav_status -- numpy array describing the navigation status of the auv [x, y, yaw]
    width -- swath width of the side-scan sonar
    length -- swath length of the side-scan sonar
    """
    vertexes = np.array([[nav_status[0] + width/2, nav_status[1] + length/2], # UR Vertex
                         [nav_status[0] - width/2, nav_status[1] + length/2], # UL Vertex
                         [nav_status[0] - width/2, nav_status[1] - length/2], # BL Vertex
                         [nav_status[0] + width/2, nav_status[1] - length/2]   # BR Vertex
                       ])
    r_vertexes = np.apply_along_axis(rotate2D, 1, vertexes, nav_status[2])
    print vertexes
    print vertexes.shape
    print r_vertexes
    return mplPath.Path(r_vertexes)
