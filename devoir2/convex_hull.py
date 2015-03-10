# http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain

import numpy as np

def convex_hull(points):
    """Computes the convex hull of a set of 2D points."""

    if len(points) <= 1:
            return points

    points = np.sort(points)

    # 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
 
    # Build lower hull 
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
 
    # Build upper hull
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
 
    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the beginning of the other list.
    result = np.array(lower[:-1] + upper[:-1], dtype=[("x", float), ("y", float), ('id', int)])
    return result
 
 
# Example: convex hull of a 10-by-10 grid.
#assert convex_hull([(i/10, i%10) for i in range(100)]) == [(0, 0), (9, 0), (9, 9), (0, 9)]


def find_best_permutation(points, distance_matrix):
    """find the best permutation possible, taking into account the convex hull"""
    # only allow permutation of points not belonging to the convex hull
    pass
