import numpy as np
import scipy as sp

"""
ICP inputs: data - (i, j) array from gridPre (gridPre[i][j])
model - (i, j) array from gridPost (gridPost[i][j])
pNorm, method, and maxIter are defined in 

ICP outputs: TT - measuring translation around two points
TR - measuring rotation around two points

KDTreeSearcher - nearest neighbor search
USE POINT TO PLANE IN EVERYTHING ()
"""

method = "Point2Plane"
# dummy data - for testing

# data will be 3x1, model will always have 3 rows but varying amount of columns
def icp(data, model, pNorm, method, maxIter):
    data_mean = np.mean(data, axis=1)
    model_mean = np.mean(model, axis=1)
