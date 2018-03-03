"""
ICP inputs: data - (i, j) array from gridPre (gridPre[i][j])
model - (i, j) array from gridPost (gridPost[i][j])
pNorm, method, and maxIter are defined in 

ICP outputs: TT - measuring translation around two points
TR - measuring rotation around two points

KDTreeSearcher - nearest neighbor search
USE POINT TO PLANE IN EVERYTHING ()

@author: Jack Biscupski
"""

import numpy as np
import scipy as sp

import h5py


# load in variables icpTest for easier testing

# using hdf5 file format to store individual grid cells
# can't load in entire gridPost (irregularly shaped grid)
f = h5py.File('icpTest.hdf5', 'r')
d_test1 = f['gridPre']['data1']
d_test2 = f['gridPre']['data2']
m_test1 = f['gridPost']['model1']
m_test2 = f['gridPost']['model2']


method = "Point2Plane"
maxIter = 500
pNorm = 10


# data will be 3x1, model will always have 3 rows but varying amount of columns
def icp(data, model, pNorm, method, maxIter):
    data_mean = np.reshape(np.mean(data, axis=1), (3, 1))  
    model_mean = np.reshape(np.mean(model, axis=1), (3, 1)) # question: pre/post grids will always have 3 rows, correct?
    # repeat matrices after reshaping them
    data = data - np.repeat(data_mean, len(data[0]), axis=1) # Is this always going to be 0? Every cell in gridPre is 3x1, so average along axis 1 is just the values
    model = model - np.repeat(data_mean, len(model[0]), axis=1)
    
    # set up Kd-Tree to be used in nearest neighbor search
    tree = sp.spatial.KDTree(np.transpose(model))
    
    # TO DO - need to find default bucket size for KDTree in MATLAB 
    # in sciPy, default leaf size is 10
    thres = 1e-5
    minIter = 10
    minRest = 9e99
    
    # Total Rotation
    TR = np.eye(np.size(model, axis=0)) # make sure size axes are right for TR and TT
    # Total Translation
    TT = np.zeros((np.size(model, axis=0), 1))
