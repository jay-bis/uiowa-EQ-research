"""
ICP inputs: data - (i, j) array from gridPre (gridPre[i][j])
model - (i, j) array from gridPost (gridPost[i][j])
pNorm, method, and maxIter are defined in 

@author: Jack Biscupski
"""

import numpy as np
from scipy import spatial
from pyntcloud import PyntCloud
import pandas as pd
import math as m

import h5py



# load in variables from icpTest for easier testing

# using hdf5 file format to store individual grid cells
# can't load in entire gridPost (irregularly shaped grid)
# delete all of this and get rid of h5py import for final version
f = h5py.File('icpTest.hdf5', 'r')
d_test1 = np.asarray(f['data1'])                     ##[0][0]
m_test1 = np.asarray(f['model1'])                    ##[0][0]
d_test2 = np.asarray(f['data2'])                     ## [45][120]
m_test2 = np.asarray(f['model2'])                    ## [45][120]
d_test3 = np.asarray(f['data3'])                     ## [217][89]
m_test3 = np.asarray(f['model3'])                    ## [217][89]
f.close()

# normally 500, set to 10 for easier testing purposes
maxIter = 10
pNorm = 2



def icp(data, model, pNorm, maxIter):
    data_mean = np.reshape(np.mean(data, axis=1), (3, 1))  
    model_mean = np.reshape(np.mean(model, axis=1), (3, 1)) 
    # repeat matrices after reshaping them
    data = data - np.repeat(data_mean, len(data[0]), axis=1) 
    model = model - np.repeat(data_mean, len(model[0]), axis=1) 
    
    
    # set up Kd-Tree to be used in nearest neighbor search
    tree = spatial.KDTree(np.transpose(model))

    thres = 1e-5
    minIter = 10
    res = 9e99
    
    # Total Rotation
    TR = np.eye(np.size(model, axis=0)) # make sure size axes are right for TR and TT
    # Total Translation
    TT = np.zeros((np.size(model, axis=0), 1))
    
    data01 = data
    data02 = data
    minRes = res
    # using this to try to figure out why residual increases seemingly exponentially through each iter in icp
    #origData = data
    
    normals = computeNormals(model)
    for i in range(maxIter):
        # second element in I tuple is the 'location' we will use to slice normals and model. Occasionally we will get more
        # than one nearest neighbor, so the first result from the array is taken. -1 because of 0-based indexing
        I = tree.query(np.transpose(data))[1][0]-1
        # @ preferred for matrix multiplication over dot()
        resI = abs((normals[:,I]) @ (data-np.reshape(model[:,I], (3, 1))))
        wghs = computeWeights(resI, pNorm)[0]
        # use shapes[0] for R and shapes[1] for T
        shapes = RigidShapeMatching(data, model[:,I], normals[:,I], wghs)
        # with this if/else block, I'm trying to replicate what is happening in MATLAB: if the lstsq method needs to be
        # used, that means matrix is singular, and in MATLAB this returns essentially all NaN's for values. Way the code was
        # set up before, rigidShapeMatching would occasionally throw a ValueError for trying to use too large of a value
        # in a trigonometric function
# alternative:
#        if shapes[2]:
#            data = data*np.nan
#            TR = np.nan*TR
#            TT = np.nan*TT
#        else:
        data = shapes[0].dot(data) + shapes[1]
        TR = shapes[0].dot(TR)
        TT = shapes[0].dot(TT) + shapes[1]
        res = m.sqrt(max(sum((data-data01)**2)))
        if minRes > res:
            minRes = res
            data02 = data
        if i > minIter and res < thres:
            data02 = data
            break        
        data01 = data # originally: data01 = data
        #print("Residual: {}".format(res))
    data = data02
    data = data + data_mean
    print('1 iteration complete. Last residual: {}'.format(res))
    return (TR, TT, data)
            

# preload array 'template' so we're not creating this array every call
masterOrig = np.array([[0], [0], [1]])

"""
    INPUTS: single cell-matrix from gridPost (model)
    OUTPUTS: normals of gridPost matrix

    @author: Jack Biscupski
"""

def computeNormals(model):
    cloudable = pd.DataFrame(data=np.transpose(model), columns=['x', 'y', 'z'])
    # calculates a pointcloud of the input model using a pandas DataFrame
    cloud = PyntCloud(cloudable)
    # use neighbors to get normals from pointcloud
    neighbors = cloud.get_neighbors(k=10)
    cloud.add_scalar_field('normals', k_neighbors=neighbors)
    # extract normals from the altered DataFrame
    normals = np.transpose(np.asarray(cloudable.loc[:, 'nx(10)':'nz(10)']))
    master = np.repeat(masterOrig, len(normals[0]), axis=1)
    # taking the dot product column-wise; the 'ij,ij->' notation is saying:
    # take dot product of ith row, jth column of normals and ith row, jth column of master
    # then, create boolean mask array for normal comparison
    I = np.einsum('ij,ij->j', normals, master) < 0
    # flip all values in column if I is true at that column (dot prod < 0)
    normals[:,I] = -normals[:,I]
    
    return(normals)


"""
INPUTS: dist - magnitude of closest points (float)
        p - pNorm (integer)
OUTPUTS: weights 

@author: Jack Biscupski
"""

def computeWeights(dist, p):
    reg = 1e-8
    return(p/(dist**(2-p) + reg))  
    

"""
INPUTS: X, Y, N, and weights
OUTPUTS: R and T

@author: Jack Biscupski
"""
def RigidShapeMatching(X, Y, N, wghs):
    X = np.transpose(X)
    Y = np.transpose(Y)
    N = np.transpose(N)
    Cn = np.reshape(np.append((np.cross(X, N)), N), (1, 6))
    A = (np.tile(wghs, (6, 1))*np.transpose(Cn)).dot(Cn)
    # should these elements be put into a list? idk
    b = []
    for i in range(0, 6):
        b.append(sum(np.sum(((Y-X)*np.tile(Cn[:,i], (1, 3))*N), axis=1)))
    try:
        X = np.linalg.solve(A, b)
        #flag = False
    except:
        #flag = True
        # rcond set here to decompose and reconstruct A
        X = np.linalg.lstsq(A, b, rcond=1e-15)[0]
    cx = m.cos(X[1])
    cy = m.cos(X[2])
    cz = m.cos(X[3])
    sx = m.sin(X[1])
    sy = m.sin(X[2])
    sz = m.sin(X[3])
    R = np.asarray([[cy*cz, sx*sy*cz-cx*sz, cx*sy*cz+sx*sz], [cy*sz, cx*cz+sx*sy*sz, cx*sy*sz-sx*cz], [-sy, sx*sy, cx*cy]])
    T = np.transpose(np.atleast_2d(X[3:6]))
    # tuple of R and T
    return (R, T)
