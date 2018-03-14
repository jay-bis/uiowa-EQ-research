from pyntcloud import PyntCloud
import pandas as pd
import numpy as np
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
