"""
This is meant to be the 'driver' program for the ICP code. The main function of this program
is: load in pre/post event DEMs, process them and extract pertinent data, create grids for re-
sampling, populate the grids with resampled data, run the ICP code over the entire grid, and
finally create new, shape-matched post event DEM.

@author: Jack Biscupski
"""

from osgeo import gdal
import numpy as np
import time
from ICP import icp


#preEQ = "C:/Users/Jack Biscupski/Desktop/Research Materials/Code/WV01_13OCT01n14JAN24_10200100262F9000_1020010029seg1_2m_dem_1.tif"
#postEQ = "C:/Users/Jack Biscupski/Desktop/Research Materials/Code/WWV01_12MAY04n13APR01_102001001CAFAC00_1020010020seg1_2m_dem_1.tif"

preEQ = "C:/Users/Jack Biscupski/Desktop/Research Materials/Code/WV01_small.tif"
postEQ = "C:/Users/Jack Biscupski/Desktop/Research Materials/Code/WWV01_small.tif"

# initializing resampling variables
method = "Point2Plane"
slidingWindow = 10
windowSize = 2
postEQclamp = 30

# normally 500, set to 5 for testing
maxIter = 5
pNorm = 2


# loading and preprocessing data
ds_pre = gdal.Open(preEQ)
ds_post = gdal.Open(postEQ)
# convert geotiffs to arrays for manipulation
arr_pre, arr_post = ds_pre.ReadAsArray(), ds_post.ReadAsArray()
pre_transform, post_transform = ds_pre.GetGeoTransform(), ds_post.GetGeoTransform() # provides info about cell resolution and x/y extent
preShape, postShape = np.shape(arr_pre), np.shape(arr_post)                                  

# initialize x/y extents for ease of use
xExtentMinPre = pre_transform[0]
xExtentMaxPre = xExtentMinPre + (pre_transform[1]*ds_pre.RasterXSize)
yExtentMaxPre = pre_transform[3]
yExtentMinPre = yExtentMaxPre + (pre_transform[5]*ds_pre.RasterYSize)
pre = tuple([xExtentMinPre, xExtentMaxPre, yExtentMinPre, yExtentMaxPre])
xExtentMinPost = post_transform[0]
xExtentMaxPost = xExtentMinPost + (post_transform[1]*ds_post.RasterXSize)
yExtentMaxPost = post_transform[3]
yExtentMinPost = yExtentMaxPost + (post_transform[5]*ds_post.RasterYSize)
post = tuple([xExtentMinPost, xExtentMaxPost, yExtentMinPost, yExtentMaxPost])


# initiates the dimensions of preM/postM
preM = np.ndarray((3, preShape[0]*preShape[1]))
postM = np.ndarray((3, postShape[0]*postShape[1]))
                                         

# create meshgrids
preX, preY = np.meshgrid(np.arange(pre[0]+0.5, pre[1], pre_transform[1]), np.arange(pre[2]+0.5, pre[3], pre_transform[1]))
postX, postY = np.meshgrid(np.arange(post[0]+0.5, post[1], post_transform[1]), np.arange(post[2]+0.5, post[3], post_transform[1]))

# fill in values for pre/post matrices
# NumPy library automatically casts missing/Nan values as -9999
preM[0,:], postM[0,:] = np.reshape(np.transpose(preX), preShape[0]*preShape[1]), np.reshape(np.transpose(postX), postShape[0]*postShape[1])
preM[1,:], postM[1,:] = np.reshape(np.transpose(preY), preShape[0]*preShape[1]), np.reshape(np.transpose(postY), postShape[0]*postShape[1])
preM[2,:], postM[2,:] = np.reshape(np.transpose(arr_pre), preShape[0]*postShape[1]), np.reshape(np.transpose(arr_post), postShape[0]*postShape[1])  

# used to create max/min extents of resampled grid
minX = min(pre[0], post[0])
maxX = max(pre[1], post[1])
maxY = max(pre[3], post[3])
minY = min(pre[2], post[2])

# resolution of resampled grid
x_Range = np.linspace(minX, maxX+1, int((maxX-minX)/slidingWindow)+1)
y_Range = np.linspace(minY, maxY+1, int((maxY-minY)/slidingWindow)+1)

gridPre, gridPost, gridPreAligned = [], [], []

# time grid resampling - testing purposes only
print("Timing resampling...")
start = time.time()
# for [:,I_pre] and [:,J_pre], etc - we are going through each matrix prefixed (for example, preM/postM or temps)
# and picking out the values by starting at the beginning of I_ and J_ and saying: if True, keep this 3x1 element, if not, discard
for i in range(len(x_Range)):
    # initiate empty list for each entry
    gridPre.append([])
    gridPost.append([])
    gridPreAligned.append([])
    # resample x coordinates
    I_pre = np.logical_and((preM[0,:] >= x_Range[i]), (preM[0,:] < (x_Range[i] + windowSize)))    
    I_post = np.logical_and((postM[0,:] >= (x_Range[i] - postEQclamp)), (postM[0,:] < (x_Range[i] + windowSize + postEQclamp)))
    
    tempData_pre = preM[:,I_pre]
    tempData_post = postM[:,I_post]
    for j in range(len(y_Range)):
        gridPre[i].append([])
        gridPost[i].append([])
        gridPreAligned[i].append([])
        # resample y coordinates
        J_pre = np.logical_and((tempData_pre[1,:] >= y_Range[j]), (tempData_pre[1,:] < (y_Range[j] + windowSize)))
        J_post = np.logical_and((tempData_post[1,:] >= (y_Range[j] - postEQclamp)), (tempData_post[1,:] < (y_Range[j] + windowSize + postEQclamp)))
        # create each grid point/cell
        gridPre[i][j] = tempData_pre[:,J_pre]
        gridPost[i][j] = tempData_post[:,J_post]
# based on personal tests so far, takes anywhere from 37-42 mins        
end = time.time()    
print("Resampling elapsed time: {} seconds".format(end-start))

# create multiple zero-instantiated arrays with same shape as gridPreAligned
# plot Tx and Ty (base grid) for final output, showing TTx, TTy, TTz displacements
Tx, TTx, Ty, TTy, Tz, TTz = (np.zeros((len(x_Range)-1, len(y_Range)-1)) for i in range(6)) 
counter = -1
Tr = np.zeros((9, (len(x_Range)-1)*(len(y_Range)-1)))
for i in range(len(x_Range)-1):
    for j in range(len(y_Range)-1):
        data = gridPre[i][j]
        model = gridPost[i][j]
        # take advantage of empty containers evaluating to false
        if not data.size or not model.size:
            counter += 1
            gridPreAligned[i][j] = data
            continue
        TR, TT, dataOut = icp(data, model, pNorm, maxIter)
        delta = np.mean(dataOut-data, axis=1)
        TTx[i,j] = delta[0]
        TTy[i,j] = delta[1]
        TTz[i,j] = delta[2]
        Tx[i,j] = TT[0]
        Ty[i,j] = TT[1]
        Tz[i,j] = TT[2]
        Tr[:,counter] = TR.flatten()
        gridPreAligned[i][j] = dataOut
        counter += 1

x = np.zeros((len(y_Range)-1, len(x_Range)-1))
y = np.zeros((len(y_Range)-1, len(x_Range)-1))
for i in range(len(x_Range)-1):
    x[:,i] = (2*x_Range[i] + windowSize)/2
for i in range(len(y_Range)-1):
    y[i] = (2*y_Range[i] + windowSize)/2

x = np.transpose(x)
y = np.transpose(y)
n = (TTx**2 + TTy**2)**0.5
