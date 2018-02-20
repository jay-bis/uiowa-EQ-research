from osgeo import gdal
import matplotlib.pyplot as plt



DEM1 = "C:/Users/Jack Biscupski/Desktop/Research Materials/Code/WV01_13OCT01n14JAN24_10200100262F9000_1020010029seg1_2m_dem_1.tif"
DEM2 = "C:/Users/Jack Biscupski/Desktop/Research Materials/Code/WWV01_12MAY04n13APR01_102001001CAFAC00_1020010020seg1_2m_dem_1.tif"

# to set how different graphics/plots display in spyder, go to:
# tools -> preferences -> Ipython console -> Graphics -> Backend
# Automatic makes new window with plot pop up


# two geo coordinate systems: spherical (lat and lon)
# what we will be using, UTM: in meters, instead of degrees

####################################################################

ds1 = gdal.Open(DEM1)
ds2 = gdal.Open(DEM2)
data1 = ds1.ReadAsArray()
data2 = ds2.ReadAsArray()



plt.figure()
plt.suptitle('2 meter x 2 meter pixel resolution')
plt.title('DEM 1')
plt.imshow(data1, cmap='jet', vmin=500, vmax=1076)
plt.colorbar()

plt.figure()
plt.suptitle('2 meter x 2 meter pixel resolution')
plt.title('DEM 2')
plt.imshow(data2, cmap='jet', vmin=500, vmax=1076)
plt.colorbar()

plt.figure()
plt.suptitle('2 meter x 2 meter pixel resolution')
plt.title('DEM Differential')
plt.imshow(data1-data2, vmin=-10, vmax=10)
plt.colorbar()
