# uiowa-EQ-research Status
Building code to analyze differences in elevation post-earthquake

# Introduction
This software will process pre/post earthquake DEMs, and build an output suitable for differencing the changes in elevation. Calculating elevation differences between DEMs is quite an easy task, so why am I building code to do this?

# The Problem
If we were to simply difference pre/post event DEMS, our resulting DEM would not be useful. Issues arise due to the fact that earthquakes change the topography of the surrounding area. Natural formations (aka 'shapes') are translated laterally - which means the pre event DEM and post event DEM are not lined up. The differences calculated by simply doing (postDEM-preDEM) would also be taking into account the horizontal movement of land.

# The Solution
This code will be utilizing an iterative closest point (ICP) algorithm, along with shape matching techniques. Shapes between pre/post event DEMs will be matched according to their point cloud, and subsequently lined up.


