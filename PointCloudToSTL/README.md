
  PointCloudToSTL
  ===============
  
  Reconstruction of parametric surface from point cloud
  
  - Purpose : reconstruct a parametric surface given unordered point cloud (only X,Y,Z coordinates).
 
  - Input (<I>points.readXYZ()</I>).<br />
Read simple point cloud in .xyz file (three X,Y,Z coordinates in text, separated by spaces not tabs).
The point cloud must be dense enough, some 0.5 .. 10 millions of points uniformly distributes over 
the model surface. Any area of missing points may be considered as a hole which may not let the code
define three groups of octree cells : raster (octree cells which have at least one could point inside),
inside (cells whole inside the body), and outside (whole cells outside the body).

  - Octree refinement (<I>refineCellsForPointCloud()</I>).<br />
Starting from uniform octree 
background (the largest cells) octree cells are subdivided into 8 sub-cells until the number of points
in a cell becomes less than POINTS_PER_CELL. <I>pointsInCells</I> contains list of points inside 
every cell if any.

  - Mark raster cells (<I>makeRasterCells()</I>).<br /> 
Raster cells are those which contain at least one point inside.
They are marked in octree by data value of "3" and the whole list of them <I>rastercells</I> is generated.
  - Classify cells of the whole octree as inside and outside of surface (<I>classifyCells()</I>).<br /> 
This is needed to define sign of signed distance function to create a parametric geometry. Cells which are 
outside will have "1" in cell data; inside have "2". The algorithm is recursive : starting from an outside cell,
all its neighbours are checked for having a raster cell cell as neighbour. Marked are all neighbours which 
do not penetrate through raster cells. If the raster cells have a hole (e.g. due to insifficient number of
points, see POINTS_PER_CELL), the algorithm moves "inside" the body and the process fails (inside part is 
void). If this happen, increase POINTS_PER_CELL. At the same time, with POINTS_PER_CELL increase, small details
disappear.
  - Calculate distance field vectors at all nodes of raster cells (<I>addDistanceFields(()</I>).<br /> 
These are vectors pointing from a node to the closest point on the surface. The code tries to build a LSQ plane over the
points in a cell and project node positions onto this plane; if this fails, the plane is built over three 
first points inside the cell; if this fails as well, just a minimum distance to a point is calculated. The minimum 
distance for a node among all neighbouring cells is selected as distance filed vector.

  - <I>calcSignedDistance()</I> calculates calculates signed distance function values at 
raster cell nodes. Signed distnace function is the length of distance field vector with a sign 
taken from the inside/ouside step : values of the function for inside nodes are negative, outside are positive.
  - And the last step : triangulate signed distance function values at raster cell nodes by marching cubes
  <I>polygoniseDistanceFields()</I>. The algorithm is modified from a traditional marching cubes (taken from 
http://paulbourke.net/geometry/polygonise/). The problem is to provide continuous approximation for
the signed distance function across cells of different sizes. This is done by using conforming finite
elements (https://github.com/AndreyKoudr/FiniteElements). 
Please keep in mind that a perfect output for this process in NOT GUARANTEED because a point cloud does not 
contain ALL data about the surface; we only try to <B>reconstruct</B> it.

  Side products
  =============
  Sets of inner and outer cells can be considered as non-conformal meshes for conforming finite element 
calculations.
