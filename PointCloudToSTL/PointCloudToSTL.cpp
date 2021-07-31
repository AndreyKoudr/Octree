#include "../OOctree.h"
#include "../Triangles.h"
#include "../CloudPoints.h"

#include <iostream>


/**
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
This is needed 
to define sign of signed distance function to create a parametric geometry. Cells which are 
outside will have "1" in cell data; inside have "2". The algorithm is recursive : starting from an outside cell,
all its neighbours are checked for having a raster cell cell as neighbour. Marked are all neighbours which 
do not penetrate through raster cells. If the raster cells have a hole (e.g. due to insifficient number of
points, see POINTS_PER_CELL), the algorithm moves "inside" the body and the process fails (inside part is 
void). If this happen, increase POINTS_PER_CELL. At the same time, with POINTS_PER_CELL increase, small details
disappear.

  - Calculate distance field vectors at all nodes of raster cells (<I>addDistanceFields(()</I>).<br /> These are
vectors pointing from a node to the closest point on the surface. The code tries to build a LSQ plane over the
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
*/

int main(int argc, char* argv[])
{
/**
  VERY IMPORTANT : number of points per cell
  (1) defines mesh refinement (the smaller number the mesh is finer near the surface and 
the resulting surface with more detail) : octree cells are refined until a cell contains 
less or equal number of points
  (2) too low number produces "holes" in the resulting surface
  (3) too low number may produce "holes" in raster cells (octree cells with some points inside)
which leads to problems with in/out cell classification and void resulting geometry
*/
  int POINTS_PER_CELL = 150;

  printf("\n  Converter from point cloud XYZ file (text space separated) to parametric geometry STL.\n");
  printf("Call : >PointCloudToSTL XYZ_file_name points_per_cell\n");
  printf("e.g.\n");
  printf(">PointCloudToSTL shuttle.0.0005.xyz 200\n");
  printf("will generate shuttle.0.0005.200.stl\n");
  printf("  with side products as \n");
  printf("- raster.stl raster cells (those with points inside) saved as cell faces divided into two triangles\n");
  printf("- inside.stl cells inside geometry\n");
  printf("- outside.stl cells outside geometry\n");

  if (argc != 3) 
    return 1;

  std::string infilename = argv[1];
  POINTS_PER_CELL = atoi(argv[2]);
  LIMIT(POINTS_PER_CELL,3,10000);
  std::string outfilename = forceExtension(infilename,"") + "." + argv[2] + ".stl";

  printf("\nOutput : %s\n\n",outfilename.c_str());

  // these are cloud points
  CloudPoints<OREAL> points;

  printf("Reading file...");

  if (points.readXYZ(infilename))
  {
    printf(" %zd points read\n",points.size());

    // get body size
    std::pair<OVECTOR,OVECTOR> minmax = points.calcMinMax();

    // extend box size
    minmax.first *= 1.5;
    minmax.second *= 1.5;
    OVECTOR d = minmax.second - minmax.first;

    // we need something like 5 x 3 x 2 for background
    OREAL minsize = std::min<OREAL>(d.X,std::min<OREAL>(d.Y,d.Z)) * 0.1;

    LINT Imax = LINT(d.X / minsize);
    LINT Jmax = LINT(d.Y / minsize);
    LINT Kmax = LINT(d.Z / minsize);

    // last parameter is max octree level (max number of cell subdivisions during
    // refinement)
    OOctree octree(minmax.first,minmax.second,Imax,Jmax,Kmax,20);

    printf("Octree has %zd background cells, %zd x %zd x %zd\n",octree.numCells(),Imax,Jmax,Kmax);

    printf("Refining octree...");

    // refine cells with points; pointsInCells contains list of points inside every cell if any
    OPOINTS pointsInCells;
    octree.refineCellsForPointCloud(points,POINTS_PER_CELL,pointsInCells);

    printf("Refined octree has %zd cells\n",octree.numCells());

    // make raster cells (cells intersected by surface) that is which have some points inside
    OCELLS rastercells;
    octree.makeRasterCells(pointsInCells,rastercells);

    // save raster cells
    {
      TTriangles<OREAL> tris;
      if (octree.cellsToTriangles(rastercells,tris))
      {
        tris.saveSTL("raster.stl","Raster",true);
      }
    }

    printf("Classifying cells (in/out)...\n");

    octree.classifyCells(rastercells);

    // save outside cells (marked by "1" in octree cell data)
    {
      TTriangles<OREAL> tris;
      if (octree.cellsToTriangles(octree.cells(),1,tris))
      {
        tris.saveSTL("outside.stl","Inside",true);
      }
    }

    // save inside cells (marked by "2" in data)
    {
      TTriangles<OREAL> tris;
      if (octree.cellsToTriangles(octree.cells(),2,tris))
      {
        tris.saveSTL("inside.stl","Outside",true);
      }
    }

    printf("Calculating distance field vectors at %zd octree nodes...\n",octree.numNodes());

    // add distance fields to nodes
    octree.addDistanceFields(points,pointsInCells);

    printf("Calculating signed distance\n");

    // calculate signed siatnce at nodes
    octree.calcSignedDistance();

    printf("Polygonising surface with marching cubes...\n");

    // save generated geometry
    {
      TTriangles<OREAL> tris;
      octree.polygoniseDistanceFields(tris,rastercells);
      tris.saveSTL(outfilename,"Geometry",true);
    }
  }

  return 0;
}
