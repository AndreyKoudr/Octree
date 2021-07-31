/*
BSD 2-Clause License

Copyright (c) 2020, Andrey Kudryavtsev (andrewkoudr@hotmail.com)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once
                      
#include "Types.h"
#include "OCell.h"
#include "ONode.h"
#include "OBackground.h"
#include "Vector.h"
#include "Triangles.h"

#include <map>
#include <set>
#include <vector>
#include <array>
#include <limits>

using namespace std;

/**
  Octree
  ======

  Background.<br /><br />

  Background is a set of largest cells obtained by uniform division of the whole cuboid region.
  Background cells are those of level 0. Every cell of level 0 can be subdivided (refined)
to make a hierarchy of cells inside a background cell. 

  Search.<br /><br />

  A search of every 3D point inside octree is two-step :
  (1) which background cell? 
  (2) starting from this background cell, search within a hierarchy
Both operations are fast. This is the way an octree provides a kind of indexing for 3D points.<br />
  Use <I>findCellCentreAndCheck()</I> to find a cell for a point position.

  Cells and nodes.<br /><br />

  Octree is a collection of cells (<I>OCELLS cells_;</I>) and nodes (ONODES nodes_;) both wholly 
defined by their single integer coordinates (cell centre for cells) <I>IPosition</I>. 
Position of cell also uniquely defines its level (number of cell subdivisions from background 
cell) (use <I>level()</I> for that) and its size (<I>intCellSize()</I>).<br />
  
  Each cell and every node carries data in OCELL_DATA and ONODE_DATA type variables. 

  Cell neighbours.<br /><br />

  No information about neighbours is in memory inside a cell : all info is generated on the 
fly by search of neighbour cells in <I>OCELLS</I> map. This makes the code very simple, reliable and 
saves lots of memory.

  Refine/derefine cells.<br /><br />

  Cells can be refined and de-refined. Max refinement level (for the the smallest posiible cell) is 
defined by <I>maxLevel</I> in the background. If you specify maxLevel as 20, it means that the smallest 
cell may have size of 

  background_cell_size / 2^20

that is, ~0.001mm for background size 1m.

  Call <I>refineCell()</I> to subidvide a cell into 8 sub-cells; it is not ALWAYS possible because all newly generated
cells and their neighbours must satisfy the <B>balanced</B> ocree condition : difference in levels of any two 
neigbour cells must not be higher than 1.<br />
  Call <I>derefineCell()</I> on a parent cell to delete all its sub-cells. <I>isLeaf()</I> function tells if
a cell has children or not (is a leaf). DE-refinement may not be successful as well due to a failing 
balanced condition.
*/


/** Comparison function to compare two integer positions */
class IPosCompare {
public:
  bool operator()(const IPosition& lhs, const IPosition& rhs) const
  {
    if (lhs.X < rhs.X)
    {
      return true;
    } else if (lhs.X > rhs.X)
    {
      return false;
    } else
    {
      if (lhs.Y < rhs.Y)
      {
        return true;
      } else if (lhs.Y > rhs.Y)
      {
        return false;
      } else
      {
        if (lhs.Z < rhs.Z)
        {
          return true;
        } else if (lhs.Z > rhs.Z)
        {
          return false;
        } else
        {
          return false;
        }
      }
    }
  }
};

/** Compare integer positions */
bool IPosComp(const IPosition& lhs, const IPosition& rhs);

/** Octree real type */
#define OREAL double

/** Octree vector type */
#define OVECTOR TVector<OREAL>

/** Data carried inside a cell. First value is cell type (raster, inner, outer), second -
  rasterised triangle number for raster cells. You can change this to your needs. */
#define OCELL_DATA std::pair<char,LINT>

/** Data carried inside a node. You can change this to your needs. This data is for distance field
  vector (shortest vector from a node to the surface) */
#define ONODE_DATA TVector<OREAL>

/** Two central data : list of cells and list of nodes defined by their integer positions. 

  IMPORTANT. If to change std::map to the Google btree map (cpp-btree-1.0.1), any map would take
  TWICE less memory with the same speed.
*/
#define OCELLS std::map<IPosition,OCell<OCELL_DATA>,IPosCompare>
#define ONODES std::map<IPosition,ONode<ONODE_DATA>,IPosCompare>

/** Temp data : list of indexes on integer coordinate. */
#define OPOSITIONS std::set<IPosition,IPosCompare>
#define OPOINTS std::map<IPosition,std::set<LINT>,IPosCompare>


class OOctree : public OBackground<OREAL> {
public:
  
  /** Default constructor. */
  OOctree() = delete;

  /** Constructor. 
    /p boxMin min box coordinate
    /p boxMax max box coordinate
    /p Imax number of cells in I direction
    /p Jmax number of cells in J direction
    /p Kmax number of cells in K direction
    /p maxLevel max octree level PLUS 1 (for 8 coordinates of finest cell nodes). OOctree itself is level 0.
  */
  OOctree(const OVECTOR &min, const OVECTOR &max, 
    const LINT Imax, const LINT Jmax, const LINT Kmax, const int maxOctreeLevel);

  //===== Statistics ===========================================================

  /** Number of cells. */
  size_t numCells();

  /** Number of nodes. */
  size_t numNodes();

  /** Access to cells. */
  OCELLS& cells();

  /** Access to nodes. */
  ONODES& nodes();

  //===== Cell sizes =========================================================

  /** Get integer size for cell of /p level. */
  inline LINT intCellSize(const int level);

  /** Get real sizes for cell of /p level. */
  inline OVECTOR realCellSizes(const int level);

  //===== Coordinates ==========================================================

  /** Convert int to real coordinate. */
  inline OVECTOR intToRealCoord(const IPosition& coord);

  /** Convert real to int coordinate. May return odd-valued integer coordinate. 
    Not quite clear why. */
  inline IPosition realToIntCoord(const OVECTOR& coord);

  /** Convert to closest cell (of size level) node coordinate */
  IPosition realToIntCoord(const OVECTOR &coord, const int level);

  /** Get 8 integer coordinates of cell corners. */
  std::array<IPosition,8> cell8Coordinates(const IPosition& ocentre);

  /** Get 8 real coordinates of cell corners. */
  std::array<OVECTOR,8> cell8RealCoordinates(const IPosition& ocentre);

  /** Get 27 (3 layers by 9 nodes) child coordinates. */
  std::array<IPosition,27> childrenNodeCoordinates(const IPosition& ocentre);

  /** Get min/max coordinates of 8 cell corners. */
  std::pair<IPosition,IPosition> minMaxCoordinates(const IPosition& ocentre);

  /** Get min/max real coordinates for hex nodes. Box extended by tolerance. */
  std::pair<OVECTOR,OVECTOR> minMaxCoordinates(const IPosition& ocentre, const OREAL tolerance);

  /** Get min/max real coordinates for hex nodes. Box extended by extension, like 1.1 */
  std::pair<OVECTOR,OVECTOR> minMaxCoordinatesExtended(const IPosition& ocentre, const OREAL extension);

  /** Min/max for all cells with data value. */
  std::pair<OVECTOR,OVECTOR> minmax(const OCELLS &cells, const char value);

  /** Extract node from cells with specified value in data. */
  void extractNodes(const OCELLS &cells, const char value, OPOSITIONS &nodes);

  /** Go through all cells marked by value to register all faces which have cells with another
    value on the other side, facecells being centres of cells with these faces,
    tris being numbers of closest triangles of body surface. */
  void extractBoundaryFaces(const OCELLS &cells, const char value, const char rastervalue,
    std::vector<std::array<IPosition,4>> &faces, std::vector<IPosition> &facecells, 
    std::vector<LINT> &tris, std::map<IPosition,std::set<LINT>,IPosCompare>& postotri);


  //===== Children/parents =====================================================

  /** 
    Level from central coordinate. 
    If position is an cell centre, full cell size is intCellSize().
    If position is a node coordinate, distance to the same level node is 
      intCellSize(level(coord)) / 2.
  */
  int level(const IPosition& position);

  /** Leaf does not have any children. */
  bool isLeaf(const IPosition& ocentre);

  /** Get cell parent centre. Use subParentDirs to find Cell. */
  IPosition parentCentre(const IPosition& ocentre, const OCell<OCELL_DATA> &cell);

  /** Get 8 children centres. Use subParentDirs to find Cells. */
  std::array<IPosition,8> childrenCentres(const IPosition& ocentre);

  /** Get 1 child centre. Use subParentDirs to find Cells. */
  IPosition childCentre(const IPosition& ocentre, const int subIndex);

  /** Convert list of centre coordinates into pointers to cells which were found. */
  size_t centresToCells(const IPosition centres[], const size_t numCentres, std::vector<OCell<OCELL_DATA> *>& cells);

  /** Get children. Returns true if at least one found. */
  bool getChildren(const IPosition& ocentre, std::array<IPosition,8> &centres, 
    std::vector<OCell<OCELL_DATA> *> &cells);

  /** Get smaller 24 neighbours across faces. */
  bool getSmallerNeighbours(const IPosition& ocentre, 
  std::array<std::array<IPosition,4>,6> &centres, std::array<std::array<OCell<OCELL_DATA> *,4>,6> &cells);


  //===== Neighbours ===========================================================

  /** Get cell neighbours for 6 faces - same size if present, otherwise one 
    level higher. A neighbour can be same size or bigger.

 _______ _______
|       |       |
| Cell  | Neigh |
|       |       |
|_______|_______|      

 _______ _______
|   |   |       |
|___|___|   N   |
|   | O |       |
|___|___|_______|      

 _________ _________
|         |    |    |
|         |    |    |
|    O    |----N----|  (parent)
|         |    |    |
|_________|____|____|      

*/

  /** Get 6 same-size neightbour centres. */
  std::array<IPosition,6> faceNeiCentres(const IPosition& ocentre);

  /** Get 6 same-size neighbours. */
  bool getSameSizeNeighbours(const IPosition& ocentre, std::array<IPosition,6>& centres, 
    std::array<OCell<OCELL_DATA> *,6> &cells);

  /** Get max 26 neighbours of same size. */
  void getAllSameSizeNeighbours(const IPosition& ocentre, std::vector<IPosition>& centres);

  /** Get two layers of same size neighbours */
  void getAllSameSizeNeighboursTwoLayers(const IPosition& ocentre, std::vector<IPosition>& centres);

  /** Get same or bigger size neighbours. */
  void getNeighbours(const IPosition& ocentre, const OCell<OCELL_DATA> &cell, 
    std::array<IPosition,6> &centres, std::array<OCell<OCELL_DATA> *,6> &cells);

  /** Has not only same-size, it has bigger neighbours. */
  bool hasBiggerNeighbours(const IPosition& ocentre, const OCell<OCELL_DATA> &cell);

  /** Get centres of smaller neighbours */
  void getFaceChildNeiCentres(const IPosition& ocentre, const int level, 
    std::array<std::array<IPosition,4>,6> &centres);

  //===== Find... ==============================================================

  /** Find a cell centre for a position of real coordinates. level must be from 0 
    (background) to maxLevel - 1 or greater (max level found). This cell may not exist at all -
    check by cells_.find(). */
  bool findCellCentre(const OVECTOR& position, const int level, IPosition& centre);

  /** Find out inside which cell a point is located and check if this cell really exists in octree. 
    Follow cell hierarchy till max maximum level specified; if level is maxLevel, a cell with max 
    possible level (min size) is returned. */
  bool findCellCentreAndCheck(const OVECTOR& position, const int level, IPosition& centre);

  /** Find back cell number */
  size_t findBackCell(const IPosition& ocentre);


  //===== Modifications ========================================================

  /** Add new cell to the list. Returns false if such node exists. It sets
    zero data to all added nodes. */
  bool addCell(const IPosition& coord, const OCell<OCELL_DATA>& cell);

  /** Add new node if it does not exist, otherwise node ref count incremented. */
  bool addNode(const IPosition &coord, const ONODE_DATA& data);

  /** Delete node */
  bool deleteNode(const IPosition &coord);

  /** Delete cell */
  bool deleteCell(const IPosition& ocentre, const OCell<OCELL_DATA> &cell);

  /** Refine cell. No interpolation for new nodes. */
  bool refineCell(const IPosition& ocentre, const OCell<OCELL_DATA> &cell);

  /** Refine cell, find it by coordinate. No interpolation for new nodes. */
  bool refineCell(const IPosition& coord);
  
  /** Refine list of cells. */
  void refineCells(const OCELLS &rastercells);

  /** Refine cells point cloud; refine octree intil cell conrin no more than numPointsPerCell. 
    Result is pointsInCells - a distribution of points over cells. */
  void refineCellsForPointCloud(const std::vector<OVECTOR> &points, const int numPointsPerCell,
    OPOINTS &pointsInCells);

  /** De-refine cell. This cell must contain 8 children every one of which
    cannot have neighbours of higher level. */
  bool derefineCell(const IPosition& ocentre, const OCell<OCELL_DATA>& cell, std::array<IPosition,8> &deleted);

  /** Derefine cell, find it by coordinate. */
  bool derefineCell(const IPosition& coord, std::array<IPosition,8> &deleted);

  /** Make list of rester cells (leaves only) and assign value 3 for all cell data. */
  void makeRasterCells(OPOINTS &pointsInCells, OCELLS &rastercells);

  /** Extend list of points in every cell with it's neighbours points. */
  void extendPointsInCells(OPOINTS &pointsInCells, OPOINTS &pointsInCellsExtended);

  /** Classify cells : data = 3 : raster, 1 : outside, 2 : inside; 0 - not leaves */
  void classifyCells(const OCELLS &rastercells);

  /** Undefined distance-field vector */
  const OVECTOR UndefinedVF = OVECTOR(
    std::numeric_limits<OREAL>::max(),
    std::numeric_limits<OREAL>::max(),
    std::numeric_limits<OREAL>::max());

  /** Add distance filed vectors to nodes. */
  void addDistanceFields(const std::vector<OVECTOR> &points, const OPOINTS &pointsInCells);

  /** Calculate signed distnace at nodes. Ssigned distance is stored in node data.second.W value. */
  void calcSignedDistance();

  /** Polyginise surface from distance field vectors in nodes. */
  void polygoniseDistanceFields(TTriangles<OREAL> &tris, OCELLS &rastercells);

  /** Mark all neighbours by mark. */
  int markAllNeighbours(OCELLS &cells, const IPosition &startcell, const char mark);

  /** Add same size neighbours around cells (to make refinement after rasterisation easier). */
  void addSameSizeNeighbours(OCELLS& cells);

  /** Convert cells for triangles for output */
  bool cellsToTriangles(const OCELLS &cells, TTriangles<OREAL>& tris);

  /** Convert cells with specified data value for triangles for output */
  bool cellsToTriangles(const OCELLS &cells, const char value, TTriangles<OREAL>& tris);


  /** Local node numbers to make triangles from faces */
  static const std::array<std::array<int,3>,12> faceTris;

  //===== Conformal finite-element approximation ===============================

  /** Get 27 integer coordinates for hex nodes for conforming FE approximation.
    Including central node. */
  void getConformingCoordinates(const IPosition& ocentre, IPosition coord[27]);

  /** Get 27 real coordinates for hex nodes for conforming FE approximation.
    Including central node. */
  void getRealConformingCoordinates(const IPosition& ocentre, OVECTOR coord[27]);

private:

  /** List of all cells. */
  OCELLS cells_;

  /** List of all nodes. */
  ONODES nodes_;

  /** To limit recursive depth in calling markAllNeighbours(). */
  const static int maxRecursiveDepth_ = 200;
  int recursiveDepth_ = 0;

  /** Empty all arrays. */
  void clearAll();

  /** Create background cells. 8 mln cells in 7 sec. in Release. */
  bool createBackground();

 };
