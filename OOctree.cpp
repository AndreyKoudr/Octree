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

                              // class header
#include "OOctree.h"
#include "MarchingCubes.h"
#include "Plane.h"
#include "Conforming3D.h"

#include <algorithm>
#include <limits>
#include <vector>
#include <set>
#include <assert.h>
#include <iostream>

using namespace std;

// Debug search routines
//#define DEBUG_FIND
#ifdef NDEBUG
  #undef DEBUG_FIND
#endif


bool IPosComp(const IPosition& lhs, const IPosition& rhs)
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


OOctree::OOctree(const OVECTOR &min, const OVECTOR &max, 
  const LINT Imax, const LINT Jmax, const LINT Kmax, const int pmaxLevel) : 
  OBackground(min,max,Imax,Jmax,Kmax,pmaxLevel)
{
  createBackground();
}

void OOctree::clearAll()
{
  cells_.clear();
  nodes_.clear();
}

size_t OOctree::numCells()
{
  return cells_.size();
}

size_t OOctree::numNodes()
{
  return nodes_.size();
}

OCELLS& OOctree::cells()
{
  return cells_;
}

ONODES& OOctree::nodes()
{
  return nodes_;
}


LINT OOctree::intCellSize(const int level)
{
  assert(level >= 0 && level < maxLevel);

  LINT shift = maxLevel - level - 1;
  LINT is = 2 << shift;
  return is;
}

OVECTOR OOctree::realCellSizes(const int level)
{
  assert(level >= 0 && level < maxLevel);

  LINT shift = maxLevel - level - 1;
  LINT is = 2 << shift;
  OVECTOR res = intUnitSizes * static_cast<OREAL>(is);
  return res;
}

bool OOctree::createBackground()
{
  // clear lists
  clearAll();

  // cell size for level 0
  LINT size = intCellSize(0);
  assert(size % 2 == 0);

#ifdef _DEBUG
  cout << "Creating background cells " << IJKnumCells.XYZ[0] << " x " << IJKnumCells.XYZ[1] <<
    " x " << IJKnumCells.XYZ[2] << endl;
#endif

  LINT x = size / 2;
  for (LINT i = 0; i < IJKnumCells.XYZ[0]; i++)
  {
#ifdef _DEBUG
    cout << i << " of " << IJKnumCells.XYZ[0] << "\r";
#endif
    LINT y = size / 2;
    for (LINT j = 0; j < IJKnumCells.XYZ[1]; j++)
    {
      LINT z = size / 2;
      for (LINT k = 0; k < IJKnumCells.XYZ[2]; k++)
      {
        OCell<OCELL_DATA> o(0,0,OCELL_DATA(0,-1));
        addCell(IPosition(x,y,z),o);

        z += size;
      }
      y += size;
    }
    x += size;
  }

#ifdef _DEBUG
  cout << endl;
#endif

  return true;
}

OVECTOR OOctree::intToRealCoord(const IPosition &coord)
{
  OVECTOR res;

  res.XYZ[0] = boxMin.XYZ[0] + coord.XYZ[0] * intUnitSizes.XYZ[0];
  res.XYZ[1] = boxMin.XYZ[1] + coord.XYZ[1] * intUnitSizes.XYZ[1];
  res.XYZ[2] = boxMin.XYZ[2] + coord.XYZ[2] * intUnitSizes.XYZ[2];

#ifndef NDEBUG
  // not quite right
  OREAL tol = std::numeric_limits<OREAL>::epsilon() * static_cast<OREAL>(10.0);

  assert(res.XYZ[0] >= boxMin.XYZ[0] - tol && res.XYZ[0] <= boxMax.XYZ[0] + tol);
  assert(res.XYZ[1] >= boxMin.XYZ[1] - tol && res.XYZ[1] <= boxMax.XYZ[1] + tol);
  assert(res.XYZ[2] >= boxMin.XYZ[2] - tol && res.XYZ[2] <= boxMax.XYZ[2] + tol);
#endif

  return res;
}

IPosition OOctree::realToIntCoord(const OVECTOR &coord)
{
  assert(coord.XYZ[0] >= boxMin.XYZ[0] && coord.XYZ[0] <= boxMax.XYZ[0]);
  assert(coord.XYZ[1] >= boxMin.XYZ[1] && coord.XYZ[1] <= boxMax.XYZ[1]);
  assert(coord.XYZ[2] >= boxMin.XYZ[2] && coord.XYZ[2] <= boxMax.XYZ[2]);

//!!!!!  OREAL tolerance = std::numeric_limits<OREAL>::epsilon() * static_cast<OREAL>(10.0);

  IPosition res;

  res.XYZ[0] = static_cast<LINT>((coord.XYZ[0] - boxMin.XYZ[0]) / intUnitSizes.XYZ[0]);
  res.XYZ[1] = static_cast<LINT>((coord.XYZ[1] - boxMin.XYZ[1]) / intUnitSizes.XYZ[1]);
  res.XYZ[2] = static_cast<LINT>((coord.XYZ[2] - boxMin.XYZ[2]) / intUnitSizes.XYZ[2]);

//!!!!!!
  //res.XYZ[0] = static_cast<LINT>((coord.XYZ[0] + tolerance - boxMin.XYZ[0]) / intUnitSizes.XYZ[0]);
  //res.XYZ[1] = static_cast<LINT>((coord.XYZ[1] + tolerance - boxMin.XYZ[1]) / intUnitSizes.XYZ[1]);
  //res.XYZ[2] = static_cast<LINT>((coord.XYZ[2] + tolerance - boxMin.XYZ[2]) / intUnitSizes.XYZ[2]);

  assert(res.XYZ[0] >= 0 && res.XYZ[0] <= IJKintSizes.XYZ[0]);
  assert(res.XYZ[1] >= 0 && res.XYZ[1] <= IJKintSizes.XYZ[1]);
  assert(res.XYZ[2] >= 0 && res.XYZ[2] <= IJKintSizes.XYZ[2]);

  return res;
}

IPosition OOctree::realToIntCoord(const OVECTOR &coord, const int level)
{
  assert(coord.XYZ[0] >= boxMin.XYZ[0] && coord.XYZ[0] <= boxMax.XYZ[0]);
  assert(coord.XYZ[1] >= boxMin.XYZ[1] && coord.XYZ[1] <= boxMax.XYZ[1]);
  assert(coord.XYZ[2] >= boxMin.XYZ[2] && coord.XYZ[2] <= boxMax.XYZ[2]);

  IPosition res;

  OVECTOR rSize = realCellSizes(level);
  LINT iSize = intCellSize(level);

  for (int j = 0; j < 3; j++)
  {
    res.XYZ[j] = static_cast<LINT>(ROUND((coord.XYZ[j] - boxMin.XYZ[j]) / rSize.XYZ[j])) * iSize;
  }

  assert(res.XYZ[0] >= 0 && res.XYZ[0] <= IJKintSizes.XYZ[0]);
  assert(res.XYZ[1] >= 0 && res.XYZ[1] <= IJKintSizes.XYZ[1]);
  assert(res.XYZ[2] >= 0 && res.XYZ[2] <= IJKintSizes.XYZ[2]);

  return res;
}

int OOctree::level(const IPosition &position)
{
  unsigned long index = 0;
  bool ok = _BitScanForward64(&index,position.XYZ[0]);
  assert(ok);

  if (!ok)
    return -1;

  int l = maxLevel - index - 1;
  assert(l <= maxLevel); //??? maybe <

  return l;
}

/**
  Cell node numeration :
                           7-------18--------6 
                          /|                /|    
                        19 |     25       17 |     
                        /  |              /  |          
                       4--------16-------5   14   
       K               |  15      23     |   |     
                       |   |             | 22|    
       |    J          |24 |    21       |   |
       |              12   3------10----13---2 
       |  /            |  /              |  / 
       | /             | 11      20      | 9
       |/              |/                |/
       *------ I       0--------8--------1
                                          
*/

/** Offsets from the centre. */
const std::array<IPosition,27> cellCentreOffsets = {
  IPosition(-1,-1,-1),
  IPosition(+1,-1,-1),
  IPosition(+1,+1,-1),
  IPosition(-1,+1,-1),
  IPosition(-1,-1,+1),
  IPosition(+1,-1,+1),
  IPosition(+1,+1,+1),
  IPosition(-1,+1,+1),

  IPosition( 0,-1,-1),
  IPosition(+1, 0,-1),
  IPosition( 0,+1,-1),
  IPosition(-1, 0,-1),

  IPosition(-1,-1, 0),
  IPosition(+1,-1, 0),
  IPosition(+1,+1, 0),
  IPosition(-1,+1, 0),

  IPosition( 0,-1,+1),
  IPosition(+1, 0,+1),
  IPosition( 0,+1,+1),
  IPosition(-1, 0,+1),

  IPosition( 0, 0,-1),
  IPosition( 0,-1, 0),
  IPosition(+1, 0, 0),
  IPosition( 0,+1, 0),
  IPosition(-1, 0, 0),
  IPosition( 0, 0,+1),

  IPosition( 0, 0, 0)
};

std::array<IPosition,8> OOctree::cell8Coordinates(const IPosition& ocentre)
{
  std::array<IPosition,8> coord;

  LINT size = intCellSize(level(ocentre));
  assert(size > 1);
  assert(size % 2 == 0);

  assert(ocentre.XYZ[0] % 2 == 0);
  assert(ocentre.XYZ[1] % 2 == 0);
  assert(ocentre.XYZ[2] % 2 == 0);

  LINT halfSize = size / 2;

  std::transform(cellCentreOffsets.begin(),cellCentreOffsets.begin() + 8,coord.begin(),
    [&ocentre,&halfSize](auto v) { return ocentre + v * halfSize; });

  return coord;
}

std::array<OVECTOR,8> OOctree::cell8RealCoordinates(const IPosition& ocentre)
{
  std::array<OVECTOR,8> coord;

  std::array<IPosition,8> icoord = cell8Coordinates(ocentre);

  std::transform(icoord.begin(),icoord.end(),coord.begin(),[&](auto v) { return intToRealCoord(v); } );

  return coord;
}

/**
  Child nodes :
                           24------25--------26 
                          /|                /|    
                        21 |     22       23 |     
                        /  |              /  |          
                      18--------19-------20  17   
       K(Z)            |  /|             |  /|     
                       | 12|             | 14|    
       |    J(Y)       |/  |             |/  |
       |               9   6-------7---- 11--8 
       |  /            |  /              |  / 
       | /             | 3       4       | 5
       |/              |/                |/
       *------ I(X)    0--------1--------2 
                                          
*/

/** Offsets from the centre. */
const std::array<IPosition,27> nodeChildDirs = {
  IPosition(-1,-1,-1),
  IPosition( 0,-1,-1),
  IPosition(+1,-1,-1),

  IPosition(-1, 0,-1),
  IPosition( 0, 0,-1),
  IPosition(+1, 0,-1),

  IPosition(-1,+1,-1),
  IPosition( 0,+1,-1),
  IPosition(+1,+1,-1),

  IPosition(-1,-1, 0),
  IPosition( 0,-1, 0),
  IPosition(+1,-1, 0),

  IPosition(-1, 0, 0),
  IPosition( 0, 0, 0),
  IPosition(+1, 0, 0),

  IPosition(-1,+1, 0),
  IPosition( 0,+1, 0),
  IPosition(+1,+1, 0),

  IPosition(-1,-1,+1),
  IPosition( 0,-1,+1),
  IPosition(+1,-1,+1),

  IPosition(-1, 0,+1),
  IPosition( 0, 0,+1),
  IPosition(+1, 0,+1),

  IPosition(-1,+1,+1),
  IPosition( 0,+1,+1),
  IPosition(+1,+1,+1)
};

std::array<IPosition,27> OOctree::childrenNodeCoordinates(const IPosition& ocentre)
{
  std::array<IPosition,27> coord;

  LINT size = intCellSize(level(ocentre));
  assert(size > 1);
  assert(size % 2 == 0);

  assert(ocentre.XYZ[0] % 2 == 0);
  assert(ocentre.XYZ[1] % 2 == 0);
  assert(ocentre.XYZ[2] % 2 == 0);

  LINT halfSize = size / 2;

  std::transform(nodeChildDirs.begin(),nodeChildDirs.end(),coord.begin(),
    [&ocentre,&halfSize](auto v) { return ocentre + v * halfSize; });

  return coord;
}

/**
  Sub cell numeration :
            
                           *-----------------* 
                          /   3    /   7    /|    
                         /-----------------/ |     
                        /   1    /   5    /|7|          
                       *-----------------* | |    
       K               |        |        |5 /|     
                       |   1    |   5    | / |    
       |    J          |        |        |/|6|
       |               |--------|--------| | * 
       |  /            |        |        |4|/ 
       | /             |   0    |   4    | /
       |/              |        |        |/
       *------ I       *-----------------*
                                          
*/

// Vectors pointing to parent centre depending on sub-cell index. 
const std::array<IPosition,8> subParentDirs = {
  IPosition(+1,+1,+1),
  IPosition(+1,+1,-1),
  IPosition(+1,-1,+1),
  IPosition(+1,-1,-1),
  IPosition(-1,+1,+1),
  IPosition(-1,+1,-1),
  IPosition(-1,-1,+1),
  IPosition(-1,-1,-1)
};


IPosition OOctree::parentCentre(const IPosition& ocentre, const OCell<OCELL_DATA> &cell)
{
  IPosition centre;

  assert(cell.level > 0);
  
  // get coordinates of parent centre
  LINT size = intCellSize(cell.level) / 2;
  centre = ocentre + subParentDirs[cell.subIndex] * size;

  return centre;
}

std::array<IPosition,8> OOctree::childrenCentres(const IPosition& ocentre)
{
  std::array<IPosition,8> centres;

  LINT csize = intCellSize(level(ocentre)) / 4;

  // "-" is correct - see subParentDirs[]
  std::transform(subParentDirs.begin(),subParentDirs.end(),
    centres.begin(),[&ocentre,&csize](IPosition v){ return ocentre - v * csize; });

  return centres;
}

IPosition OOctree::childCentre(const IPosition& ocentre, const int subIndex)
{
  assert(subIndex >= 0 && subIndex <= 7);
//!!!!!!!!!!!!!!!!!!!!!
  if (!(subIndex >= 0 && subIndex <= 7))
  {
    cout << "ERROR" << endl;
    int gsgsgs = 0;
  }



  LINT size = intCellSize(level(ocentre)) / 4;
  IPosition centre = ocentre - subParentDirs[subIndex] * size;
  return centre;
}

bool OOctree::isLeaf(const IPosition& ocentre)
{
  if (level(ocentre) >= maxLevel - 1)
    return true;

  // take a single child
  IPosition cCentre = childCentre(ocentre,0);

  OCELLS::iterator i = cells_.find(cCentre);

  // this is leaf if not found
  return (i == cells_.end());
}

size_t OOctree::centresToCells(const IPosition centres[], const size_t numCentres, 
  std::vector<OCell<OCELL_DATA> *> &cells)
{
  for (size_t i = 0; i < numCentres; i++)
  {
    OCELLS::iterator iter = cells_.find(centres[i]);
    if (iter != cells_.end())
    {
      cells.push_back(const_cast<OCell<OCELL_DATA> *>(&(iter->second)));
    }
  }

  return cells.size();
}

bool OOctree::getChildren(const IPosition& ocentre, std::array<IPosition,8> &centres, 
  std::vector<OCell<OCELL_DATA> *> &cells)
{
  centres = childrenCentres(ocentre);
  size_t count = centresToCells(&centres[0],8,cells);
  assert(count == 0 || count == 8);
  
  return (count > 0);
}

bool OOctree::addCell(const IPosition& ocentre, const OCell<OCELL_DATA> &cell)
{
  OCELLS::iterator it = cells_.find(ocentre);
  if (it == cells_.end())
  {
    // insert cell
    cells_.insert(std::pair<IPosition,OCell<OCELL_DATA> >(ocentre,cell));

    // insert 8 nodes
    std::array<IPosition,8> coord = cell8Coordinates(ocentre);
    for (size_t i = 0; i < 8; i++)
    {
      // does not set a node value
      addNode(coord[i],ONODE_DATA());
    }

    return true;
  } else
  {
    return false;
  }
}

bool OOctree::addNode(const IPosition &coord, const ONODE_DATA &data)
{
  ONODES::iterator it = nodes_.find(coord);

  if (it == nodes_.end())
  {
    ONode<ONODE_DATA> n(data);
    n.refCount = 1;

    nodes_.insert(std::pair<IPosition,ONode<ONODE_DATA> >(coord,n));

    return true;
  } else
  {
    it->second.refCount++;
    return false;
  }
}

std::pair<IPosition,IPosition> OOctree::minMaxCoordinates(const IPosition& ocentre)
{
  int l = level(ocentre);

  LINT size = intCellSize(l);
  assert(size > 1);
  assert(size % 2 == 0);

  assert(ocentre.XYZ[0] % 2 == 0);
  assert(ocentre.XYZ[1] % 2 == 0);
  assert(ocentre.XYZ[2] % 2 == 0);

  LINT halfSize = size / 2;

  IPosition min = ocentre - IPosition(halfSize,halfSize,halfSize);
  IPosition max = ocentre + IPosition(halfSize,halfSize,halfSize);

  return std::pair<IPosition,IPosition>(min,max);
}

std::pair<OVECTOR,OVECTOR> OOctree::minMaxCoordinates(const IPosition& ocentre, const OREAL tolerance)
{
  OVECTOR centre = intToRealCoord(ocentre);

  OVECTOR sizes = realCellSizes(level(ocentre));
  OVECTOR halfSizes = sizes * OREAL(0.5) + OVECTOR(tolerance,tolerance,tolerance);

  OVECTOR min = centre - halfSizes;
  OVECTOR max = centre + halfSizes;

  return std::pair<OVECTOR,OVECTOR>(min,max);
}

std::pair<OVECTOR,OVECTOR> OOctree::minMaxCoordinatesExtended(const IPosition& ocentre, const OREAL extension)
{
  assert(extension >= OREAL(1.0));

  OVECTOR centre = intToRealCoord(ocentre);

  OVECTOR sizes = realCellSizes(level(ocentre));
  OVECTOR halfSizes = sizes * (OREAL(0.5) * extension);

  OVECTOR min = centre - halfSizes;
  OVECTOR max = centre + halfSizes;

  return std::pair<OVECTOR,OVECTOR>(min,max);
}

/**
  Cell face numeration :
            
                           *-----------------*
                          /|                /|    
                         / |      5        / |     
                        /  |              /  |          
                       *------------3----*   |    
       K(Z)            |   |             |   |     
                       | 0 |             | 1 |    
       |    J(Y)       |   |             |   |
       |               |   *-----2-------|---*
       |  /            |  /              |  / 
       | /             | /        4      | /
       |/              |/                |/
       *------ I(X)    *-----------------*
*/

/** Possible same-size nighbour face directions for an cell. Neighbours are 
  numbered as faces. */
const std::array<IPosition,6> faceNeiDirs = {
  IPosition(-1, 0, 0),
  IPosition(+1, 0, 0),
  IPosition( 0,-1, 0),
  IPosition( 0,+1, 0),
  IPosition( 0, 0,-1),
  IPosition( 0, 0,+1)
};

std::array<IPosition,6> OOctree::faceNeiCentres(const IPosition& ocentre)
{
  std::array<IPosition,6> centres;

  int l = level(ocentre);
  LINT size = intCellSize(l);

  std::transform(faceNeiDirs.begin(),faceNeiDirs.end(),centres.begin(),
    [&ocentre,&size](auto v) { return ocentre + v * size; } );

  return centres;
}

bool OOctree::getSameSizeNeighbours(const IPosition& ocentre, std::array<IPosition,6> &centres, 
  std::array<OCell<OCELL_DATA> *,6> &cells)
{
  // check if such cell exists
  auto iter = cells_.find(ocentre);
  if (iter == cells_.end())
  {
    return false;
  }

  // get centres of neighbours
  centres = faceNeiCentres(ocentre);

  for (size_t i = 0; i < 6; i++)
  {
    assert(level(ocentre) == level(centres[i]));

    OCELLS::iterator iter = cells_.find(centres[i]);
    if (iter != cells_.end())
    {
      OCell<OCELL_DATA> *o = const_cast<OCell<OCELL_DATA> *>(&(iter->second));
      cells[i] = o;
    } else
    {
      // out of bounding box
      cells[i] = nullptr;
    }
  }

  return true;
}

void OOctree::getAllSameSizeNeighbours(const IPosition& ocentre, std::vector<IPosition>& centres)
{
  centres.clear();

  LINT size = intCellSize(level(ocentre));
  for (int i = -1; i <= +1; i++)
  {
    for (int j = -1; j <= +1; j++)
    {
      for (int k = -1; k <= +1; k++)
      {
        if (i == 0 && j == 0 && k == 0)
          continue;

        IPosition c = ocentre + IPosition(i,j,k) * size;

        OCELLS::iterator iter = cells_.find(c);
        if (iter != cells_.end())
        {
          centres.push_back(c);
        }
      }
    }
  }
}

void OOctree::getAllSameSizeNeighboursTwoLayers(const IPosition& ocentre, std::vector<IPosition>& centres)
{
  OPOSITIONS allcentres;

  std::vector<IPosition> centres0;
  getAllSameSizeNeighbours(ocentre,centres0);
  allcentres.insert(centres0.begin(),centres0.end());

  for (auto c : centres0)
  {
    std::vector<IPosition> centres1;
    getAllSameSizeNeighbours(c,centres1);
    allcentres.insert(centres1.begin(),centres1.end());
  }

  std::copy(allcentres.begin(),allcentres.end(),std::back_inserter(centres));
}

void OOctree::getNeighbours(const IPosition& ocentre, const OCell<OCELL_DATA> &cell, 
  std::array<IPosition,6> &centres, std::array<OCell<OCELL_DATA> *,6> &cells)
{
  LINT l = level(ocentre);
  if (l > 0)
  {
    std::array<IPosition,6> sameSizeCentres = faceNeiCentres(ocentre);
    IPosition parCentre = parentCentre(ocentre,cell);
    std::array<IPosition,6> parentSameSizeCentres = faceNeiCentres(parCentre);

    for (size_t i = 0; i < 6; i++)
    {
      OCELLS::iterator iter = cells_.find(sameSizeCentres[i]);
      if (iter != cells_.end())
      {
        OCell<OCELL_DATA> *o = const_cast<OCell<OCELL_DATA> *>(&(iter->second));
        centres[i] = sameSizeCentres[i];
        cells[i] = o;
      } else
      {
        OCELLS::iterator iter = cells_.find(parentSameSizeCentres[i]);
        if (iter != cells_.end())
        {
          OCell<OCELL_DATA> *o = const_cast<OCell<OCELL_DATA> *>(&(iter->second));
          centres[i] = parentSameSizeCentres[i];
          cells[i] = o;
        } else
        {
          // out of bounding box
          cells[i] = nullptr;
        }
      }
    }
  } else
  {
    // this is background cell of level 0, it cannot have parent
    getSameSizeNeighbours(ocentre,centres,cells);
  }
}


// Possible child-size nighbour directions for an cell in 1/4s of cell size. Neighbours are numbered as faces. 
const std::array<std::array<IPosition,4>,6> faceChildNeiDirs = {{
  // face 0
  {
    IPosition(-3,-1,-1),
    IPosition(-3,+1,-1),
    IPosition(-3,+1,+1),
    IPosition(-3,-1,+1)
  },
  {
    IPosition(+3,-1,-1),
    IPosition(+3,+1,-1),
    IPosition(+3,+1,+1),
    IPosition(+3,-1,+1)
  },
  {
    IPosition(-1,-3,-1),
    IPosition(+1,-3,-1),
    IPosition(+1,-3,+1),
    IPosition(-1,-3,+1)
  },
  {
    IPosition(-1,+3,-1),
    IPosition(+1,+3,-1),
    IPosition(+1,+3,+1),
    IPosition(-1,+3,+1)
  },
  {
    IPosition(-1,-1,-3),
    IPosition(+1,-1,-3),
    IPosition(+1,+1,-3),
    IPosition(-1,+1,-3)
  },
  {
    IPosition(-1,-1,+3),
    IPosition(+1,-1,+3),
    IPosition(+1,+1,+3),
    IPosition(-1,+1,+3)
  }
}};

void OOctree::getFaceChildNeiCentres(const IPosition& ocentre, const int level, 
  std::array<std::array<IPosition,4>,6> &centres)
{
  LINT size = intCellSize(level) / 4;
  for (size_t i = 0; i < 6; i++)
  {
    for (size_t j = 0; j < 4; j++)
    {
      centres[i][j] = ocentre + faceChildNeiDirs[i][j] * size;
    }
  }
}

bool OOctree::getSmallerNeighbours(const IPosition& ocentre, 
  std::array<std::array<IPosition,4>,6> &centres, std::array<std::array<OCell<OCELL_DATA> *,4>,6> &cells)
{
  getFaceChildNeiCentres(ocentre,level(ocentre),centres);

  bool found = false;
  for (size_t i = 0; i < 6; i++)
  {
    for (size_t j = 0; j < 4; j++)
    {
      auto iter = cells_.find(centres[i][j]);
      if (iter != cells_.end())
      {
        OCell<OCELL_DATA> *o = const_cast<OCell<OCELL_DATA> *>(&(iter->second));
        cells[i][j] = o;
        found = true;
      } else
      {
        cells[i][j] = nullptr;
      }
    }
  }

  return found;
}

bool OOctree::hasBiggerNeighbours(const IPosition& ocentre, const OCell<OCELL_DATA> &cell)
{
#if 0
  OCELLS::iterator i = cells_.find(ocentre);
  assert(i != cells_.end());

  IPosition parent = parentCentre(ocentre,i->second);

  std::array<IPosition,6> centres;
  std::array<OCell<OCELL_DATA> *,6> neighbours;
  getSameSizeNeighbours(parent,centres,neighbours);

  for (size_t i = 0; i < 6; i++)
  {
    if (neighbours[i] != nullptr)
    {
      if (isLeaf(centres[i]))
      {
        return true;
      }
    }
  }

  return false;

#else
  std::array<IPosition,6> centres;
  std::array<OCell<OCELL_DATA> *,6> neighbours;
  getNeighbours(ocentre,cell,centres,neighbours);

  bool has = false;

  int ol = level(ocentre);
  assert(ol == cell.level); 

  for (size_t i = 0; i < 6; i++)
  {
    if (neighbours[i] != nullptr)
    {
      if (isLeaf(centres[i]) && level(centres[i]) < cell.level)
      {
        has = true;
        break;
      }
    }
  }

  return has;
#endif
}

bool OOctree::refineCell(const IPosition& ocentre, const OCell<OCELL_DATA> &cell)
{
  if (
    // must be no children
    isLeaf(ocentre) &&
    // otherwise it is gonna be unbalanced 
    !hasBiggerNeighbours(ocentre,cell))
  {
    std::array<IPosition,8> centres = childrenCentres(ocentre);

    for (size_t i = 0; i < 8; i++)
    {
      // new cells
      OCell<OCELL_DATA> o(cell.level + 1,char(i),OCELL_DATA(0,cell.data.second));
//      OCell<OCELL_DATA> o(cell.level + 1,char(i),OCELL_DATA(0,-1));

      // create new nodes with interpolated values; unchanged if node
      // already exists
      addCell(centres[i],o);
    }

    return true;
  } else
  {
    return false;
  }
}

bool OOctree::refineCell(const IPosition &coord)
{
  OCELLS::iterator i = cells_.find(coord);
  if (i != cells_.end())
  {
    return refineCell(coord,i->second);
  } else
  {
    return false;
  }
}

bool OOctree::deleteNode(const IPosition &coord)
{
  ONODES::iterator it = nodes_.find(coord);
  if (it == nodes_.end())
  {
    return false;
  } else
  {
    it->second.refCount--;
    if (it->second.refCount <= 0)
    {
      // remove this node
      nodes_.erase(it);
      return true;
    } else
    {
      return false;
    }
  }
}

bool OOctree::deleteCell(const IPosition& ocentre, const OCell<OCELL_DATA> &cell)
{
  OCELLS::iterator it = cells_.find(ocentre);
  if (it != cells_.end())
  {
    // delete 8 nodes
    std::array<IPosition,8> coord = cell8Coordinates(ocentre);
    for (size_t i = 0; i < 8; i++)
    {
      deleteNode(coord[i]);
    }

    // delete cell
    cells_.erase(ocentre);

    return true;
  } else
  {
    return false;
  }
}

bool OOctree::derefineCell(const IPosition& ocentre, const OCell<OCELL_DATA> &cell, 
  std::array<IPosition,8> &deleted)
{
  // first, check if cells to be deleted do not have smaller neighbours
  // (to keep balancing)
  std::vector<OCell<OCELL_DATA> *> children;
  
  // no children, nothing to do
  if (!getChildren(ocentre,deleted,children))
    return false;

  assert(children.size() == 8);

  // cells to be deleted must not have smaller neighbours
  for (size_t i = 0; i < 8; i++)
  {
    std::array<IPosition,6> cs;
    std::array<OCell<OCELL_DATA> *,6> cells;
    getSameSizeNeighbours(deleted[i],cs,cells);

    for (size_t j = 0; j < 6; j++)
    {
      if (cells[j] != nullptr)
      { 
        if (!isLeaf(deleted[j]))
          return false;
      }
    }

    //if (hasSmallerNeighbours(deleted[i],*children[i]))
    //  return false;  
  }

  // ok, go on
  for (size_t i = 0; i < 8; i++)
  {
    deleteCell(deleted[i],*children[i]);
  }

  return true;
}

bool OOctree::derefineCell(const IPosition &coord, std::array<IPosition,8> &deleted)
{
  OCELLS::iterator i = cells_.find(coord);
  if (i != cells_.end())
  {
    return derefineCell(coord,i->second, deleted);
  } else
  {
    return false;
  }
}

size_t OOctree::findBackCell(const IPosition& ocentre)
{
  // cell size
  LINT size = intCellSize(0);

  // get background cell number
  IPosition backCell;

  for (size_t j = 0; j < 3; j++)
  {
    backCell.XYZ[j] = ocentre.XYZ[j] / size;
    assert(backCell.XYZ[j] >= 0 && backCell.XYZ[j] <= (IJKnumCells.XYZ[j] - 1));
    LIMIT(backCell.XYZ[j],0,IJKnumCells.XYZ[j] - 1);
  }

  // outside
  assert(backCell.XYZ[0] >= 0 && backCell.XYZ[0] < IJKnumCells.XYZ[0]);
  assert(backCell.XYZ[1] >= 0 && backCell.XYZ[1] < IJKnumCells.XYZ[1]);
  assert(backCell.XYZ[2] >= 0 && backCell.XYZ[2] < IJKnumCells.XYZ[2]);

  // get starting background cell
  size_t cellNo = positionToBackCellIndex(backCell);

  return cellNo;
}

bool OOctree::findCellCentre(const OVECTOR &position, const int plevel, IPosition &centre)
{

#ifdef DEBUG_FIND
  cout << "findCellCentre() for position " << position.XYZ[0] << " " << position.XYZ[1] << " " <<
    position.XYZ[2] << ", level " << plevel << endl;
#endif

  OREAL tolerance = std::numeric_limits<OREAL>::epsilon() * static_cast<OREAL>(10.0);

  // outside box
  for (size_t j = 0; j < 3; j++)
  {
    if (position.XYZ[j] < boxMin[j] - tolerance || position.XYZ[j] > boxMax[j] + tolerance)
      return false;
  }

  // get background cell number
  IPosition backCell;

  for (size_t j = 0; j < 3; j++)
  {
    backCell.XYZ[j] = static_cast<LINT>((position.XYZ[j] - boxMin[j]) / cellSize[j]);
    LIMIT(backCell.XYZ[j],0,IJKnumCells.XYZ[j] - 1);
  }

#ifdef DEBUG_FIND
  cout << "backCell " << backCell.XYZ[0] << " " << backCell.XYZ[1] << " " << backCell.XYZ[2] << endl;
#endif

  // outside
  if (backCell.XYZ[0] < 0 || backCell.XYZ[0] >= IJKnumCells.XYZ[0])
    return false;
  if (backCell.XYZ[1] < 0 || backCell.XYZ[1] >= IJKnumCells.XYZ[1])
    return false;
  if (backCell.XYZ[2] < 0 || backCell.XYZ[2] >= IJKnumCells.XYZ[2])
    return false;

  // get starting background cell
  size_t cellNo = positionToBackCellIndex(backCell);

  // supply some preliminary info on integer coordinates even if 
  // false is returned
  centre = IPosition(intCellSize(0) / 2,intCellSize(0) / 2,intCellSize(0) / 2) + backCell * intCellSize(0); 
  assert(cells_.find(centre) != cells_.end());

  // rough integer coordinates
  IPosition ipos = realToIntCoord(position);

  // look up sub-cell tree
  while (1)
  {
    int l = level(centre);
    if (l >= plevel)
      break; 

    // get sub cell
    LINT size = intCellSize(l) / 2;
    std::pair<IPosition,IPosition> minmax = minMaxCoordinates(centre);

    IPosition subNums = (ipos - minmax.first) / size;
    int index = static_cast<int>(subNums[0] * 4 + subNums[1] * 2 + subNums[2]);
    IPosition cCentre = childCentre(centre,index);

    centre = cCentre;
  }

  return true;
}

bool OOctree::findCellCentreAndCheck(const OVECTOR &position, const int plevel, IPosition &centre)
{
  bool ok = findCellCentre(position,plevel,centre);
  if (!ok)
    return false;

  auto iter = cells_.find(centre);
  return (iter != cells_.end());
}

static void lineBresenham(const LINT x0, const LINT y0, const LINT z0,
  const LINT x1, const LINT y1, const LINT z1, std::set<IPosition,IPosCompare> &points)
{
  //points.clear();

  LINT dx = x1 - x0;
  LINT dy = y1 - y0;
  LINT dz = z1 - z0;

  LINT ax = std::abs(dx) << 1;
  LINT ay = std::abs(dy) << 1;
  LINT az = std::abs(dz) << 1;

  LINT sx = (dx >= 0) ? +1 : -1;
  LINT sy = (dy >= 0) ? +1 : -1;
  LINT sz = (dz >= 0) ? +1 : -1;

  LINT x = x0;
  LINT y = y0;
  LINT z = z0;

  LINT xd,yd,zd;

  // along x
  if (ax >= std::max(ay,az))
  {
    yd = ay - (ax >> 1);
    zd = az - (ax >> 1);
    for (;;)
    {
      points.insert(IPosition(x,y,z));
      if (x == x1)
      {
        return;
      }

      if (yd >= 0)
      {
        y += sy;
        yd -= ax;
      }

      if (zd >= 0)
      {
        z += sz;
        zd -= ax;
      }

      x += sx;
      yd += ay;
      zd += az;
    }
  // along y
  } else if (ay >= std::max(ax,az)) 
  {
    xd = ax - (ay >> 1);
    zd = az - (ay >> 1);
    for (;;)
    {
      points.insert(IPosition(x,y,z));
      if (y == y1)
      {
        return;
      }

      if (xd >= 0)
      {
        x += sx;
        xd -= ay;
      }

      if (zd >= 0)
      {
        z += sz;
        zd -= ay;
      }

      y += sy;
      xd += ax;
      zd += az;
    }
  // along z
  } else if (az >= std::max(ax,ay))  
  {
    xd = ax - (az >> 1);
    yd = ay - (az >> 1);
    for (;;)
    {
      points.insert(IPosition(x,y,z));
      if (z == z1)
      {
        return;
      }

      if (xd >= 0)
      {
        x += sx;
        xd -= az;
      }

      if (yd >= 0)
      {
        y += sy;
        yd -= az;
      }

      z += sz;
      xd += ax;
      yd += ay;
    }
  }
}

void OOctree::refineCells(const OCELLS &rastercells)
{
  for (auto v : rastercells)
  {
    // pass tri number to this cell
    OCELLS::iterator i = cells_.find(v.first);
    if (i != cells_.end())
    {
      refineCell(i->first,i->second);
    }
  }
}

int OOctree::markAllNeighbours(OCELLS &cells, const IPosition &ocentre, const char mark)
{
  OCELLS::iterator it = cells_.find(ocentre); 
  assert(it != cells_.end());

  if (++recursiveDepth_ > maxRecursiveDepth_)
    return 0;

  // mark count
  int count = 0;

  // mark all same size and bigger neighbours
  {
    std::array<IPosition,6> ncentres;
    std::array<OCell<OCELL_DATA> *,6> ncells;
    getNeighbours(ocentre,it->second,ncentres,ncells);

    for (int i = 0; i < 6; i++)
    {
      if (ncells[i] != nullptr && isLeaf(ncentres[i]) && ncells[i]->data.first == 0)
      {
        ncells[i]->data.first = mark;
        count++;
        count += markAllNeighbours(cells,ncentres[i],mark);
      }
    }
  }

  // mark smaller neighbours
  {
    std::array<std::array<IPosition,4>,6> ncentres;
    std::array<std::array<OCell<OCELL_DATA> *,4>,6> ncells;
    getSmallerNeighbours(ocentre,ncentres,ncells);
    for (int i = 0; i < 6; i++)
    {
      for (int j = 0; j < 4; j++)
      {
        if (ncells[i][j] != nullptr && isLeaf(ncentres[i][j]) && ncells[i][j]->data.first == 0)
        {
          ncells[i][j]->data.first = mark;
          count++;
          count += markAllNeighbours(cells,ncentres[i][j],mark);
        }
      }
    }
  }

  return count;
}


/*
  Cell face numeration :
                          (7)---------------(6)
                          /|                /|    
                         / |      5        / |     
                        /  |              /  |          
                      (4)-----------3---(5)  |    
       K(Z)            |   |             |   |     
                       | 0 |             | 1 |    
       |    J(Y)       |   |             |   |
       |               |  (3)----2-------|--(2)
       |  /            |  /              |  / 
       | /             | /        4      | /
       |/              |/                |/
       *------ I(X)   (0)---------------(1)
                                          
*/

const std::array<std::array<int,3>,12> OOctree::faceTris = {{
  {3,0,4},
  {3,4,7},
  {1,2,6},
  {1,6,5},
  {0,1,5},
  {0,5,4},
  {2,3,7},
  {2,7,6},
  {3,2,1},
  {3,1,0},
  {4,5,6},
  {4,6,7}
}};

std::pair<OVECTOR,OVECTOR> OOctree::minmax(const OCELLS &cells, const char value)
{
  OVECTOR min = +std::numeric_limits<OREAL>::max();
  OVECTOR max = -std::numeric_limits<OREAL>::max();

  for (auto c : cells)
  {
    if (isLeaf(c.first))
    {
      if (c.second.data.first == value)
      {
        std::pair<OVECTOR,OVECTOR> mima = minMaxCoordinates(c.first,OREAL(0.0));
        for (size_t j = 0; j < 3; j++)
        {
          min.XYZ[j] = std::min(mima.first.XYZ[j],min.XYZ[j]);
          max.XYZ[j] = std::max(mima.second.XYZ[j],max.XYZ[j]);
        }
      }
    }
  }

  return std::pair<OVECTOR,OVECTOR>(min,max);
}

void OOctree::addSameSizeNeighbours(OCELLS& cells)
{
  OCELLS newCells;

  for (auto it = cells.begin(); it != cells.end(); it++)
  {
    std::array<IPosition,6> centres; 
    std::array<OCell<OCELL_DATA> *,6> cells;
    if (getSameSizeNeighbours(it->first,centres,cells))
    {
      for (size_t i = 0; i < 6; i++)
      {
        if (cells[i] != nullptr)
        {
          // copy data with tri number; level and subindex remain wrong
          newCells.insert(std::pair<IPosition,OCell<OCELL_DATA> >(centres[i],it->second));
        }
      }
    }
  }

  for (OCELLS::iterator it = newCells.begin(); it != newCells.end(); it++)
  {
    cells.insert(*it);
  }  
}

void OOctree::extractNodes(const OCELLS &cells, const char value, OPOSITIONS &nodes)
{
  for (auto c : cells)
  {
    if (isLeaf(c.first))
    {
      if (c.second.data.first == value)
      {
        std::array<IPosition,8> coords = cell8Coordinates(c.first);

        for (auto c : coords)
        {
          nodes.insert(c);
        }
      }
    }
  }
}

/**
  Face nodes :
            
                           7-----------------6
                          /|                /|    
                         / |     (5)       / |     
                        /  |              /  |          
                       4-----------(3)---5   |    
       K(Z)            |   |             |   |     
                       |(0)|             |(1)|    
       |    J(Y)       |   |             |   |
       |               |   3----(2)------|---2
       |  /            |  /              |  / 
       | /             | /       (4)     | /
       |/              |/                |/
       *------ I(X)    0-----------------1
*/

const std::array<std::array<LINT,4>,6> faceNodes = {{
  {3,0,4,7},
  {1,2,6,5},
  {0,1,5,4},
  {2,3,7,6},
  {1,0,3,2},
  {4,5,6,7}
}};

void OOctree::extractBoundaryFaces(const OCELLS &cells, const char value,  
  const char rastervalue, std::vector<std::array<IPosition,4>> &faces, 
  std::vector<IPosition> &facecells, std::vector<LINT> &tris, 
  std::map<IPosition,std::set<LINT>,IPosCompare>& postotri)
{
  for (auto c : cells)
  {
    if (isLeaf(c.first))
    {
      if (c.second.data.first == value)
      {
        std::array<IPosition,6> centres;
        std::array<OCell<OCELL_DATA> *,6> neighbours;

        // it is important : neighbours between raster (value 3) and inner cells (value 1)
        // must be same size, normally of highest level
        getSameSizeNeighbours(c.first,centres,neighbours);

        for (size_t j = 0; j < 6; j++)
        {
          if (neighbours[j] != nullptr && neighbours[j]->data.first == rastervalue)
          {
            assert(isLeaf(centres[j]));
            assert(level(c.first) == level(centres[j]));

            std::array<IPosition,8> nodes = cell8Coordinates(c.first);
            std::array<IPosition,4> face({nodes[faceNodes[j][0]],nodes[faceNodes[j][1]],
              nodes[faceNodes[j][2]],nodes[faceNodes[j][3]]});

            faces.push_back(face);
            facecells.push_back(c.first);
            tris.push_back(neighbours[j]->data.second);

            LINT tri = neighbours[j]->data.second;
            assert(tri >= 0);
            for (int k = 0; k < 4; k++)
            {
              auto it = postotri.find(nodes[faceNodes[j][k]]);
              if (it == postotri.end())
              {
                postotri.insert(std::pair<IPosition,std::set<LINT>>(nodes[faceNodes[j][k]],
                  std::set<LINT>(initializer_list<LINT>{ tri })));
              } else
              {
                it->second.insert(tri);
              }
            }
          }
        }
      }
    }
  }
}

void OOctree::refineCellsForPointCloud(const std::vector<OVECTOR> &points, const int numPointsPerCell,
  OPOINTS &pointsInCells)
{
  assert(numPointsPerCell > 0);

  // refine cells
  //for (int level = 0; level < 4; level++)
  for (int level = 0; level < maxLevel; level++)
  {
    // find ALL cells which contain at least one point
    for (size_t i = 0; i < points.size(); i++)
    {
      IPosition centre;
      if (findCellCentreAndCheck(points[i],level,centre))
      {
        assert(OOctree::level(centre) == level);
        pointsInCells[centre].insert(i);
      }
    }

    bool found = false;
    for (auto &list : pointsInCells)
    {
      if (list.second.size() > numPointsPerCell)
      {
        refineCell(list.first);

        // refine closest neighbours as well to allow more freedom for cell refinement
        std::vector<IPosition> centres;
//!!!
#if 1
        // refine only immediate neighbours
        getAllSameSizeNeighbours(list.first,centres);
#else
        // refine two layers of neighbours, a liitle bit better result
        getAllSameSizeNeighboursTwoLayers(list.first,centres);
#endif

        for (auto c : centres)
        {
          refineCell(c);
        }
        
        // keep points inside all cells, do not clear parent
        //list.second.clear();
        found = true;
      }
    }

    if (!found)
    {
      // mark ALL cells which contain a point as raster cells
      for (size_t i = 0; i < points.size(); i++)
      {
        IPosition centre;
        if (findCellCentreAndCheck(points[i],level + 1,centre))
        {
          pointsInCells[centre].insert(i);
        }
      }

      // stop
      break;
    }
  }


#if 0
  // too small cells make hole in the surface; we need to get rid of them
  int minLevel = 10000;
  int maxLevel = -10000;

  for (auto c : rastercells)
  {
    int l = level(c.first);
    minLevel = std::min<int>(minLevel,l);
    maxLevel = std::max<int>(maxLevel,l);
  }

//!!!

//  printf("min level %d max level %d",minLevel,maxLevel);

  if (maxLevel > minLevel)
  {
    for (int lev = maxLevel; lev > minLevel; lev--)
    {
      OPOSITIONS todelete;
      for (auto c : rastercells)
      {
        int l = level(c.first);
        if (l == lev)
        {
          std::array<IPosition,8> deleted;
          if (derefineCell(parentCentre(c.first,c.second),deleted))
          {
            for (auto d : deleted)
            {
              todelete.insert(d);
            }
          }
        }
      }

      for (auto d : todelete)
      {
        pointsInCells.erase(d);
        rastercells.erase(d);
      }
    }
  }
#endif

}

void OOctree::makeRasterCells(OPOINTS &pointsInCells, OCELLS &rastercells)
{
  // mark all cells as "unknown"
  for (auto &c : cells_)
  {
    c.second.data.first = 0;
  }

  // mark all leaf cells which contain at least one point as "raster"
  for (auto &pc : pointsInCells)
  {
    cells_[pc.first].data.first = 3;

    OCELLS::iterator it = cells_.find(pc.first);
    if (it != cells_.end())
    {
      if (isLeaf(it->first))
        rastercells.insert(std::pair<IPosition,OCell<OCELL_DATA> >(it->first,it->second));
    }
  }
}

void OOctree::extendPointsInCells(OPOINTS &pointsInCells, OPOINTS &pointsInCellsExtended)
{
  // extend list of points for each cell
  for (auto pc : pointsInCells)
  {
    IPosition ocentre = pc.first;

    // extend point list for cell ocentre by its neighbours

    std::vector<IPosition> neighbours;
    getAllSameSizeNeighbours(ocentre,neighbours);

    pointsInCellsExtended[ocentre].insert(pointsInCells[ocentre].begin(),pointsInCells[ocentre].end());

    for (auto nei : neighbours)
    {
      if (pointsInCells.find(nei) != pointsInCells.end() && pointsInCells[nei].size())
      {
        pointsInCellsExtended[ocentre].insert(pointsInCells[nei].begin(),pointsInCells[nei].end());
      }
    }

    // take all its neighbours (if they even have no points inside) and add all THEIR points
    // - this is done to avoid "holes" and "slits" between generated triangles

    for (auto nei : neighbours)
    {
      std::vector<IPosition> neineighbours;
      getAllSameSizeNeighbours(nei,neineighbours);

      // insert its own points if any
      if (pointsInCells.find(nei) != pointsInCells.end())
        pointsInCellsExtended[nei].insert(pointsInCells[nei].begin(),pointsInCells[nei].end());

      // take all neighbour' neighbours
      for (auto neinei : neineighbours)
      {
        if (pointsInCells.find(neinei) != pointsInCells.end() && pointsInCells[neinei].size())
        {
          pointsInCellsExtended[nei].insert(pointsInCells[neinei].begin(),pointsInCells[neinei].end());
        }
      }
    }
  }
}

void OOctree::classifyCells(const OCELLS &rastercells)
{
  int rasterlevel = level(rastercells.begin()->first);

  // raster cells are marked by "3" 
  for (auto &v : cells_)
  {
    auto r = rastercells.find(v.first);
    if (r != rastercells.end())
    {
      v.second.data.first = 3;
      v.second.data.second = r->second.data.second;
    } else
    {
      v.second.data.first = 0;
    }
  }

  // find first cell NOT marked as raster, it maybe both "outside (1)" or "inside (2)"
  for (auto &v : cells_)
  {
    if (isLeaf(v.first) && level(v.first) == rasterlevel && v.second.data.first == 0)
    {
      // mark this cell as first "not raster"
      v.second.data.first = 1;
      break;
    }
  }

  // find all neighbour cells (leafs) 
  bool found = false;
  do {
    found = false;
    for (auto &v : cells_)
    {
      if (isLeaf(v.first) && v.second.data.first == 1)
      {
        recursiveDepth_ = 0;
        if (markAllNeighbours(cells_,v.first,1) > 0)
          found = true;
      }
    }
  } while (found);

  // mark all remaining leaves as 2 (opposite to 1 - inside/outside)
  for (auto &v : cells_)
  {
    if (isLeaf(v.first) && v.second.data.first != 1 && v.second.data.first != 3)
    {
      // mark this cell as first "not raster"
      v.second.data.first = 2;
    }
  }

  // assign tri number to 1 and 2 cells; not used in point clouds
  for (auto &v : cells_)
  {
    if (isLeaf(v.first) && (v.second.data.first == 1 || v.second.data.first == 2))
    {
      std::array<IPosition,6> ncentres;
      std::array<OCell<OCELL_DATA> *,6> ncells;
      getSameSizeNeighbours(v.first,ncentres,ncells);

      for (int i = 0; i < 6; i++)
      {
        if (ncells[i] != nullptr && isLeaf(ncentres[i]) && ncells[i]->data.first == 3)
        {
          v.second.data.second = ncells[i]->data.second;
        }
      }
    }
  }
}

void OOctree::addDistanceFields(const std::vector<OVECTOR> &points, const OPOINTS &pointsInCells)
{
  // fill all nodes data with unknown value
  for (auto &node : nodes_)
  {
    node.second.data = UndefinedVF;
  }

#if 0
  for (auto &node : nodes_)
  {
    OVECTOR c = intToRealCoord(node.first);

    OREAL minDist = std::numeric_limits<OREAL>::max();
    for (auto p : points)
    {
      OVECTOR df = p - c;
      OREAL dist = !df;
      if (dist < minDist)
      {
        minDist = dist;
        node.second.data = df;
      }
    }
  }

int gsggsgs = 0;
#else
  for (auto cp : pointsInCells)
  {
    // this is a cell, get its corners (nodes)
    IPosition centre = cp.first;
    std::array<IPosition,8> icoord = cell8Coordinates(centre);

    // get real coordinates
    std::array<OVECTOR,8> coord;
    std::transform(icoord.begin(),icoord.end(),coord.begin(),[&](auto v) { return intToRealCoord(v); } );

#if 1
    // make plane
    if (cp.second.size() >= 3)
    {
      bool ok = true;
      std::vector<LINT> cpv(cp.second.begin(),cp.second.end()); 

      //cpv.resize(3);
      //std::vector<OVECTOR> vpoints = {points[cpv[0]],points[cpv[1]],points[cpv[2]]};
      std::vector<OVECTOR> vpoints;
      for (auto i : cpv)
      {
        vpoints.push_back(points[i]);
      }

      Plane<OREAL> plane(vpoints,ok);
      
      if (!ok)
      {
        plane = Plane<OREAL>(points[cpv[0]],points[cpv[1]],points[cpv[2]],ok);
      } else
      {
        //Plane<OREAL> plane1 = Plane<OREAL>(points[cpv[0]],points[cpv[1]],points[cpv[2]],ok);
        int gsgsgs = 0;
      }
      
      if (ok)
      {
        for (size_t i = 0; i < 8; i++)
        {
          // project node on plane
          ONode<ONODE_DATA> &node = nodes_[icoord[i]];
          OVECTOR p = plane.makeProjection(coord[i]);
          OVECTOR distFieldVec = p - coord[i];
          OREAL dist = !distFieldVec;

          // not yet inited
          if (node.data == UndefinedVF)
          {
            node.data = distFieldVec;
          } else
          {
            // assign node a shortest distance field vector
            OREAL oldDist = !node.data;
            if (dist < oldDist)
            {
              node.data = distFieldVec;
            }
          }
        }
      } else
      {
        for (size_t i = 0; i < 8; i++)
        {
          ONode<ONODE_DATA> &node = nodes_[icoord[i]];

          // check all cloud points in this cell
          for (auto ip : cp.second)
          {
            // cloud point
            OVECTOR p = points[ip];

            // distance field vector from node to cloud point
            OVECTOR distFieldVec = p - coord[i];
            OREAL dist = !distFieldVec;

            // not yet inited
            if (node.data == UndefinedVF)
            {
              node.data = distFieldVec;
            } else
            {
              //// assign node a shortest distance field vector
              //OREAL oldDist = !node.data;
              //if (dist < oldDist)
              //{
              //  node.data = distFieldVec;
              //}
            }
          }
        }
      }
    } else
    {
      // less than three points; just calc min distance to a point inside
      for (size_t i = 0; i < 8; i++)
      {
        ONode<ONODE_DATA> &node = nodes_[icoord[i]];

        // check all cloud points in this cell
        for (auto ip : cp.second)
        {
          // cloud point
          OVECTOR p = points[ip];

          // distance field vector from node to cloud point
          OVECTOR distFieldVec = p - coord[i];
          OREAL dist = !distFieldVec;

          // not yet inited
          if (node.data == UndefinedVF)
          {
            node.data = distFieldVec;
          } else
          {
            //// assign node a shortest distance field vector
            //OREAL oldDist = !node.data;
            //if (dist < oldDist)
            //{
            //  node.data = distFieldVec;
            //}
          }
        }
      }
    }
#else

    for (size_t i = 0; i < 8; i++)
    {
      ONode<ONODE_DATA> &node = nodes_[icoord[i]];

      // check all cloud points in this cell
      for (auto ip : cp.second)
      {
        // cloud point
        OVECTOR p = points[ip];

        // distance field vector from node to cloud point
        OVECTOR distFieldVec = p - coord[i];
        OREAL dist = !distFieldVec;

        // not yet inited
        if (node.data == UndefinedVF)
        {
          node.data = distFieldVec;
        } else
        {
          // assign node a shortest distance field vector
          OREAL oldDist = !node.data;
          if (dist < oldDist)
          {
            node.data = distFieldVec;
          }
        }
      }
    }
#endif
  }
#endif

#ifdef _DEBUG
  // count number of nodes with unassigned distance field vector
  {
    size_t count = 0;
    for (auto node : nodes_)
    {
      if (node.second.data == UndefinedVF)
        count++;
    }

    int gsgsgsg = 0;
  }
#endif
}

void OOctree::calcSignedDistance()
{
  // get list of outside nodes
  std::set<IPosition,IPosCompare> nodeOutside;

  for (auto cell : cells_)
  {
    IPosition centre = cell.first;

    // this is a cell, get its corners (nodes)
    std::array<IPosition,8> icoord = cell8Coordinates(centre);

    if (cell.second.data.first == 1)
    {
      for (auto i : icoord)
      {
        nodeOutside.insert(i);
      }
    }
  }

  // calc signed distance from distance field vector length and inside/outside info
  // and store it in node data.second.W
  for (auto &cell : cells_)
  {
    IPosition centre = cell.first;

    if (isLeaf(centre) && cell.second.data.first == 3) //???
    {
      // this is a cell, get its corners (nodes)
      std::array<IPosition,8> icoord = cell8Coordinates(centre);

      // all nodes must have distance field vector
      for (auto nodeCentre : icoord)
      {
        if (!(nodes_[nodeCentre].data == UndefinedVF))
        {
          if (nodeOutside.find(nodeCentre) != nodeOutside.end())
          {
            nodes_[nodeCentre].data.W = !nodes_[nodeCentre].data;
          } else
          {
            nodes_[nodeCentre].data.W = -(!nodes_[nodeCentre].data);
          }
        } else
        {
          // this is incorrect and should never happen
          nodes_[nodeCentre].data.W = 0.0;
        }
      }
    }
  }
}

/**
  Sub cell numeration :
            
                           *-----------------* 
                          /   3    /   7    /|    
                         /-----------------/ |     
                        /   1    /   5    /|7|          
                       *-----------------* | |    
       K               |        |        |5 /|     
                       |   1    |   5    | / |    
       |    J          |        |        |/|6|
       |               |--------|--------| | * 
       |  /            |        |        |4|/ 
       | /             |   0    |   4    | /
       |/              |        |        |/
       *------ I       *-----------------*
                                          
  Conforming node numeration :
            
                           7-------18--------6 
                          /|                /|    
                        19 |     25       17 |     
                        /  |              /  |          
                       4--------16-------5   14   
       K               |  15      23     |   |     
                       |   |             | 22|    
       |    J          |24 |    21       |   |
       |              12   3------10----13---2 
       |  /            |  /              |  / 
       | /             | 11      20      | 9
       |/              |/                |/
       *------ I       0--------8--------1
                                          
*/

/** Local coordinates of conforming nodes, [0..1], last 26-th node is not used but added to specify 
  its local coordinate in sub cells */
const OVECTOR conformingNodeLocalCoords[27] = {
  OVECTOR(0,0,0),
  OVECTOR(1,0,0),
  OVECTOR(1,1,0),
  OVECTOR(0,1,0),

  OVECTOR(0,0,1),
  OVECTOR(1,0,1),
  OVECTOR(1,1,1),
  OVECTOR(0,1,1),

  OVECTOR(0.5,0,0),
  OVECTOR(1,0.5,0),
  OVECTOR(0.5,1,0),
  OVECTOR(0,0.5,0),

  OVECTOR(0,0,0.5),
  OVECTOR(1,0,0.5),
  OVECTOR(1,1,0.5),
  OVECTOR(0,1,0.5),

  OVECTOR(0.5,0,1),
  OVECTOR(1,0.5,1),
  OVECTOR(0.5,1,1),
  OVECTOR(0,0.5,1),

  OVECTOR(0.5,0.5,0),
  OVECTOR(0.5,0,0.5),
  OVECTOR(1,0.5,0.5),
  OVECTOR(0.5,1,0.5),
  OVECTOR(0,0.5,0.5),
  OVECTOR(0.5,0.5,1),

  OVECTOR(0.5,0.5,0.5)
};

// 8 sub cells defined by local node numbers. 
const LINT octreeSubCells[8][8] =
{
  {0,8,20,11,12,21,26,24},
  {12,21,26,24,4,16,25,19},
  {11,20,10,3,24,26,23,15},
  {24,26,23,15,19,25,18,7},

  {8,1,9,20,21,13,22,26},
  {20,9,2,10,26,22,14,23},
  {21,13,22,26,16,5,17,25},
  {26,22,14,23,25,17,6,18}
};

void OOctree::polygoniseDistanceFields(TTriangles<OREAL> &tris, OCELLS &rastercells)
{
  // tolerance
  OREAL tol = scale * 0.00001;

  for (auto &cell : rastercells)
  {
    IPosition centre = cell.first;

    // this is a cell, get its corners (nodes)
    std::array<IPosition,8> icoord = cell8Coordinates(centre);

    // all nodes must have distance field vector
    bool ok = true;
    for (auto nodeCentre : icoord)
    {
      if (nodes_[nodeCentre].data == UndefinedVF)
      {
        ok = false;
        break;
      }
    }

    if (ok)
    {
      // apply conforming finite elements to provide continuous approximation for signed distance function
      IPosition confcoord[27];
      getConformingCoordinates(centre,confcoord);

      OVECTOR rconfcoord[27];
      for (int j = 0; j < 27; j++)
      {
        rconfcoord[j] = intToRealCoord(confcoord[j]);
      }

      std::array<LINT,26> confnodes;
      OREAL cellvalues[26] = {0};

      for (LINT i = 0; i < 26; i++)
      {
        auto it = nodes_.find(confcoord[i]);
#if 0
        if (it != nodes_.end())
#else
        if (it != nodes_.end() && !(nodes_[confcoord[i]].data == UndefinedVF))
#endif
        {
          confnodes[i] = i;
          cellvalues[i] = it->second.data.W;
        } else
        {
          confnodes[i] = -1;
          cellvalues[i] = 0.0;
        }
      }

#if 0
      // all nodes must have distance field vector
      bool ok1 = true;
      for (int i = 0; i < 26; i++)
      {
        if (confnodes[i] >= 0)
        {
          if (nodes_[confcoord[i]].data == UndefinedVF)
          {
            ok1 = false;
            break;
          }
        }
      }

      if (!ok1)
      {
        int result = Polygonise<OREAL>(&rconfcoord[0],cellvalues,0,tris,0,0,0,tol);
        continue;
      }  
#endif

      Conforming3D<OREAL,OVECTOR> element(confnodes);

      if (element.hasHangingNodes())
      {
#if 0
        int result = Polygonise<OREAL>(&rconfcoord[0],cellvalues,0,tris,0,0,0,tol);
#else
        // polygonise sub-cells
        for (size_t i = 0; i < 8; i++)
        {
          // local 8 coordinates for every sub cell
          OVECTOR localCoords[8];
          for (size_t j = 0; j < 8; j++)
          {
            localCoords[j] = conformingNodeLocalCoords[octreeSubCells[i][j]];
          }

          OVECTOR coords[8];
          OREAL values[8];

          for (size_t j = 0; j < 8; j++)
          {
            std::array<OREAL,26> shapefuncs;
            element.shapeFunctions(localCoords[j],shapefuncs);
            coords[j] = OVECTOR(0,0,0);
            for (size_t k = 0; k < 8; k++)
            {
              coords[j] += rconfcoord[k] * shapefuncs[k];
            }

            element.shapeFunctionsConforming(localCoords[j],shapefuncs);
            values[j] = OREAL(0.0);
            for (size_t k = 0; k < 26; k++)
            {
              values[j] += cellvalues[k] * shapefuncs[k];
            }
          }

          int result = Polygonise<OREAL>(&coords[0],values,0,tris,0,0,0,tol);
        }
#endif
      } else
      {
        int result = Polygonise<OREAL>(&rconfcoord[0],cellvalues,0,tris,0,0,0,tol);
      }
    }
  }
}

#if 0
void OOctree::polygoniseDistanceFields(TTriangles<OREAL> &tris, std::vector<OVECTOR> &points, 
  OPOINTS &pointsInCells, OCELLS &rastercells)
{
  // tolerance
  OREAL tol = scale * 0.00001;

TTriangles<OREAL> etris;

//  for (auto &cell : cells_)
  for (auto &cell : rastercells)
  {
    IPosition centre = cell.first;

 //   if (isLeaf(centre) && cell.second.data.first == 3)
    {
      // this is a cell, get its corners (nodes)
      std::array<IPosition,8> icoord = cell8Coordinates(centre);

      // all nodes must have distance field vector
      bool ok = true;
      for (auto nodeCentre : icoord)
      {
        if (nodes_[nodeCentre].data == UndefinedVF)
        {
//if (result == 0)
{
 // cell.second.data.first = 5;
}

          ok = false;
          break;
        }
      }

      if (ok)
      {
        // apply conforming finite elements to provide continuous approximation for signed distance function
        IPosition confcoord[27];
        getConformingCoordinates(centre,confcoord);

        OVECTOR rconfcoord[27];
        for (int j = 0; j < 27; j++)
        {
          rconfcoord[j] = intToRealCoord(confcoord[j]);
        }

        std::array<LINT,26> confnodes;
        OREAL cellvalues[26] = {0};

        for (LINT i = 0; i < 26; i++)
        {
          auto it = nodes_.find(confcoord[i]);
#if 0
          if (it != nodes_.end())
#else
          if (it != nodes_.end() && !(nodes_[confcoord[i]].data == UndefinedVF))
#endif
          {
            confnodes[i] = i;
            cellvalues[i] = it->second.data.W;
          } else
          {
            confnodes[i] = -1;
            cellvalues[i] = 0.0;
          }
        }

//!!!!
//for (int i = 0; i < 8; i++)
//{
//  printf("%f %f %f  ",rconfcoord[i].X,rconfcoord[i].Y,rconfcoord[i].Z);
//}
//printf("\n");

        int result = Polygonise<OREAL>(&rconfcoord[0],cellvalues,0,tris,0,0,0,tol);
        Polygonise<OREAL>(&rconfcoord[0],cellvalues,0,etris,0,0,0,tol);
if (result == 2)
{
  for (int i = 0; i < 8; i++)
  {
    printf("%f ",cellvalues[i]);
  }
  printf("result %d\n",result);
}

if (result == 0)
{
  cell.second.data.first = 5;
}
if (result == 1)
{
  cell.second.data.first = 6;
}
if (result == 2)
{
  cell.second.data.first = 7;
}
if (result == 3)
{
  cell.second.data.first = 8;
}

        continue;

#if 0
        // all nodes must have distance field vector
        bool ok1 = true;
        for (int i = 0; i < 26; i++)
        {
          if (confnodes[i] >= 0)
          {
            if (nodes_[confcoord[i]].data == UndefinedVF)
            {
              ok1 = false;
              break;
            }
          }
        }

        if (!ok1)
        {
          int result = Polygonise<OREAL>(&rconfcoord[0],cellvalues,0,tris,0,0,0,tol);
if (result == 0)
{
  cell.second.data.first = 6;
}
          continue;
        }  
#endif

        Conforming3D<OREAL,OVECTOR> element(confnodes);

        if (element.hasHangingNodes())
        {
#if 0
          int result = Polygonise<OREAL>(&rconfcoord[0],cellvalues,0,tris,0,0,0,tol);

#else
          for (size_t i = 0; i < 8; i++)
          {
            // local 8 coordinates for every sub cell
            OVECTOR localCoords[8];
            for (size_t j = 0; j < 8; j++)
            {
              localCoords[j] = conformingNodeLocalCoords[octreeSubCells[i][j]];
            }

            OVECTOR coords[8];
            OREAL values[8];

            for (size_t j = 0; j < 8; j++)
            {
              std::array<OREAL,26> shapefuncs;
              element.shapeFunctions(localCoords[j],shapefuncs);
              coords[j] = OVECTOR(0,0,0);
              for (size_t k = 0; k < 8; k++)
              {
                coords[j] += rconfcoord[k] * shapefuncs[k];
              }

              element.shapeFunctionsConforming(localCoords[j],shapefuncs);
              values[j] = OREAL(0.0);
              for (size_t k = 0; k < 26; k++)
              {
                values[j] += cellvalues[k] * shapefuncs[k];
              }
            }

            int result = Polygonise<OREAL>(&coords[0],values,0,tris,0,0,0,tol);
if (result == 0)
{
  cell.second.data.first = 7;
}

          }
#endif
        } else
        {
          int result = Polygonise<OREAL>(&rconfcoord[0],cellvalues,0,tris,0,0,0,tol);
if (result == 0)
{
  cell.second.data.first = 8;
}

        }
      }
    }
  }

//!!!
etris.saveSTL("wrongtris.stl","Wrong",true);
}
#endif


bool OOctree::cellsToTriangles(const OCELLS &cells, TTriangles<OREAL>& tris)
{
  for (auto c : cells)
  {
    if (isLeaf(c.first))
    {
      std::array<OVECTOR,8> coords = cell8RealCoordinates(c.first);

      for (auto t : faceTris)
      {
        tris.addTri(coords[t[0]],coords[t[1]],coords[t[2]],OREAL(0.0));
      }
    }
  }

  return true;
}

bool OOctree::cellsToTriangles(const OCELLS &cells, const char value, TTriangles<OREAL>& tris)
{
  for (auto c : cells)
  {
    if (isLeaf(c.first))
    {
      if (c.second.data.first == value)
      {
        std::array<OVECTOR,8> coords = cell8RealCoordinates(c.first);

        for (auto t : faceTris)
        {
          tris.addTri(coords[t[0]],coords[t[1]],coords[t[2]],OREAL(0.0));
        }
      }
    }
  }

  return true;
}

void OOctree::getConformingCoordinates(const IPosition& ocentre, IPosition coord[27])
{
  LINT size = intCellSize(level(ocentre));
  assert(size > 1);
  assert(size % 2 == 0);

  assert(ocentre.XYZ[0] % 2 == 0);
  assert(ocentre.XYZ[1] % 2 == 0);
  assert(ocentre.XYZ[2] % 2 == 0);

  LINT halfSize = size / 2;

  for (size_t i = 0; i < 27; i++)
  {
    coord[i] = ocentre + cellCentreOffsets[i] * halfSize;
  }
}

void OOctree::getRealConformingCoordinates(const IPosition& ocentre, OVECTOR coord[27])
{
  IPosition icoords[27];
  getConformingCoordinates(ocentre,icoords);
  for (int j = 0; j < 27; j++)
  {
    coord[j] = intToRealCoord(icoords[j]);
  }
}
