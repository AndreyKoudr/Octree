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
#include "Vector.h"

#include <limits>
#include <algorithm>
#include <array>
#include <assert.h>

using namespace std;


/** Very important class to define integer positions of cells and nodes within octree. */
#define IPosition TVector<LINT>

/** 
  Background parameters for a balanced octree.

  Background is a set of largest cells obtained by uniform division of the whole cuboid region.
  Background cells are cells of level 0. Every cell of level 0 can be subdivided (refined)
to make a hierarchy of cells inside a background cell. 

  A search of every 3D point inside octree is two-step :
  (1) which background cell? 
  (2) starting from this background cell, search within a hierarchy
Both operations are fast. This is the way an octree provides a kind of indexing for 3D points.

  T is float type (float or double) for real coordinates.
*/

template <typename T> class OBackground {
public:

  OBackground() = delete;

  /** Contructor */
  OBackground(const TVector<T> &min, const TVector<T> &max, 
    const LINT Im, const LINT Jm, const LINT Km, const int maxOctreeLevel);

  /** Return 8 boundary nodes */
  std::array<TVector<T>,8> boundaryNodes();

  /** Convert back cell I,J,K (not normal integer!) position into back cell number */
  size_t positionToBackCellIndex(const IPosition &position);

  /** Convert back cell number into I,J,K (not normal integer!) position */
  IPosition backCellIndexToPosition(const size_t index);

  /** Num background cells in I,J,K directions. */
  IPosition IJKnumCells = IPosition(0,0,0);

  /** Size of whole region in I,J,K directions in integer units. */
  IPosition IJKintSizes = IPosition(0,0,0);

  /** Box X,Y,Z min coordinate. */
  TVector<T> boxMin = TVector<T>(T0,T0,T0);
  /** Box X,Y,Z max coordinate. */
  TVector<T> boxMax = TVector<T>(T0,T0,T0);

  /** Max octree level. */
  int maxLevel = 0;

  /** Real sizes of integer unit in 3 directions. */
  TVector<T> intUnitSizes = TVector<T>(T0,T0,T0);

  /** Cell size in I,J,K directions. */
  TVector<T> cellSize = TVector<T>(T0,T0,T0);

  /** Max box size. */
  T scale = T0;
};

template <typename T> OBackground<T>::OBackground(const TVector<T> &min, const TVector<T> &max, 
  const LINT Im, const LINT Jm, const LINT Km, const int maxOctreeLevel)
{
  assert(Im > 0);
  assert(Jm > 0);
  assert(Km > 0);
  assert(max.XYZ[0] > min.XYZ[0]);
  assert(max.XYZ[1] > min.XYZ[1]);
  assert(max.XYZ[2] > min.XYZ[2]);

  IJKnumCells = IPosition(Im,Jm,Km);
  maxLevel = maxOctreeLevel;

  LINT l = 1i64 << maxLevel;
  IJKintSizes = IJKnumCells * l;

  assert(IJKintSizes.XYZ[0] > 0);
  assert(IJKintSizes.XYZ[1] > 0);
  assert(IJKintSizes.XYZ[2] > 0);

  TVector<T> d = max - min;

  cellSize.XYZ[0] = d.XYZ[0] / static_cast<T>(Im);
  cellSize.XYZ[1] = d.XYZ[1] / static_cast<T>(Jm);
  cellSize.XYZ[2] = d.XYZ[2] / static_cast<T>(Km);

#ifndef NDEBUG
  T tolerance = std::numeric_limits<T>::epsilon() * static_cast<T>(10.0);

  assert(cellSize.XYZ[0] > tolerance);
  assert(cellSize.XYZ[1] > tolerance);
  assert(cellSize.XYZ[2] > tolerance);
#endif
  
  boxMin = min;
  boxMax = max;

  intUnitSizes.XYZ[0] = cellSize.XYZ[0] * pow(T05,maxLevel);
  intUnitSizes.XYZ[1] = cellSize.XYZ[1] * pow(T05,maxLevel);
  intUnitSizes.XYZ[2] = cellSize.XYZ[2] * pow(T05,maxLevel);

  scale = std::max<T>(d.XYZ[0],std::max<T>(d.XYZ[1],d.XYZ[2]));
}

template <typename T> std::array<TVector<T>,8> OBackground<T>::boundaryNodes()
{
  std::array<TVector<T>,8> nodes = {
    TVector<T>(boxMin.XYZ[0],boxMin.XYZ[1],boxMin.XYZ[2]),
    TVector<T>(boxMax.XYZ[0],boxMin.XYZ[1],boxMin.XYZ[2]),
    TVector<T>(boxMax.XYZ[0],boxMax.XYZ[1],boxMin.XYZ[2]),
    TVector<T>(boxMin.XYZ[0],boxMax.XYZ[1],boxMin.XYZ[2]),
    TVector<T>(boxMin.XYZ[0],boxMin.XYZ[1],boxMax.XYZ[2]),
    TVector<T>(boxMax.XYZ[0],boxMin.XYZ[1],boxMax.XYZ[2]),
    TVector<T>(boxMax.XYZ[0],boxMax.XYZ[1],boxMax.XYZ[2]),
    TVector<T>(boxMin.XYZ[0],boxMax.XYZ[1],boxMax.XYZ[2])
  };

  return nodes;
}

template <typename T> size_t OBackground<T>::positionToBackCellIndex(const IPosition &position)
{
  assert(position.XYZ[0] >= 0 && position.XYZ[0] < IJKnumCells.XYZ[0]);
  assert(position.XYZ[1] >= 0 && position.XYZ[1] < IJKnumCells.XYZ[1]);
  assert(position.XYZ[2] >= 0 && position.XYZ[2] < IJKnumCells.XYZ[2]);

  size_t index = static_cast<size_t>(position.XYZ[0]) * 
      (static_cast<size_t>(IJKnumCells.XYZ[1]) * static_cast<size_t>(IJKnumCells.XYZ[2])) +
    static_cast<size_t>(position.XYZ[1]) * 
      static_cast<size_t>(IJKnumCells.XYZ[2]) +
    static_cast<size_t>(position.XYZ[2]);

  return index;
}

template <typename T> IPosition OBackground<T>::backCellIndexToPosition(const size_t index)
{
  assert(index < IJKnumCells.XYZ[0] * IJKnumCells.XYZ[1] * IJKnumCells.XYZ[2]);

  size_t rest = index;
  size_t x = index / (IJKnumCells.XYZ[1] * IJKnumCells.XYZ[2]);
  rest -= x * (IJKnumCells.XYZ[1] * IJKnumCells.XYZ[2]);
  size_t y = rest / IJKnumCells.XYZ[2];
  rest -= y * IJKnumCells.XYZ[2];
  size_t z = rest;

  return IPosition(x,y,z);
}
