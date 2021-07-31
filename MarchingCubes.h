#pragma once

#include <vector>
#include "Vector.h"
#include "Triangles.h"

/*

  Taken from
  http://paulbourke.net/geometry/polygonise/

*/

                              // static data for marching cubes
extern int edgeTable[256];
extern int triTable[256][16];

/*
   Linearly interpolate the position where an isosurface cuts
   an edge between two vertices, each with their own scalar value
*/

template <class T> TVector<T> VertexInterp(T isolevel, TVector<T> p1, TVector<T> p2, 
  T valp1, T valp2, T tolerance)
{
#if 1
  assert(isolevel == isolevel);

  assert(valp1 == valp1);

  assert(valp2 == valp2);
  

  if (std::abs(isolevel - valp1) < tolerance)
    return(p1);
  if (std::abs(isolevel - valp2) < tolerance)
    return(p2);
  if (std::abs(valp1 - valp2) < tolerance)
    return((p1 + p2) * 0.5); //!!!
//          return(p1);

  //assert(std::abs(valp2 - valp1) > tolerance);
  T coef = (isolevel - valp1) / (valp2 - valp1);

  assert(coef == coef);

  TVector<T> p = p1 + (p2 - p1) * coef;

  assert(p[0] == p[0]);
  assert(p[1] == p[1]);
  assert(p[2] == p[2]);

  return p;
#else
   if (std::abs(isolevel - valp1) < tolerance)
      return(p1);
   if (std::abs(isolevel - valp2) < tolerance)
      return(p2);
   if (std::abs(valp1 - valp2) < tolerance)
      return(p1);

   T mu = (isolevel - valp1) / (valp2 - valp1);

   TVector<T> p = p1 + (p2 - p1) * mu;

   return p;
#endif
}

/*
   Given a grid cell and an isolevel, calculate the triangular
   facets required to represent the isosurface through the cell.
   Return the number of triangular facets, the array "triangles"
   will be loaded up with the vertices at most 5 triangular facets.
	 0 will be returned if the grid cell is either totally above
   of totally below the isolevel.
*/
template <class T> int Polygonise(TVector<T> cellcoords[8], T cellvalues[8], 
  T isolevel, TTriangles<T> &triangles, T limit0, T limit1, T limit2, T limittolerance)
{
#if 0
printf("v %e %fe %e %e %e %e %e %e\n",cellvalues[0],cellvalues[1],cellvalues[2],cellvalues[3],
  cellvalues[4],cellvalues[5],cellvalues[6],cellvalues[7]);
#endif


  TVector<T> vertlist[12];

//  T tolerance = TOLERANCE;
  T tolerance = limittolerance;

  /*
    Adjust cell values to avoid holes at surface where potential == limit0,
      limit1 or limit2
  */
#if 1
  for (int i = 0; i < 8; i++)
  {
    T d = std::abs(cellvalues[i] - limit0);
    if (d < limittolerance)
    {
      cellvalues[i] = limit0;
      continue;
    }
    d = std::abs(cellvalues[i] - limit1);
    if (d < limittolerance)
    {
      cellvalues[i] = limit1;
      continue;
    }
    d = std::abs(cellvalues[i] - limit2);
    if (d < limittolerance)
    {
      cellvalues[i] = limit2;
      continue;
    }
  }
#endif

  /*
    Determine the index into the edge table which
    tells us which vertices are inside of the surface
  */
  int cubeindex = 0;
#if 0
  if (cellvalues[0] < isolevel ) cubeindex |= 1;
  if (cellvalues[1] < isolevel ) cubeindex |= 2;
  if (cellvalues[2] < isolevel ) cubeindex |= 4;
  if (cellvalues[3] < isolevel ) cubeindex |= 8;
  if (cellvalues[4] < isolevel ) cubeindex |= 16;
  if (cellvalues[5] < isolevel ) cubeindex |= 32;
  if (cellvalues[6] < isolevel ) cubeindex |= 64;
  if (cellvalues[7] < isolevel ) cubeindex |= 128;
#else
  if (cellvalues[0] <= isolevel + tolerance) cubeindex |= 1;
  if (cellvalues[1] <= isolevel + tolerance) cubeindex |= 2;
  if (cellvalues[2] <= isolevel + tolerance) cubeindex |= 4;
  if (cellvalues[3] <= isolevel + tolerance) cubeindex |= 8;
  if (cellvalues[4] <= isolevel + tolerance) cubeindex |= 16;
  if (cellvalues[5] <= isolevel + tolerance) cubeindex |= 32;
  if (cellvalues[6] <= isolevel + tolerance) cubeindex |= 64;
  if (cellvalues[7] <= isolevel + tolerance) cubeindex |= 128;
#endif

  /* Cube is entirely in/out of the surface */
  if (edgeTable[cubeindex] == 0)
    return 0;

  /* Find the vertices where the surface intersects the cube */
  if (edgeTable[cubeindex] & 1)
    vertlist[0] =
       VertexInterp<T>(isolevel,cellcoords[0],cellcoords[1],cellvalues[0],cellvalues[1],tolerance);
  if (edgeTable[cubeindex] & 2)
    vertlist[1] =
       VertexInterp<T>(isolevel,cellcoords[1],cellcoords[2],cellvalues[1],cellvalues[2],tolerance);
  if (edgeTable[cubeindex] & 4)
    vertlist[2] =
       VertexInterp<T>(isolevel,cellcoords[2],cellcoords[3],cellvalues[2],cellvalues[3],tolerance);
  if (edgeTable[cubeindex] & 8)
    vertlist[3] =
       VertexInterp<T>(isolevel,cellcoords[3],cellcoords[0],cellvalues[3],cellvalues[0],tolerance);
  if (edgeTable[cubeindex] & 16)
    vertlist[4] =
       VertexInterp<T>(isolevel,cellcoords[4],cellcoords[5],cellvalues[4],cellvalues[5],tolerance);
  if (edgeTable[cubeindex] & 32)
    vertlist[5] =
       VertexInterp<T>(isolevel,cellcoords[5],cellcoords[6],cellvalues[5],cellvalues[6],tolerance);
  if (edgeTable[cubeindex] & 64)
    vertlist[6] =
       VertexInterp<T>(isolevel,cellcoords[6],cellcoords[7],cellvalues[6],cellvalues[7],tolerance);
  if (edgeTable[cubeindex] & 128)
    vertlist[7] =
       VertexInterp<T>(isolevel,cellcoords[7],cellcoords[4],cellvalues[7],cellvalues[4],tolerance);
  if (edgeTable[cubeindex] & 256)
    vertlist[8] =
       VertexInterp<T>(isolevel,cellcoords[0],cellcoords[4],cellvalues[0],cellvalues[4],tolerance);
  if (edgeTable[cubeindex] & 512)
    vertlist[9] =
       VertexInterp<T>(isolevel,cellcoords[1],cellcoords[5],cellvalues[1],cellvalues[5],tolerance);
  if (edgeTable[cubeindex] & 1024)
    vertlist[10] =
       VertexInterp<T>(isolevel,cellcoords[2],cellcoords[6],cellvalues[2],cellvalues[6],tolerance);
  if (edgeTable[cubeindex] & 2048)
    vertlist[11] =
       VertexInterp<T>(isolevel,cellcoords[3],cellcoords[7],cellvalues[3],cellvalues[7],tolerance);

                              // add triangles
  int ntriang = 0;
  for (int i = 0; triTable[cubeindex][i] != -1; i += 3) {
    triangles.addTri(
      vertlist[triTable[cubeindex][i]],
      vertlist[triTable[cubeindex][i + 1]],
      vertlist[triTable[cubeindex][i + 2]],
      tolerance);
    ntriang++;
  }

  return ntriang;
}


