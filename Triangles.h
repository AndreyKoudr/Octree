/**
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
#include "Strings.h"

#include <assert.h>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>

using namespace std;

// output all nodes and indices to console in DEBUG
//#define DEBUG_ALL

/**
  Class TTriangles to load/save STL files
  ---------------------------------------
  https://en.wikipedia.org/wiki/STL_(file_format)

  Simple templated class to keep a surface as a set of 3D triangles represented by node coordinates <I>coords</I>
and connectivity array <I>corners</I> (3 coordinate indices for every triangle (face)). So, the
three coordinates for an i-th triangle can be extracted as

  coord0 = coords[corners[i * 3]];
  coord1 = coords[corners[i * 3 + 1]];
  coord2 = coords[corners[i * 3 + 2]];
  
  Triangles can be loaded and saved from/into STL files. No checks about node numeration :
it is supposed that face normal is defined by counter-clockwise numeration of nodes when
looking from the normal sharp end (see faceNormal()).

  When loading, duplicate nodes are excluded with a tolerance and face corners renumbered.
The tolerance specified in loadSTL() is important, it must be large enough - it is
used in exclusion of duplicate coordinates. If tolerance too large, degenerated triangles 
may appear; if too small - duplicate nodes may not be all excluded.

  STL files in the binary form may contain multiple parts by simply glueing multiple STL files
together; there is a provision in the code for that; not well tested though. Binary STL
files contain floats in only 4-byte format.

  In text STL files the code selects proper number of digits in text representation for
templated float and double classes.

*/

/** To compare edges. */
bool edgeComp(const std::pair<LINT,LINT> &a, const std::pair<LINT,LINT> &b);

/** Comparison function for this edges. */
class EdgeCompare {
public:
  bool operator()(const std::pair<LINT,LINT> &a, const std::pair<LINT,LINT> &b) const
  {
    std::pair<LINT,LINT> aa = a;
    std::pair<LINT,LINT> bb = b;

    if (aa.second < aa.first)
    {
      LINT temp = aa.first; aa.first = aa.second; aa.second = temp;
    }

    if (bb.second < bb.first)
    {
      LINT temp = bb.first; bb.first = bb.second; bb.second = temp;
    }

    if (aa.first < bb.first)
    {
      return true;
    } else if (aa.first > bb.first)
    {
      return false;
    } else
    {
      if (aa.second < bb.second)
      {
        return true;
      } else if (aa.second > bb.second)
      {
        return false;
      } else
      {
        return false;
      }
    }
  }
};

template <class T> class TTriangles {
public:
                              // node coordinates
  std::vector<TVector<T> > coords;
                              // corner indices into vectors, 3 per face (triangle)
  std::vector<LINT> corners;
                              // constructor(s)
  TTriangles() = default;
                              // destructor
  ~TTriangles() = default;
                              // number of triangles (faces)
  size_t numFaces() const;
                              // add flat triangle of three corners
                              // like that from STL file; call coords.reserve()
                              // before adding.
                              // NB to avoid checks (much faster), set tolerance to zero.
                              // angledegtolerance is set to throw off sliver tris
  bool addTri(TVector<T> v0, TVector<T> v1, TVector<T> v2, T tolerance, T angletolerancedeg = T(0.0));
                              // build connectivity array (with unique corners) by exclusion of
                              // duplicate coords and renumbering node indices (corners)
  bool buildConnectivityArray(const T tolerance);
                              // build connectivity array (with unique corners) by exclusion of
                              // duplicate coords and renumbering node indices (corners) and 
                              // reshuffe and cutting array of node data
  template <typename TI> bool buildConnectivityArray(const T tolerance, std::vector<TI>& coorddata)
  {
    if (numFaces() == 0)
      return false;

    assert(coords.size() > 0);
    assert(corners.size() > 0);
    assert(coorddata.size() == coords.size());
                              // make two passes
    for (int j = 0; j < 2; j++)
    {
                              // exclude duplicate nodes
      std::vector<LINT> replacement;
      if (!removeDupNodes(coords,replacement,tolerance))
        return false;
                              // renumber all indices
      for (size_t i = 0; i < corners.size(); i++)
      {
        corners[i] = replacement[corners[i]];
      }

      // this is very important, do not remove
      auto coorddata1 = coorddata;

      for (size_t i = 0; i < replacement.size(); i++)
      {
        coorddata[replacement[i]] = coorddata1[i];
      }
    }

    coorddata.resize(coords.size());

    return true;
  }

                              // get min/max of node coordinates
  std::pair<TVector<T>,TVector<T>> minmax();
                              // get triangles around every corner
  void getCornerTris(std::vector<std::vector<LINT> >& cornerTris);
                              // find neighbour tris in specified number of layers 
                              // around tris already in triNumbers list (must not be empty)
  void getTriNeighbours(const int numLayersAround, 
    const std::vector<std::vector<LINT> >& cornerTris, std::set<LINT> &triNumbers);
                              // get edges and triangles on sides
  void getEdgeTris(std::map<std::pair<LINT,LINT>, std::vector<LINT>,EdgeCompare>& edgeTris);
                              // coords
  std::array<TVector<T>,3> faceCoords(size_t faceNo) const;
                              // get tri normal
  TVector<T> faceNormal(size_t faceNo) const;
                              // get tri corners
  std::array<TVector<T>,3> threeCorners(size_t faceNo) const;
                              // apply Laplace smooth to nodes; coef = 0 - old nodes;
                              // 0.5 - mean; 1.0 - pure Laplace (dangerous, may produce degenarate
                              // tris)
  void smoothLaplace(const T coef);
                             // load tris from STL file, exclude duplicate nodes with
                              // tolerance and renumber corners, throw off sliver tris
  bool loadSTL(const std::string filename, std::string &partname, bool &binary, T tolerance, 
    T angledegtolerance = T(0.0));
                              // save as STL file
  bool saveSTL(const std::string filename, const std::string partname, bool binary) const;
                              // save triangles into OBJ file
  bool saveOBJ(const std::string filename, const std::string partname) const;

                              // generate points from triangle surfaces. pointsPerTri - 
                              // #points per traingle (approximate). Points may be duplicate.
  bool generatePoints(int pointsPerTri, bool excludeDuplicates, T tolerance, 
    std::vector<TVector<T> > &points);
                              // generate points from triangle surfaces. pointsPerTri - 
                              // #points per traingle (approximate). Points may be duplicate.
  bool generatePoints1(T fraction, bool excludeDuplicates, T tolerance, 
    std::vector<TVector<T> > &points);

private:
                              // save as STL file
  bool saveSTL(FILE *fp, std::string partname, bool binary) const;
                              // save triangles into OBJ file
  bool saveOBJ(FILE *fp, const std::string partname) const;
};

template<class T> size_t TTriangles<T>::numFaces() const
{
  assert(corners.size() % 3 == 0);
  return corners.size() / 3;
}

template<class T> bool TTriangles<T>::addTri(TVector<T> v0, TVector<T> v1, TVector<T> v2, T tolerance,
  T angletolerancedeg)
{
if (corners.size() / 3 == 4049)
{
  int gsgsgs = 0;
}
                              // avoid checks
  if (tolerance > T(0.0) || angletolerancedeg > T(0.0))
  {
    TVector<T> v01 = v1 - v0;
    TVector<T> v12 = v2 - v1;
    TVector<T> v20 = v0 - v2;
                              // check degeneration
    T len0 = !v01;
    T len1 = !v12;
    T len2 = !v20;

    if (len0 <= tolerance)
      return false;
    if (len1 <= tolerance)
      return false;
    if (len2 <= tolerance)
      return false;

    T area = !(v01 ^ v12) * T(0.5);
    if (area < tolerance * tolerance)
      return false;

    T angle = (!(v01 ^ v12)) / ((!v01 * !v12)) * PCI;
    if (angle < angletolerancedeg)
      return false;
  }
                              // connectivity array (with coord duplicates)
  corners.push_back(coords.size());
  corners.push_back(coords.size() + 1);
  corners.push_back(coords.size() + 2);
                            // coordinates contain all three nodes
  coords.push_back(v0);
  coords.push_back(v1);
  coords.push_back(v2);

  assert(corners.size() % 3 == 0);

  return true;
}

template<class T> bool TTriangles<T>::saveSTL(FILE *fp, std::string partname, bool binary) const
{
  if (binary)
  {
    size_t error = 0;

    char header[80] = {0};
                            // remove "solid" if any
    std::string name = upCase(partname);
    if (name.substr(0,5) == "SOLID")
    {
                            // spoil first character by replacing by space
      partname[0] = ' ';
    }
                            // limit name length by 80 characters
    if (partname.length() > 80)
      partname.erase(80,std::string::npos);
                            // move string to header
    assert(partname.length() <= 80);
    memmove(header,partname.c_str(),partname.length());
                            // save header
    error = fwrite(header,sizeof(header),1,fp);
    if (error == 0) return false;
                            // save number of triangles
    int numtriangles = static_cast<int>(numFaces());
    assert(sizeof(numtriangles) == 4);
    error = fwrite(&numtriangles,sizeof(numtriangles),1,fp);
    if (error == 0) return false;
                            // save triangles (vectors must be 4-byte floats)
    float v[3] = {0};
    for (size_t i = 0; i < numFaces(); i++)
    {
      size_t i3 = i * 3;
                            // save normal
      TVector<T> normal = faceNormal(i);
      v[0] = static_cast<float>(normal.X);
      v[1] = static_cast<float>(normal.Y);
      v[2] = static_cast<float>(normal.Z);
      error = fwrite(v,sizeof(v),1,fp);
      if (error == 0) return false;
                            // save three triangle nodes
      TVector<T> co[3];
      co[0] = coords[corners[i3]];
      co[1] = coords[corners[i3 + 1]];
      co[2] = coords[corners[i3 + 2]];

      for (int j = 0; j < 3; j++)
      {
        v[0] = static_cast<float>(co[j].X);
        v[1] = static_cast<float>(co[j].Y);
        v[2] = static_cast<float>(co[j].Z);
        error = fwrite(v,sizeof(v),1,fp);
        if (error == 0) return false;
      }
                            // save attribute
      short int attr = 0;
      assert(sizeof(attr) == 2);
      error = fwrite(&attr,sizeof(attr),1,fp);
      if (error == 0) return false;
    }

    return true;
  } else
  {
    int error = 0;

    error = fprintf(fp,"solid %s\n",partname.c_str());
    if (error < 0) return false;

    //assert(coords.size() % 3 == 0);

                            // format string depends only on current float type
    int numdigits = std::numeric_limits<T>::digits10;
    std::string digitsstr = to_string(numdigits);
    std::string formatstr = std::string("vertex ") + std::string("%.") + digitsstr + "e %." + digitsstr + "e %." + digitsstr + "e\n";
    std::string normalstr = std::string("facet normal ") + std::string("%.") + digitsstr + "e %." + digitsstr + "e %." + digitsstr + "e\n";

    for (size_t i = 0; i < numFaces(); i++)
    {
      size_t i3 = i * 3;

      TVector<T> normal = faceNormal(i);

      error = fprintf(fp,normalstr.c_str(),normal.X,normal.Y,normal.Z);
      if (error < 0) return false;
      error = fprintf(fp,"outer loop\n");
      if (error < 0) return false;
                            // save three triangle nodes
      TVector<T> co[3];
      co[0] = coords[corners[i3]];
      co[1] = coords[corners[i3 + 1]];
      co[2] = coords[corners[i3 + 2]];

      error = fprintf(fp,formatstr.c_str(),co[0].X,co[0].Y,co[0].Z);
      if (error < 0) return false;
      error = fprintf(fp,formatstr.c_str(),co[1].X,co[1].Y,co[1].Z);
      if (error < 0) return false;
      error = fprintf(fp,formatstr.c_str(),co[2].X,co[2].Y,co[2].Z);
      if (error < 0) return false;

      error = fprintf(fp,"endloop\n");
      if (error < 0) return false;
      error = fprintf(fp,"endfacet\n");
      if (error < 0) return false;
    }

    error = fprintf(fp,"endsolid %s\n",partname.c_str());
    if (error < 0) return false;

    return true;
  }
}

template<class T> bool TTriangles<T>::saveSTL(const std::string filename, const std::string partname, 
  bool binary) const
{
  FILE *fp = nullptr;
  if (fopen_s(&fp,filename.c_str(),"wb") != 0)
    return false;

  bool res = saveSTL(fp,partname,binary);

  fclose(fp);

  return res;
}

template<class T> bool TTriangles<T>::saveOBJ(FILE *fp, const std::string partname) const
{
                            // we save only corners
  assert(coords.size() > 0);

  int error = 0;
                            // header
  error = fprintf(fp,"# %s\n",partname.c_str());
  if (error < 0) return false;
                            // save coordinates
  error = fprintf(fp,"\n# coordinates\n\n");
  if (error < 0) return false;
                            // format string depends only on current float type
  int numdigits = std::numeric_limits<T>::digits10;
  std::string digitsstr = to_string(numdigits);
  std::string formatstr = std::string("v ") + std::string("%.") + digitsstr + "e %." + digitsstr + "e %." + digitsstr + "e\n";

  for (size_t i = 0; i < coords.size(); i++)
  {
    error = fprintf(fp,formatstr.c_str(),coords[i].X,coords[i].Y,coords[i].Z);
    if (error < 0) return false;
  }
                            // save coordinates
  error = fprintf(fp,"\n# faces\n\n");
  if (error < 0) return false;

  for (size_t i = 0; i < numFaces(); i++)
  {
    size_t i3 = i * 3;
    size_t i0 = corners[i3];
    size_t i1 = corners[i3 + 1];
    size_t i2 = corners[i3 + 2];

    error = fprintf(fp,"f %zd %zd %zd\n",i0 + 1,i1 + 1,i2 + 1);
    if (error < 0) return false;
  }

  return true;
}

template<class T> bool TTriangles<T>::saveOBJ(const std::string filename, const std::string partname) const
{
  FILE *fp = nullptr;
  if (fopen_s(&fp,filename.c_str(),"wb") != 0)
    return false;

  bool res = saveOBJ(fp,partname);

  fclose(fp);

  return res;
}

template<class T> bool TTriangles<T>::loadSTL(const std::string filename, std::string &partname, 
  bool &binary, T tolerance, T angledegtolerance)
{
  FILE *fp = nullptr;
  if (fopen_s(&fp,filename.c_str(),"rb") != 0)
    return false;
                            // result
  bool OK = false;
                            // read first bytes to find out if the file is binary or text
  binary = false;

  char header[200];
  size_t bytesread = fread(header,1,sizeof(header),fp);
  OK = (bytesread > 0);
  if (!OK) 
  {
    fclose(fp);
    return false;
  }
                            // text file will contain "solid" and a number of LFs in 200 bytes
                            // truncate header and use std::string to find "solid"
  header[199] = 0;
  auto sheader = std::string(header);
                          // "solid" starting from zero byte
  bool solidfound = (sheader.find("solid",0,5) == 0) ||
    (sheader.find("SOLID",0,5) == 0);
                          // count number of LFs
  int LFcount = 0;
  for (size_t i = 0; i < bytesread; i++)
  {
    if (header[i] == 10) LFcount++;
  }
                          // 100% the file is ASCII
  binary = !(solidfound && LFcount > 0);
                          // file is binary just seek its beginning
  if (binary)
  {
    partname = std::string(header);
                          // go to start of file
    OK = (fseek(fp,0,SEEK_SET) == 0);
                          // read non-standard file as many files glued together
    while (!feof(fp))
    {
                          // all int-s are 32-bit in STL file
      assert(sizeof(int) == 4);
                          // skip header
      OK = OK && (fseek(fp,80,SEEK_CUR) == 0);
                          // read number of triangles
      int numtriangles = 0;
      OK = OK && (fread(&numtriangles,1,sizeof(numtriangles),fp) == sizeof(numtriangles));
      if (!OK) 
      {
        fclose(fp);
                            // OK, if some triangles already loaded
        bool ok = buildConnectivityArray(tolerance);
        return ok;
      }
      OK = (numtriangles > 0);
      if (!OK) 
      {
        fclose(fp);
                            // OK, if some triangles already loaded
        bool ok = buildConnectivityArray(tolerance);
        return ok;
      }
                          // read all triangles
      assert(sizeof(float) == 4);
      float xyz[3];

      for (int i = 0; i < numtriangles; i++)
      {
        //TVector<T> normal;

        OK = (fread(xyz,1,sizeof(xyz),fp) == sizeof(xyz));
        if (!OK) 
        {
          fclose(fp);
          return false;
        }

        //normal.X = static_cast<T>(xyz[0]);
        //normal.Y = static_cast<T>(xyz[1]);
        //normal.Z = static_cast<T>(xyz[2]);
        //normal = +normal;
                            // three vertices
        TVector<T> v[3];

        for (int j = 0; j < 3; j++)
        {
          OK = (fread(xyz,1,sizeof(xyz),fp) == sizeof(xyz));
          if (!OK) 
          {
            fclose(fp);
            return false;
          }

          v[j].X = static_cast<T>(xyz[0]);
          v[j].Y = static_cast<T>(xyz[1]);
          v[j].Z = static_cast<T>(xyz[2]);
        }

        addTri(v[0],v[1],v[2],tolerance,angledegtolerance);
                            // 2 bytes of unused data
        short int temp = 0;
        OK = (fread(&temp,1,2,fp) == 2);
        if (!OK) 
        {
          fclose(fp);
          return false;
        }
      }
    }

    fclose(fp);

    if (tolerance == 0.0)
    {
      std::pair<TVector<T>,TVector<T>> mm = minmax();
      TVector<T> d = mm.second - mm.first;
      tolerance = !d * 0.000001;
    }
                            // remove duplicate nodes and renumber corners
    bool ok = buildConnectivityArray(tolerance);
    return ok;
  } else
                            // ASCII
  {
                            // close binary file and open as text
    fclose(fp);

    std::ifstream file;
    std::ios_base::openmode oMode(std::ios::in|std::ios::binary);
    file.open(filename.c_str(), oMode);

    if (!file.is_open())
    {
      return false;
    }
                            // line by line...

                            // loop theoretically may contain more than 3 nodes but
                            // it does not seem possible as binary version does not allow
                            // it at all
    int nodecount = 0;
    TVector<T> normal;
    TVector<T> v[3];

    while (!file.eof())
    {
      std::string line;
      std::getline(file,line);
                            // last line might well be without ending (CR)/LF
      if (file.fail()) {
        break;
      }
                            // trim line
      line = trim(line," \n\r\t");
                            // skip empty
      if (line.length() == 0)
        continue;
                            // original line
      std::string original = line;
                            // uppercase
      line = upCase(line);
                            // lines starting with SOLID
      if (line.substr(0,5) == "SOLID")
      {
        if (original.length() > 6)
        {
                            // store part name
          partname = original.substr(6);
          partname = trim(partname," \n\r\t");
        }
        continue;
      }
                            // lines starting with ENDSOLID
      if (line.substr(0,8) == "ENDSOLID")
      {
        continue;
      }

      if (line == "ENDFACET")
      {
      } else if (line.substr(0,5) == "FACET")
      {
                            // normal for this facet, this line is expected to be like
                            // "facet normal -0.7074 0.0 0.7074"
        int pos1[100],pos2[100];
        int numwords = parseWords(line,' ',pos1,pos2,100);
                            // skip this line
        if (numwords != 5)
        {
          normal = TVector<T>();
                            // error, issue warning, do not return
          continue;
        }
                            // do not parse normal

      } else if (line == "OUTER LOOP")
      {
        nodecount = 0;

      } else if (line == "ENDLOOP")
      {
        if (nodecount >= 3)
        {
                            // add triangle
          addTri(v[0],v[1],v[2],tolerance,angledegtolerance);
        } else
        {
                            // Wrong count of nodes
                            // error, issue warning, do not return
        }

                            // "vertex -1.0 0.5 0.33"
      } else if (line.substr(0,6) == "VERTEX")
      {
        int pos1[100],pos2[100];
        int numwords = parseWords(line,' ',pos1,pos2,100);
                            // "vertex" line contains less than 3 coordinates, take only
                            // first available coordinates, others left zeroes
        if (nodecount < 3)
        {
          if (numwords < 4)
          {
            if (numwords == 2)
            {
              v[nodecount].X = static_cast<T>(atof(line.substr(pos1[1],(pos2[1] - pos1[1] + 1)).c_str()));
              v[nodecount].Y = T(0.0);
              v[nodecount].Z = T(0.0);
            } if (numwords == 3)
            {
              v[nodecount].X = static_cast<T>(atof(line.substr(pos1[1],(pos2[1] - pos1[1] + 1)).c_str()));
              v[nodecount].Y = static_cast<T>(atof(line.substr(pos1[2],(pos2[2] - pos1[2] + 1)).c_str()));
              v[nodecount].Z = T(0.0);
            }
                              // error, issue warning, do not return
          } else
          {
            v[nodecount].X = static_cast<T>(atof(line.substr(pos1[1],(pos2[1] - pos1[1] + 1)).c_str()));
            v[nodecount].Y = static_cast<T>(atof(line.substr(pos1[2],(pos2[2] - pos1[2] + 1)).c_str()));
            v[nodecount].Z = static_cast<T>(atof(line.substr(pos1[3],(pos2[3] - pos1[3] + 1)).c_str()));
          }

          nodecount++;
        }
      } else

                            // unidentified line
      {
                            // error, issue warning, do not return
      }
    } // while loop on all lines

    file.close();

    if (tolerance == 0.0)
    {
      std::pair<TVector<T>,TVector<T>> mm = minmax();
      TVector<T> d = mm.second - mm.first;
      tolerance = !d * 0.000001;
    }

                            // remove duplicate nodes and renumber corners
    bool ok = buildConnectivityArray(tolerance);
    return ok;

  } // ASCII

  return false;
}
                              // tolerance for comparison routines
static double gtolerance2 = TOLERANCE(float) * TOLERANCE(float);

                              // comparator for sort()
template <class T> bool comparePoint(TVector<T> p1, TVector<T> p2)
{
  if (p1.X < p2.X)
  {
    return true;
  } else if (p1.X > p2.X)
  {
    return false;
  } else
  {
    if (p1.Y < p2.Y)
    {
      return true;
    } else if (p1.Y > p2.Y)
    {
      return false;
    } else
    {
      if (p1.Z < p2.Z)
      {
        return true;
      } else if (p1.Z > p2.Z)
      {
        return false;
      } else
      {
        return false;
      }
    }
  }
}
                              // comparator for equality
template <class T> bool equalPoint(TVector<T> p1, TVector<T> p2)
{
  T dist2 = (p1 - p2).length2();
  return (dist2 < gtolerance2);
}
                            // exclude duplicates from 3D array
template <class T> bool removeDupNodes(std::vector<TVector<T> > &vectors, 
  std::vector<LINT> &replacement, T tolerance)
{
                              // list is empty
  if (vectors.size() < 1) return true;
                              // just one vector
  if (vectors.size() == 1)
  {
    replacement.push_back(0);
    return true;
  }
                              // set tolerance for comparison
  gtolerance2 = static_cast<double>(tolerance * tolerance);
                              // put index number in W
  for (size_t i = 0; i < vectors.size(); i++)
  {
    vectors[i].W = static_cast<T>(i);
  }
                              // replacement after sort
  std::vector<LINT> replacement0(vectors.size(),-1);
                              // erase
  std::sort(vectors.begin(), vectors.end(), comparePoint<T>);

  for (size_t i = 0; i < vectors.size(); i++)
  {
    size_t index = ROUND(vectors[i].W);
    replacement0[i] = index;
  }
                              // all replacements must be filled up
#ifdef _DEBUG
  for (auto r : replacement0)
  {
    assert(r >= 0);
  }
#endif
                              // make replacement
  size_t count = 0;
  vectors[0].W = static_cast<T>(count);
  for (auto v = vectors.begin() + 1; v != vectors.end(); v++)
  {
    auto vprev = v - 1;
    if (!equalPoint<T>(*vprev,*v))
    {
      count++;
    }
    v->W = static_cast<T>(count);
  }

  replacement.clear();
  replacement.resize(vectors.size(),-1);
  for (size_t i = 0; i < vectors.size(); i++)
  {
    size_t index = ROUND(vectors[i].W);
    replacement[replacement0[i]] = index;
  }
                              // all replacements must be filled up
#ifdef _DEBUG
  for (auto r : replacement)
  {
    assert(r >= 0);
  }
#endif
                              // exclude
  auto unique_end = std::unique(vectors.begin(), vectors.end(), equalPoint<T>);
  vectors.erase(unique_end, vectors.end());

  return true;
}

template <class T> bool TTriangles<T>::buildConnectivityArray(const T tolerance)
{
  if (numFaces() == 0)
    return false;

  assert(coords.size() > 0);
  assert(corners.size() > 0);
  //assert(tolerance >= 0.00000001);
                              // make two passes
  for (int j = 0; j < 2; j++)
  {
                              // exclude duplicate nodes
    std::vector<LINT> replacement;
    if (!removeDupNodes(coords,replacement,tolerance))
      return false;
                              // renumber all indices
    for (size_t i = 0; i < corners.size(); i++)
    {
      corners[i] = replacement[corners[i]];
    }
  }

  // test for correct node numbering : coincident nodes indicate that
  // tolarance is TOO high
  for (size_t i = 0; i < corners.size(); i += 3)
  {
    if (corners[i] == corners[i + 1] || corners[i] == corners[i + 2] || corners[i + 1] == corners[i + 2])
      return false;
  }


#ifdef DEBUG_ALL
#ifdef _DEBUG
  printf("\nNODES\n");
  for (size_t i = 0; i < coords.size(); i++)
  {
    printf("%lu %f %f %f\n",i,coords[i].X,coords[i].Y,coords[i].Z);
  }
  printf("\nCORNERS\n");
  for (size_t i = 0; i < corners.size(); i += 3)
  {
    printf("%lu %zd %zd %zd\n",i,corners[i],corners[i + 1],corners[i + 2]);
    assert(corners[i] != corners[i + 1]);
    assert(corners[i] != corners[i + 2]);
    assert(corners[i + 1] != corners[i + 2]);
  }
#endif
#endif

  return true;
}

template <class T> std::pair<TVector<T>,TVector<T>> TTriangles<T>::minmax()
{
  assert(coords.size() > 0);
  TVector<T> min = coords[0];
  TVector<T> max = coords[0];

  for (size_t i = 1; i < coords.size(); i++)
  {
    for (size_t j = 0; j < 3; j++)
    {
      min.XYZ[j] = std::min(min.XYZ[j],coords[i].XYZ[j]);
      max.XYZ[j] = std::max(max.XYZ[j],coords[i].XYZ[j]);
    }
  }

  return std::pair<TVector<T>,TVector<T>>(min,max);
}

template <class T> void TTriangles<T>::getCornerTris(std::vector<std::vector<LINT> >& cornerTris)
{
  assert(coords.size() > 0);
  cornerTris.clear();
  cornerTris.resize(coords.size(),std::vector<LINT>());

  for (size_t i = 0; i < corners.size(); i += 3)
  {
    size_t trinumber = i / 3;

    for (size_t j = 0; j < 3; j++)
    {
      cornerTris[corners[i + j]].push_back(trinumber);
    }
  }
}

template <class T> void TTriangles<T>::getEdgeTris(std::map<std::pair<LINT,LINT>,
  std::vector<LINT>,EdgeCompare>& edgeTris)
{
  assert(coords.size() > 0);
  edgeTris.clear();

  for (size_t i = 0; i < corners.size(); i += 3)
  {
    LINT trinumber = i / 3;

    std::pair<LINT,LINT> e[3] = {
      std::pair<LINT,LINT>(corners[i],corners[i + 1]),
      std::pair<LINT,LINT>(corners[i + 1],corners[i + 2]),
      std::pair<LINT,LINT>(corners[i + 2],corners[i])
    };

    for (size_t j = 0; j < 3; j++)
    {
      // swap numbers
      if (e[j].second < e[j].first)
      {
        LINT temp = e[j].second; e[j].second = e[j].first; e[j].first = temp;
      }

      auto it = edgeTris.find(e[j]);
      if (it == edgeTris.end())
      {
        edgeTris.insert(std::pair<std::pair<LINT,LINT>,std::vector<LINT>>(e[j],
          std::vector<LINT>(std::initializer_list<LINT>{trinumber})));
      } else
      {
        it->second.push_back(trinumber);
      }
    }
  }
}

template <class T> void TTriangles<T>::smoothLaplace(const T coef)
{
  std::vector<std::vector<LINT> > cornerTris;
  getCornerTris(cornerTris);

  std::vector<TVector<T> > newCoords(coords.size());

  for (size_t i = 0; i < coords.size(); i++)
  {
    std::set<LINT> neiNodes;
    for (auto tri : cornerTris[i])
    {
      neiNodes.insert(corners[tri * 3]);
      neiNodes.insert(corners[tri * 3 + 1]);
      neiNodes.insert(corners[tri * 3 + 2]);
    }

    TVector<T> sum(0,0,0);
    int count = 0;
    for (auto n : neiNodes)
    {
      if (n != i)
      {
        sum += coords[n];
        count++;
      }
    }

    newCoords[i] = coords[i] + ((sum /= static_cast<T>(count)) - coords[i]) * coef;

#if 0
#ifdef _DEBUG
    T diff = !(coords[i] - newCoords[i]);
    cout << diff << endl;
#endif
#endif
  }

  coords = newCoords;
}

template <class T> TVector<T> TTriangles<T>::faceNormal(size_t faceNo) const
{
  LINT i0 = corners[faceNo * 3];
  LINT i1 = corners[faceNo * 3 + 1];
  LINT i2 = corners[faceNo * 3 + 2];
  TVector<T> n = +((coords[i1] - coords[i0]) ^ (coords[i2] - coords[i1]));
  return n;
}

template <class T> std::array<TVector<T>,3> TTriangles<T>::threeCorners(size_t faceNo) const
{
  LINT i0 = corners[faceNo * 3];
  LINT i1 = corners[faceNo * 3 + 1];
  LINT i2 = corners[faceNo * 3 + 2];
  return std::array<TVector<T>,3>({coords[i0],coords[i1],coords[i2]});
}

template <class T> std::array<TVector<T>,3> TTriangles<T>::faceCoords(size_t faceNo) const
{
  LINT i0 = corners[faceNo * 3];
  LINT i1 = corners[faceNo * 3 + 1];
  LINT i2 = corners[faceNo * 3 + 2];
  std::array<TVector<T>,3> r;
  r[0] = coords[i0];
  r[1] = coords[i1];
  r[2] = coords[i2];
  return r;
}

template <class T> void TTriangles<T>::getTriNeighbours(const int numLayersAround, 
  const std::vector<std::vector<LINT> >& cornerTris, std::set<LINT> &triNumbers)
{
  // trinumbers should not be empty (must conaint at least one starting tri number)
  assert(triNumbers.size() > 0);

  for (size_t i = 0; i < numLayersAround; i++)
  {
    std::set<LINT> newtrinumbers;
    for (auto tri : triNumbers)
    {
      for (size_t j = 0; j < 3; j++)
      {
        LINT corner = corners[tri * 3 + j];
        for (auto it : cornerTris[corner])
        {
          newtrinumbers.insert(it);
        }
      }
    }
    triNumbers.insert(newtrinumbers.begin(),newtrinumbers.end());
  }
}

template <class T> bool TTriangles<T>::generatePoints(int pointsPerTri, bool excludeDuplicates, 
  T tolerance, std::vector<TVector<T> > &points)
{
  int numSteps = (int) sqrt(double(pointsPerTri));
  LIMIT_MIN(numSteps,1);
                            // got through all tris
  for (size_t i = 0; i < numFaces(); i++)
  {
     size_t i3 = i * 3;
                            // three triangle nodes
    TVector<T> co[3];
    co[0] = coords[corners[i3]];
    co[1] = coords[corners[i3 + 1]];
    co[2] = coords[corners[i3 + 2]];

    T d = 1.0 / T(numSteps);
    for (int j = 0; j < numSteps + 1; j++)
    {
      T u = T(j) * d;
      for (int k = 0; k < numSteps + 1; k++)
      {
        T v = T(k) * d;
        T w = 1.0 - u - v;

        if (w >= 0.0)
        {
          TVector<T> p = co[0] * u + co[1] * v + co[2] * w;
          points.push_back(p);
        }
      }
    }
  }

  if (excludeDuplicates)
  {
    if (tolerance == 0.0)
    {
      std::pair<TVector<T>,TVector<T>> mm = minmax();
      TVector<T> d = mm.second - mm.first;
      tolerance = !d * 0.000001;
    }
                            // remove duplicate nodes and renumber corners
    std::vector<LINT> replacement;
    removeDupNodes(points,replacement,tolerance);
  }

  return true;
}


template <class T> void refineTriangle(const TVector<T> corner0, 
  const TVector<T> corner1, const TVector<T> corner2, 
  const T fraction, std::vector<TVector<T>> &points)
{
  // for every corner, two other corners
  static int sides[3][2] = {{1,2},{2,0},{0,1}};

  // corners
  TVector<T> c[3] = {corner0, corner1, corner2};

  // num divisions of all sides
  int divs[3];
  // side length in front of every corner
  T len[3];
  for (int i = 0; i < 3; i++)
  {
    len[i] = !(c[sides[i][1]] - c[sides[i][0]]);
    divs[i] = int(len[i] / fraction);
    LIMIT_MIN(divs[i],1);
  }

  // the shortest side?
  int start = 0; 
  T startLen = len[0];
  for (int i = 1; i < 3; i++)
  {
    if (len[i] < startLen)
    {
      startLen = len[i];
      start = i;
    }
  }

  // now go along two longest side in front of the shortest
  int ds = std::max<int>(divs[sides[start][0]],divs[sides[start][1]]);

  TVector<T> vstart = c[start];
  TVector<T> vend0 = c[sides[start][0]];
  TVector<T> vend1 = c[sides[start][1]];

  points.push_back(vstart);

  for (int i = 1; i <= ds; i++)
  {
    T U = T(i) / T(ds);

    TVector<T> v0 = vstart + (vend0 - vstart) * U;
    TVector<T> v1 = vstart + (vend1 - vstart) * U;

    // cross direction
    T clen = !(v1 - v0);
    int cdivs = int(clen / fraction);
    LIMIT_MIN(cdivs,1);

    for (int j = 0; j <= cdivs; j++)
    {
      T V = T(j) / T(cdivs);
      TVector<T> v = v0 + (v1 - v0) * V;
      points.push_back(v);
    }
  }
}

template <class T> bool TTriangles<T>::generatePoints1(T fraction, bool excludeDuplicates, 
  T tolerance, std::vector<TVector<T> > &points)
{
                            // got through all tris
  for (size_t i = 0; i < numFaces(); i++)
  {
     size_t i3 = i * 3;
                            // three triangle nodes
    TVector<T> co[3];
    co[0] = coords[corners[i3]];
    co[1] = coords[corners[i3 + 1]];
    co[2] = coords[corners[i3 + 2]];

    refineTriangle(co[0],co[1],co[2],fraction,points);
  }

  if (excludeDuplicates)
  {
    if (tolerance == 0.0)
    {
      std::pair<TVector<T>,TVector<T>> mm = minmax();
      TVector<T> d = mm.second - mm.first;
      tolerance = !d * 0.000001;
    }
                            // remove duplicate nodes and renumber corners
    std::vector<LINT> replacement;
    removeDupNodes(points,replacement,tolerance);
  }

  return true;
}
