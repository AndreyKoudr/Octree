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
           
#include <vector>     
#include <fstream>
#include <iostream>              
#include "Vector.h"
#include "Strings.h"
#include "Triangles.h"

template <class T> class CloudPoints : public std::vector<TVector<T> > {
public:
  /** Constructor */
  CloudPoints();

  /** Read from text XYZ file. */
  bool readXYZ(const std::string &fileName);

  /** Save points as XYZ text file. */
  bool writeXYZ(const std::string &fileName);

  /** Read coordinates from STL file. Approximate number of generated points
    per triangle is pointsPerTri. Duplicate points can be excluded. If pointsPerTri is zero,
    just nodal coordinates are copied from STL. */
  bool readSTL(const std::string &fileName, int pointsPerTri, bool excludeDuplicates);
  bool readSTL1(const std::string &fileName, T fractionSize, bool excludeDuplicates);

  /** Calculate min/max values. */
  std::pair<TVector<T>,TVector<T>> calcMinMax();
};

template<class T> CloudPoints<T>::CloudPoints() : std::vector<TVector<T> >()
{
}

template<class T> bool CloudPoints<T>::readXYZ(const std::string &fileName)
{
  std::ifstream file;
  std::ios_base::openmode oMode(std::ios::in|std::ios::binary);
  file.open(fileName.c_str(), oMode);

  if (!file.is_open())
  {
    return false;
  }

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

    int pos1[100],pos2[100];
    int numwords = parseWords(line,' ',pos1,pos2,100);
                        // skip this line
    if (numwords != 3)
    {
                        // error, skip this line
      continue;
    }

    TVector<T> v(0,0,0,0);
    v.X = static_cast<T>(atof(line.substr(pos1[0],(pos2[0] - pos1[0] + 1)).c_str()));
    v.Y = static_cast<T>(atof(line.substr(pos1[1],(pos2[1] - pos1[1] + 1)).c_str()));
    v.Z = static_cast<T>(atof(line.substr(pos1[2],(pos2[2] - pos1[2] + 1)).c_str()));

    std::vector<TVector<T> >::push_back(v);
  }

  file.close();

  return true;
}

template<class T> bool CloudPoints<T>::readSTL(const std::string &fileName, int pointsPerTri, 
  bool excludeDuplicates)
{
  TTriangles<T> tris;
  std::string partname;
  bool binary;
  if (tris.loadSTL(fileName,partname,binary,T(0.0)))
  {
    if (pointsPerTri < 1)
    {
      *((std::vector<TVector<T>> *) this) = tris.coords;
    } else
    {
      tris.generatePoints(pointsPerTri,excludeDuplicates,0.0,*((std::vector<TVector<T>> *) this));
    }
    return true;
  } else
  {
    return false;
  }
}

template<class T> bool CloudPoints<T>::readSTL1(const std::string &fileName, T fractionSize, 
  bool excludeDuplicates)
{
  TTriangles<T> tris;
  std::string partname;
  bool binary;
  if (tris.loadSTL(fileName,partname,binary,T(0.0)))
  {
    std::pair<TVector<T>,TVector<T> > minmax = tris.minmax();
    T fsize = !(minmax.second - minmax.first) * fractionSize;

    tris.generatePoints1(fsize,excludeDuplicates,0.0,*((std::vector<TVector<T>> *) this));
    return true;
  } else
  {
    return false;
  }
}

template<class T> std::pair<TVector<T>,TVector<T>> CloudPoints<T>::calcMinMax()
{
  assert(std::vector<TVector<T> >::size() > 0);

  TVector<T> vmin(0,0,0,0), vmax(0,0,0,0);

  if (std::vector<TVector<T> >::size() > 0)
  {
    vmin = vmax = std::vector<TVector<T> >::front();
    for (auto v : *this)
    {
      for (size_t j = 0; j < 3; j++)
      {
        vmin[j] = std::min(vmin[j],v[j]);
        vmax[j] = std::max(vmax[j],v[j]);
      }
    }
  }

  return std::pair<TVector<T>,TVector<T>>(vmin,vmax);
}

template<class T> bool CloudPoints<T>::writeXYZ(const std::string &fileName)
{
  std::ofstream file;
  file.open(fileName.c_str());

  if (!file.is_open())
  {
    return false;
  }

  file.setf(ios::fixed,ios::floatfield);
  file.precision(14);

  for (TVector<T> v : *this)
  {
    file << v.X << " " << v.Y << " " << v.Z << "\n"; 
  }

  file.close();

  return true;
}
