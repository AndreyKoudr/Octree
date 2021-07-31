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

using namespace std;

/**
  Sub cells :
            
                           *---------------------* 
                          /    3     /     7    /|    
                         /---------------------/ |     
                        /    1     /     5    /|7|          
                       *---------------------* | |    
                       |          |          |5 /|     
                       |     1    |    5     | / |    
                       |          |          |/|6|
                       |----------|----------| | * 
       K               |          |          |4|/ 
       | J             |     0    |    4     | /
       |/              |          |          |/
       *---I           *---------------------*
*/

/** 
  Cell (octant) identified by its integer position in a map. It can also carry a data of type T.
*/

template <typename T> class OCell {
public:

  /** Constructor. */
  OCell() = default;

  /** Constructor. */
  OCell(const char plevel, const char psubIndex, const T &pdata);

  /** Copy contructor. */
  OCell(const OCell& copy);

  /** Assignment. */
  OCell& operator = (const OCell& copy);

  /** 
    OCell level, a background cell has level 0, all others should be [1..maxLevel - 1]  
    (see background). maxLevel is for NODES of smallest octants.
  */
  char level = 0;

  /** Sub-index within its parent, [0..7], see picture above. */
  char subIndex = 0;

  /** Can also carry a data */
  T data;
};


template <typename T> OCell<T>::OCell(const char plevel, const char psubIndex, const T &pdata)
{
  level = plevel;
  subIndex = psubIndex;
  data = pdata;
}

template <typename T> OCell<T>::OCell(const OCell& copy)
{
  level = copy.level;
  subIndex = copy.subIndex;
  data = copy.data;
}

template <typename T> OCell<T>& OCell<T>::operator = (const OCell<T>& copy)
{
  level = copy.level;
  subIndex = copy.subIndex;
  data = copy.data;
  return *this;
}

