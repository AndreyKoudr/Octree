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
  A node for a finite-element linear conforming approximation within octree. Can also carry data of type T.

  When you add a cell to an octree, 8 new nodes are inserted into the map, increasing their reference counts.
  When a cell is excluded, its nodal ref counts are decreased and nodes are physically removed from total node
  map is only ref count reaches zero.
*/

template <typename T> class ONode {
public:

  ONode() = default;

  /** Constructor. There should be no default constructor as we need to specify data
    associated with the node. */
  ONode(const T& pdata);

  /** Copy contructor */
  ONode(const ONode& copy);

  /** Assignment operator. */
  ONode& operator = (const ONode& copy);

  /** Reference count to safely exclude from list. */
  LINT refCount = 0;

  /** Data to hold within this node. */
  T data;
};

template <typename T> ONode<T>::ONode(const T& pdata)
{
  data = pdata;
}

template <typename T> ONode<T>::ONode(const ONode<T>& copy)
{
  refCount = copy.refCount;
  data = copy.data;
}

template <typename T> ONode<T>& ONode<T>::operator = (const ONode<T>& copy)
{
  refCount = copy.refCount;
  data = copy.data;
  return *this;
}

