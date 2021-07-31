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

#include <assert.h>
#include <vector>
#include "Vector.h"

/**
  Simple solution to linear, over/underdetermined system with plain C, no
SIMD.
*/
                              // Solve square system of equations
template <class T> bool solveSystem(int N, T A[], T B[], T tolerance)
{
  int K,K1,J,I;
  int I0,I1,I2;
  T AKK;

  for (K = 0; K < N; K++)
  {
    K1 = K + 1;
    AKK = A[K*N+K];

    if (std::abs(AKK) < tolerance)
	  {	
	    return false;
    }

    B[K] /= AKK;
    if (K == (N-1)) break;

    for (J = K1; J < N; J++)
    {
      I0 = K*N+J;
      A[I0] /= AKK;
      for (I = K1; I < N; I++) 
      {
        I1 = I*N+J;
        I2 = I*N+K;
        A[I1] -= (A[I2] * A[I0]);
      }
      B[J] -= (A[J*N+K] * B[K]);
    }
  }

Back:
  K1 = K;
  K -= 1;
  if (K < 0) goto Success;
  for (J = K1; J < N; J++) B[K] -= A[K*N+J] * B[J];
  goto Back;

Success :
  return true;
}
                              // Reshuffle matrices A and B from row
template <class T, class Tint> 
void findPivot(Tint N, Tint row, T A[], T B[], std::vector<T> &temp)
{
  for (Tint K = row; K < N - 1; K++)
  {
    Tint K1 = K + 1;
                          
    // Find equation with the largest diagonal element
    Tint index = K; T dmax = std::abs(A[K*N+K]);
    for (Tint i = K1; i < N; i++)
    {
      Tint dindex = i * N + K;
      if (std::abs(A[dindex]) > dmax)
      {
        index = i;
        dmax = std::abs(A[dindex]);
      }
    }
                          
    // Make permutation
    if (index != K)
    {
      // Swap matrix raws
      memmove(&temp[0],&A[K * N],sizeof(T) * N);
      memmove(&A[K * N],&A[index * N],sizeof(T) * N);
      memmove(&A[index * N],&temp[0],sizeof(T) * N);

      // Swap right-hand side values
      T btemp = B[K];
      B[K] = B[index];
      B[index] = btemp;
    }
  }
}


                              // Solve linear system of equations Ax = B with 
                              // pivoting
template <class T, class Tint> 
bool 
solveSystemWithPivoting(Tint N, T A[], T B[], T tolerance)
{
  Tint K,K1,J,I;
  Tint I0,I1,I2;
  T AKK;

  // Temporary storage to swap matrix rows
  std::vector<T> temp(N);

  // temp
  std::vector<T> akk;
  T akkmin;

  // Solution itself
  for (K = 0; K < N; K++)
  {
//!!!
    //printf("%d of %d\r",K,N);

    // Reshuffle matrices A and B
    findPivot<T,Tint>(N,K,A,B,temp);

    K1 = K + 1;
    AKK = A[K*N+K];

    // Temp
    akk.push_back(AKK);
    if (K == 0)
    {
      akkmin = std::abs(AKK);
    } else
    {
      if (std::abs(AKK) < akkmin) 
        akkmin = std::abs(AKK);
    }

    // Pivot is too small
    if (std::abs(AKK) < tolerance)
    {	
      return false;
    }

    B[K] /= AKK;
    if (K == (N-1)) break;

    for (J = K1; J < N; J++)
    {
      I0 = K*N+J;
      A[I0] /= AKK;
      for (I = K1; I < N; I++) 
      {
        I1 = I*N+J;
        I2 = I*N+K;
        A[I1] -= (A[I2] * A[I0]);
      }
      B[J] -= (B[K] * A[J*N+K]);
    }
  }

  // Back substitution
  for (;;)
  {
    K1 = K;
    K -= 1;
    if (K < 0) 
      break;
    for (J = K1; J < N; J++) 
      B[K] -= B[J] * A[K*N+J];
  }

  return true;
}

 
//                              // Solve AT*A = B system, N1(number of rows) x 
//                              // N2(columns) of A; right-hand side is B - N2 
//                              // vectors of three components; solution vector 
//                              // is in B
//template <class T, class Tint> 
//bool 
//SolveATASystem(Tint N1, Tint N2, T A[], TVector<T> B[], T tolerance)
//{
//  // AT * A
//  std::vector<T> ATA(N2 * N2);
//
//  // To deallocate AT
//  {
//    // Make transpose
//    std::vector<T> AT(N1 * N2);
//
//    transpose(A,N1,N2,&AT[0]);
//
//    fastM(&AT[0],A,&ATA[0],N2,N1,N2,N1,N2,N2);
//  }
//                              
//  // Solve system
//  bool res = solveSystemWithPivoting(N2,&ATA[0],B,tolerance);
//
//  return res;
//}

