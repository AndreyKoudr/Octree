#ifndef PLANE_HPP
#define PLANE_HPP

#include "Vector.h"
#include "Systems.h"

#include <limits>
#include <assert.h>

// Build a plane defined by normal N (normalised) and scalar D. Vectors 01 and 02  should 
// not be collinear
template <class T> bool makePlaneOf3Vectors(const TVector<T> &V0, 
  const TVector<T> &V1, const TVector<T> &V2,
  TVector<T> &N, T &D, const T tolerance)
{
  N = (V1 - V0) ^ (V2 - V0);
  T len = !N;
  if (len > tolerance)
  {
    // normal normalised
    N = N * (static_cast<T>(1.0) / len); D = -(N * V0); 
    return true;
  }
  else
  {
    return false;
  }
}

// Project vector on plane
template <class T> void ProjectVectorOnPlane(TVector<T> V, TVector<T> N, T D, TVector<T> *Projection)
{
  T lambda = (-D - V * N) / (N * N);
  *Projection = V + N * lambda;
};

template <typename T>
struct Plane {

  typedef TVector<T> vector3;

  // Plane normal (normalised) - normalised together with Plane constant
  vector3 normal = vector3(1.0,0.0,0.0);
  T constant = 0.0;

  // constructors
  Plane() = default;

  // copy constructor
  Plane(const Plane& copy)
  {
    normal = copy.normal;
    constant = copy.constant;
  }

  // assignment operator
  Plane& operator = (const Plane& copy)
  {
    normal = copy.normal;
    constant = copy.constant;
    return *this;
  }

  // constructor from 4 scalars
  Plane(const T A, const T B, const T C, const T D)
  {
    normal = vector3(A,B,C);
    constant = D;
  }

  // constructor from normal and scalar
  Plane(const vector3 &N, const T D)
  {
    normal = N;
    constant = D;
  }

  // constructor form normal and point on Plane
  Plane(const vector3 &N, const vector3 &P)
  {
    normal = N;
    constant = - (P * N);
  }

  // constructor from 3 points
  Plane(const vector3 &V1, const vector3 &V2, const vector3 &V3, bool& OK)
  {
    OK = makePlaneOf3Vectors<T>(V1,V2,V3,normal,constant,V1.tolerance());
  }

  // constructor from more than 3 points 
  // (https://stackoverflow.com/questions/1400213/3d-least-squares-plane)
  Plane(const std::vector<vector3> &points, bool &ok)
  {
    assert(points.size() >= 3);

    T A[9] = {0};
    T B[3] = {0};

    for (size_t i = 0; i < points.size(); i++)
    {
      A[0 * 3 + 0] += points[i].X * points[i].X; 
      A[0 * 3 + 1] += points[i].X * points[i].Y; 
      A[0 * 3 + 2] += points[i].X; 

      A[1 * 3 + 0] += points[i].X * points[i].Y; 
      A[1 * 3 + 1] += points[i].Y * points[i].Y; 
      A[1 * 3 + 2] += points[i].Y; 

      A[2 * 3 + 0] += points[i].X; 
      A[2 * 3 + 1] += points[i].Y; 

      B[0] += points[i].X * points[i].Z;
      B[1] += points[i].Y * points[i].Z;
      B[2] += points[i].Z;
    }

    A[2 * 3 + 2] = T(points.size()); 

    T tolerance = std::numeric_limits<T>::epsilon() * static_cast<OREAL>(10.0);

    ok = solveSystemWithPivoting<T>(3,&A[0],&B[0],tolerance);

    if (ok)
    {
      vector3 n(B[0],B[1],T(-1.0));
      normal = +n;
      constant = B[2] / (!n);
    }
  }


  // tolerant equality
  bool operator == (Plane<T> &other) const
  {
    return (normal == other.normal && std::abs(constant - other.constant) < normal.tolerance());
  }

  // get distance (signed) from point to Plane
  T distance(const vector3 &V) const
  {
    return ((V * normal) + constant) / (normal * normal);
  }

  // get distance (signed) to another Plane
  T distance(const Plane& other) const
  {
    // normals must be normalised
    assert((!normal - static_cast<T>(1.0)) < normal.tolerance());
    assert((!other.normal - static_cast<T>(1.0)) < normal.tolerance());

    return (other.constant - constant);
  }

  // parallel to another Plane?
  bool is_parallel(const Plane& other) const
  {
    // normals must be normalised
    assert((!normal - static_cast<T>(1.0)) < normal.tolerance());
    assert((!other.normal - static_cast<T>(1.0)) < normal.tolerance());

    T angular_tolerance = std::numeric_limits<T>::epsilon() * static_cast<T>(100.0);

    return (!(normal ^ other.normal) < angular_tolerance);
  }

  // make projection of vector to plane
  TVector<T> makeProjection(const TVector<T> V) const
  {
    TVector<T> P;
    ProjectVectorOnPlane(V,normal,constant,&P);
    return P;
  };

};
 

#endif
