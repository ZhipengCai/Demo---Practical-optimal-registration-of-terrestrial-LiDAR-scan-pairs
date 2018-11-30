/////////////////////////////////////////////////////////////////////////////
//
//              R O T A T I O N   S E A R C H
//
// This package contains the source code which implements the
// BnB rotation search algorithm and the nested 6 DoF registration
// algorithm proposed in
//
// A. Parra Bustos, T.-J. Chin, A. Eriksson, H. Li and D. Suter
// Fast Rotation Search with Stereographic Projections for 3D Registration
// IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI)
//
// Copyright (c) 2016 Alvaro PARRA BUSTOS (aparra@cs.adelaide.edu.au.)
// School of Computer Science, The University of Adelaide, Australia
// The Australian Center for Visual Technologies
// http://cs.adelaide.edu.au/~aparra
// Please acknowledge the authors by citing the above paper in any academic
// publications that have made use of this package or part of it.
//
/////////////////////////////////////////////////////////////////////////////

#ifndef REG_DATA_INDEXATION_
#define REG_DATA_INDEXATION_

#include <vector>

namespace reg {
namespace search {

template<class SSR>
class DataIndexation
{
public:
    virtual ~DataIndexation(){}
    virtual int evalUpperBound(SSR ssr, int lwbnd) const = 0;
    virtual int evalUpperBound(SSR ssr, int lwbnd,
                               int *matchList, int matchListSize,
                               std::vector<bool> &matches) const = 0;

    virtual int evalLowerBound(SSR ssr) const = 0;
    virtual int evalLowerBound(SSR ssr, int *matchList, int matchListSize) const = 0;
    virtual int evalLowerBound(SSR ssr, int *matchList, int matchListSize,
                               std::vector<bool> &matches) const = 0;

    virtual int sweep() = 0;
    virtual int size() const = 0;

};


inline double dist_sq( Vector3 &a1, double *a2)
{
  double dist_sq = 0, diff;
  for (int i=0; i<3;i++)
  {
    diff = (a1[i] - a2[i]);
    dist_sq += diff*diff;
  }
  return dist_sq;
}


} // End namespace sarch
} // End namespace reg

#endif

