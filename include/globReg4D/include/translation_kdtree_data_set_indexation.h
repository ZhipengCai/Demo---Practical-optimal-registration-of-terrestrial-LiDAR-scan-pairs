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


#ifndef REG_TRANSLATION_KDT_DATA_INDEXATION_
#define REG_TRANSLATION_KDT_DATA_INDEXATION_

#include "reg_common.h"
#include "data_indexation.h"
#include "state.h"
#include "reg_rtree.h"
#include <nanoflann.hpp>


using namespace nanoflann;

namespace reg {
namespace search {
   

template<class SSR>
class TranslationKDTreeDataSetIndexation : public DataIndexation<SSR>
{
public:

    TranslationKDTreeDataSetIndexation(const Matrix3X &M, const Matrix3X &B, double th);

    ~TranslationKDTreeDataSetIndexation();

    int sweep();
    int size() const;

    int evalUpperBound(SSR ssr, int lwbnd) const;
    int evalUpperBound(SSR ssr, int lwbnd,
                       int *matchList, int matchListSize,
                       std::vector<bool> &matches) const;

    int evalLowerBound(SSR ssr) const;
    int evalLowerBound(SSR ssr, int *matchList, int matchListSize) const;
    int evalLowerBound(SSR ssr, int *matchList, int matchListSize,
                       std::vector<bool> &matches) const;

    void setM( Matrix4X *M);

private:

    const Matrix3X &M_in;
    const Matrix3X &B_in;
    const double th;

    Matrix4X *M;
    //pcl::KdTreeFLANN<Point> treeB;
    //flann::Index<flann::L2<double> > treeB;
    
    typedef KDTreeSingleIndexAdaptor<
    L2_Simple_Adaptor<double, PointCloud > ,
    PointCloud, 3  > KdTree;

    
    PointCloud *cloud;
    KdTree *treeB;
    
    int _size;
};

template<class SSR>
class TranslationKDTreeDataSetIndexation4Dof : public DataIndexation<SSR>
{
public:

    TranslationKDTreeDataSetIndexation4Dof(const Matrix3X &M, const Matrix3X &B, const Vector &TH);

    ~TranslationKDTreeDataSetIndexation4Dof();

    int sweep();
    int size() const;

    int evalUpperBound(SSR ssr, int lwbnd) const;
    int evalUpperBound(SSR ssr, int lwbnd,
                       int *matchList, int matchListSize,
                       std::vector<bool> &matches) const;

    int evalLowerBound(SSR ssr) const;
    int evalLowerBound(SSR ssr, int *matchList, int matchListSize) const;
    int evalLowerBound(SSR ssr, int *matchList, int matchListSize,
                       std::vector<bool> &matches) const;

    void setM( Matrix4X *M);

private:

    const Matrix3X &M_in;
    const Matrix3X &B_in;
    const Vector &th;

    Matrix4X *M;
    //pcl::KdTreeFLANN<Point> treeB;
    //flann::Index<flann::L2<double> > treeB;

    typedef KDTreeSingleIndexAdaptor<
    L2_Simple_Adaptor<double, PointCloud > ,
    PointCloud, 3  > KdTree;


    PointCloud *cloud;
    KdTree *treeB;

    int _size;
};

} // End namespace sarch
} // End namespace reg

#include "translation_kdtree_data_set_indexation.hpp"

#endif
