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

#ifndef REG_NEAREST_NEIGHBOUR_SEARCH_
#define REG_NEAREST_NEIGHBOUR_SEARCH_

#include "reg_common.h"
#include "reg_dtransf.h"
#include <pcl/kdtree/kdtree_flann.h>


namespace reg
{
namespace search
{

typedef pcl::PointXYZ Point;
typedef pcl::PointCloud<Point>::Ptr Cloud_ptr;


class NearestNeighbourSearch
{
public:
    //virtual double distance(Point p)=0;
    virtual double distance(const Eigen::Vector3f &p)  =0 ;
    virtual double distance(const Point &p)  =0;
    virtual ~NearestNeighbourSearch(){}
};


class KDTreeNearestNeighbourSearch: public NearestNeighbourSearch
{
private:
    std::vector<int> idx;
    std::vector<float> sqrdist;
    pcl::KdTreeFLANN<Point>::Ptr kdtree_ptr;

public:
    KDTreeNearestNeighbourSearch(const Cloud_ptr &B);
    ~KDTreeNearestNeighbourSearch();

    double distance(const Eigen::Vector3f &p) ;
    double distance(const Point &p) ;
};


class DTransfNearestNeighbourSearch: public NearestNeighbourSearch
{
private:
    reg::DTransf dtransf;

public:
    DTransfNearestNeighbourSearch(const Cloud_ptr B,
                                  int resolution, double borderRate=0);
    ~DTransfNearestNeighbourSearch();

    double distance(const Point &p) ;
    double distance(const Eigen::Vector3f &p) ;
};



} // End namespace geometry
} // End namespace reg
#endif
