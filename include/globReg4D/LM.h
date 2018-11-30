// ----------------------------------------------------------------------------
// -                       Fast Global Registration                           -
// ----------------------------------------------------------------------------
// The MIT License (MIT)
//
// Copyright (c) Intel Corporation 2016
// Qianyi Zhou <Qianyi.Zhou@gmail.com>
// Jaesik Park <syncle@gmail.com>
// Vladlen Koltun <vkoltun@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
// ----------------------------------------------------------------------------

#ifndef LM_H_
#define LM_H_

#include <vector>
#include <flann/flann.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <iostream>
#include <fstream>
#include <math.h>
#include <random>

//qa lib
#include "include/reg_common.h"
#include "../util/inputGen.h"

//#define DIV_FACTOR			1.1		// (1.4 by default) Division factor used for graduated non-convexity
//#define USE_ABSOLUTE_SCALE	1		// Measure distance in absolute scale (1) or in scale relative to the diameter of the model (0)
//#define MAX_CORR_DIST		0.05	// (0.025 by default) Maximum correspondence distance (also see comment of USE_ABSOLUTE_SCALE)
//#define ITERATION_NUMBER	512	// (64 by default) Maximum number of iteration
//#define TUPLE_SCALE			0.95	// (0.95 by default) Similarity measure used for tuples of feature points.
//#define TUPLE_MAX_CNT		1000	// (1000 by default) Maximum tuple numbers.


using namespace Eigen;
using namespace std;

typedef vector<Vector3f> Points;
typedef vector<VectorXf> Feature;

class CApp{
public:

    double divFactor;
    bool useAbsoluteScale;
    double maxCorrDist;
    int iterationNumber;
    double tupleScale;
    int tupleMaxCnt;

    CApp():divFactor(1.1),
        useAbsoluteScale(true),
        maxCorrDist(0.05),
        iterationNumber(512),
        tupleScale(0.95),
        tupleMaxCnt(1000){}
    
    void pclPC2CAppPoints(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, Points &out);


    void NormalizePoints();
    void AdvancedMatching();
    Matrix4f GetTrans();
    //compute the residual for current iteration
    double computeResidual(const Points &p, const Points &q, const double &mu);
    double OptimizePairwise(bool decrease_mu_, int numIter_);
    double OptimizePairwise_4Dof(bool decrease_mu_, int numIter_);

    void LM(pcl::PointCloud<pcl::PointXYZ>::Ptr cloudS, pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT, vector<corrTab> &corr, double inlTh, Transform3 &result);
private:
    // containers
    vector<Points> pointcloud_;
    vector<Feature> features_;
    Matrix4f TransOutput_;
    vector<pair<int, int>> corres_;

    // for normalization
    Points Means;
    float GlobalScale;
    float StartScale;

    void SearchFLANNTree(flann::Index<flann::L2<float>>* index,
                         VectorXf& input,
                         std::vector<int>& indices,
                         std::vector<float>& dists,
                         int nn);
};

#endif
