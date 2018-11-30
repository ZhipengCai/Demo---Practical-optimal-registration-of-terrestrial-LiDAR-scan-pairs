#ifndef INPUTGEN_H_
#define INPUTGEN_H_

#include <iostream>
#include <vector>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/search/pcl_search.h>
#include <pcl/keypoints/iss_3d.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/conversions.h>
#include <pcl/features/fpfh_omp.h>
#include <pcl/features/normal_3d.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "../globReg4D/include/reg_common.h"

using namespace std;
using namespace reg;

namespace GenIn{

//struct corrTab{ //correspondence table, including the indices in corrIdx and the upperbound
//    int idxS; //first idx in corrIdx
//    int idxT; //second idx in corrIdx
//    int upBnd; //upperbound for this correspondences
//    int lwBnd;
//    double disFeature; //distance between the feature descriptors
//};

//voxel grid filter
void VGF(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, pcl::PointCloud<pcl::PointXYZ>::Ptr cloudVG,
	 double inlTh);
//ISS keypoint extraction
void ISSExt(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, pcl::PointCloud<pcl::PointXYZ>::Ptr ISS, 
            pcl::PointIndicesPtr ISSIdx, double inlTh);
//FPFH computation
void FPFHComp(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, double inlTh, pcl::PointIndicesPtr ISSIdx, pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhOut);
//correspondence computation
void corrComp(pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs,
              pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfht,
              vector<corrTab> &corr,
              int MaxNOCorrPP, vector<int> &corrNOS, vector<int> &corrNOT);
//compute the center and radius of a point cloud
void CentAndRComp(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, Vector3 &center, double &radius);
//translate a pcl point cloud
void transPC(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, Vector3 trans);
void transPC(Matrix3X &cloud, const Matrix3X &cloudIn, Vector3 trans);
//1D rotation
void rotate1D_aroundZ(Matrix3X &x, double w);
//transform the pcl point cloud to Matrix3X
void DataTrans(pcl::PointCloud<pcl::PointXYZ>::Ptr cloudIn, Matrix3X &cloudOut);
//adding the centering operation into the final rigid transformation
Transform3 TransMatCompute(const Transform3 &T, const Vector3 &vS, const Vector3 &vT);
}

#endif
