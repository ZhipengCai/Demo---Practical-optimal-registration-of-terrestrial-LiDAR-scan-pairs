#ifndef GLOBREG4D_H_
#define GLOBREG4D_H_

#include "../util/inputGen.h"
#include "include/reg_common.h"
#include "include/registration.h"

using namespace GenIn;
using namespace reg::search;

namespace GLOBREG{
//4DOF optimization
void globReg4D(pcl::PointCloud<pcl::PointXYZ>::Ptr cloudS, pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT,
               vector<corrTab> &corr, vector<int> &corrNOS, vector<int> &corrNOT,
               double inlTh, double tCubeSize, Transform3 &result, bool useFPA = true);

int FPA(const Matrix3X &matrixS, const Matrix3X &matrixT,
         vector<corrTab> &corr, vector<int> &corrNOS, vector<int> &corrNOT,
         double inlTh, Transform3 &result);

int evalObj(Matrix3X matrixS, Matrix3X matrixT, std::vector<corrTab> corr,
            double inlTh, bool enforceOnetoOne = false);

//RANSAC
int RANSAC(pcl::PointCloud<pcl::PointXYZ>::Ptr cloudS, pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT,
           vector<corrTab> &corr, double inlTh, Transform3 &result);
}

#endif
