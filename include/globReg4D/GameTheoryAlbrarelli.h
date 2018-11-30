#ifndef GAMETHEORYALBRARELLI_H
#define GAMETHEORYALBRARELLI_H

#include <random>

//qa lib
#include "../globReg4D/include/reg_common.h"

//pcl lib
#include <pcl/registration/transformation_estimation_svd.h>

//std lib
#include <math.h>
#include <algorithm> //min max


using namespace std;
using namespace reg;

class GameTheoryAlbrarelli
{
public:
    GameTheoryAlbrarelli();

    void GTReg(pcl::PointCloud<pcl::PointXYZ>::Ptr cloudS, pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT, vector<corrTab> &corr,
               double inlTh, reg::Transform3 &result);

};

#endif // GAMETHEORYALBRARELLI_H
