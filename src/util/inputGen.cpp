#include "../../include/util/inputGen.h"

void GenIn::VGF(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, pcl::PointCloud<pcl::PointXYZ>::Ptr cloudVG,
                double inlTh){
    //format for filtering
    pcl::PCLPointCloud2::Ptr cloud2 (new pcl::PCLPointCloud2 ());
    pcl::PCLPointCloud2::Ptr cloudVG2 (new pcl::PCLPointCloud2 ());
    pcl::toPCLPointCloud2(*cloud,*cloud2);
    //set up filtering parameters
    pcl::VoxelGrid<pcl::PCLPointCloud2> sor;
    sor.setInputCloud (cloud2);
    sor.setLeafSize (inlTh, inlTh, inlTh);
    //filtering process
    sor.filter (*cloudVG2);
    pcl::fromPCLPointCloud2(*cloudVG2,*cloudVG);
}

void GenIn::ISSExt(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, pcl::PointCloud<pcl::PointXYZ>::Ptr ISS,
            pcl::PointIndicesPtr ISSIdx, double inlTh){
    double iss_salient_radius_ = 6 * inlTh;
    double iss_non_max_radius_ = 4 * inlTh;
    double iss_gamma_21_ (0.975);
    double iss_gamma_32_ (0.975);
    double iss_min_neighbors_ (5);
    int iss_threads_ (8); //switch to the number of threads in your cpu for acceleration

    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ> ());
    pcl::ISSKeypoint3D<pcl::PointXYZ, pcl::PointXYZ> iss_detector;

    iss_detector.setSearchMethod (tree);
    iss_detector.setSalientRadius (iss_salient_radius_);
    iss_detector.setNonMaxRadius (iss_non_max_radius_);
    iss_detector.setThreshold21 (iss_gamma_21_);
    iss_detector.setThreshold32 (iss_gamma_32_);
    iss_detector.setMinNeighbors (iss_min_neighbors_);
    iss_detector.setNumberOfThreads (iss_threads_);
    iss_detector.setInputCloud (cloud);
    iss_detector.compute (*ISS);
    ISSIdx->indices = iss_detector.getKeypointsIndices()->indices;
    ISSIdx->header = iss_detector.getKeypointsIndices()->header;
}

void GenIn::FPFHComp(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, double inlTh, pcl::PointIndicesPtr ISSIdx, pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhOut){
    //compute normal
    pcl::PointCloud<pcl::Normal>::Ptr normal (new pcl::PointCloud<pcl::Normal> ());
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ> ());
    pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> neSource;
    neSource.setInputCloud (cloud);
    neSource.setSearchMethod (tree);
    neSource.setRadiusSearch (3*inlTh);
    neSource.compute (*normal);

    //compute fpfh using normals
    pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfhEst;
    fpfhEst.setInputCloud(cloud);
    fpfhEst.setInputNormals(normal);
    fpfhEst.setSearchMethod(tree);
    fpfhEst.setRadiusSearch(8*inlTh);
    fpfhEst.setNumberOfThreads (4);
    fpfhEst.setIndices(ISSIdx);
    fpfhEst.compute(*fpfhOut);
}

void GenIn::corrComp(pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs,
                                   pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfht,
                                   vector<corrTab> &corr, int MaxNOCorrPP,
                                   vector<int> &corrNOS, vector<int> &corrNOT){
    int n = min(MaxNOCorrPP, (int)fpfht->size()); //maximum number of correspondences to find for each source point
    corr.clear();
    corrNOS.assign(fpfhs->size(),0);
    corrNOT.assign(fpfht->size(),0);
    // Use a KdTree to search for the nearest matches in feature space
    pcl::KdTreeFLANN<pcl::FPFHSignature33> treeS;
    treeS.setInputCloud (fpfhs);
    pcl::KdTreeFLANN<pcl::FPFHSignature33> treeT;
    treeT.setInputCloud (fpfht);
    for(size_t i=0;i<fpfhs->size();i++){
        vector<int> corrIdxTmp(n);
        vector<float> corrDisTmp(n);
        //find the best n matches in target fpfh
        treeT.nearestKSearch(*fpfhs,i,n,corrIdxTmp,corrDisTmp);
        for(size_t j=0;j<corrIdxTmp.size();j++){
            bool removeFlag = true;
            int searchIdx = corrIdxTmp[j];
            vector<int> corrIdxTmpT(n);
            vector<float> corrDisTmpT(n);
            treeS.nearestKSearch(*fpfht,searchIdx,n,corrIdxTmpT,corrDisTmpT);
            for(size_t k=0;k<n;k++){
                if(corrIdxTmpT.data()[k]==i){
                    removeFlag = false;
                    break;
                }
            }
            if(removeFlag == false){
                corrTab corrTabTmp;
                corrTabTmp.idxS = i;
                corrTabTmp.idxT = corrIdxTmp[j];
                corrTabTmp.disFeature = corrDisTmp[j];
                corrTabTmp.upBnd = fpfht->size();
                corrTabTmp.lwBnd = 1;
                corr.push_back(corrTabTmp);
                corrNOS[i]++;
                corrNOT[corrIdxTmp[j]]++;
            }
        }
    }
}

void GenIn::CentAndRComp(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, Vector3 &center, double &radius){
    pcl::PointXYZ xyzMin, xyzMax;
    //reset the bounding box
    xyzMax.x=xyzMax.y=xyzMax.z=-FLT_MAX;
    xyzMin.x=xyzMin.y=xyzMin.z=FLT_MAX;
    for(size_t i=0;i<cloud->size();i++){
        xyzMax.x=max(xyzMax.x,cloud->points[i].x);
        xyzMax.y=max(xyzMax.y,cloud->points[i].y);
        xyzMax.z=max(xyzMax.z,cloud->points[i].z);
        xyzMin.x=min(xyzMin.x,cloud->points[i].x);
        xyzMin.y=min(xyzMin.y,cloud->points[i].y);
        xyzMin.z=min(xyzMin.z,cloud->points[i].z);
    }
    center.x = (xyzMax.x+xyzMin.x)*0.5;
    center.y = (xyzMax.y+xyzMin.y)*0.5;
    center.z = (xyzMax.z+xyzMin.z)*0.5;

    double radiusX = (xyzMax.x-center.x);
    double radiusY = (xyzMax.y-center.y);
    double radiusZ = (xyzMax.z-center.z);
    radius = sqrt(radiusX*radiusX+radiusY*radiusY+radiusZ*radiusZ);
}

void GenIn::transPC(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, Vector3 trans){
    for(size_t i=0;i<cloud->size();i++){
        cloud->points[i].x += trans.x;
        cloud->points[i].y += trans.y;
        cloud->points[i].z += trans.z;
    }
}

void GenIn::transPC(Matrix3X &cloud, const Matrix3X &cloudIn, Vector3 trans){
    for (size_t i=0;i<cloud.cols(); i++)
    {
        cloud.x[i*3]   = cloudIn.x[i*3] + trans.x;
        cloud.x[i*3+1] = cloudIn.x[i*3+1] + trans.y;
        cloud.x[i*3+2] = cloudIn.x[i*3+2] + trans.z;
    }
}

void GenIn::rotate1D_aroundZ(Matrix3X &x, double w){
    double cosw = cos(w);
    double sinw = sin(w);
    for(size_t i=0;i<x.cols();i++){
        double xx = x(0,i);
        double xy = x(1,i);
        x(0,i) = xx*cosw-xy*sinw;
        x(1,i) = xx*sinw+xy*cosw;
    }
}

void GenIn::DataTrans(pcl::PointCloud<pcl::PointXYZ>::Ptr cloudIn, Matrix3X &cloudOut){
    cloudOut.setSize(cloudIn->size());
    for(size_t i=0;i<cloudIn->size();i++){
        cloudOut.x[3*i] = (double) cloudIn->points[i].x;
        cloudOut.x[3*i+1] = (double) cloudIn->points[i].y;
        cloudOut.x[3*i+2] = (double) cloudIn->points[i].z;
    }
}



Transform3 GenIn::TransMatCompute(const Transform3 &T, const Vector3 &vS, const Vector3 &vT){
    Matrix3 R;
    R.m[0] = T.x[0];
    R.m[1] = T.x[1];
    R.m[2] = T.x[2];

    R.m[3] = T.x[4];
    R.m[4] = T.x[5];
    R.m[5] = T.x[6];

    R.m[6] = T.x[8];
    R.m[7] = T.x[9];
    R.m[8] = T.x[10];

    Vector3 t(T.x[12], T.x[13], T.x[14]);
    t = t-vT+multiply(R,vS.m);

    return Transform3(R,t);
}


