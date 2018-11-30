#include "../../include/globReg4D/globReg4D.h"

void GLOBREG::globReg4D(pcl::PointCloud<pcl::PointXYZ>::Ptr cloudS, pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT,
                        vector<corrTab> &corr, vector<int> &corrNOS, vector<int> &corrNOT,
                        double inlTh, double tCubeSize, Transform3 &result, bool useFPA){
    //convert data format
    Matrix3X matrixS;
    DataTrans(cloudS, matrixS);
    Matrix3X matrixT;
    DataTrans(cloudT, matrixT);

    if (useFPA){
        cout<<"NO correspondences before FMP = "<<corr.size()<<endl;
        FPA(matrixS,matrixT,corr,corrNOS,corrNOT,inlTh,result);
        cout<<"NO correspondences after FMP = "<<corr.size()<<endl;
    }
    double tDom[6] = {-tCubeSize,-tCubeSize,-tCubeSize,tCubeSize,tCubeSize,tCubeSize};
    TranslationSearchSpaceRegion3DOF trSSR(tDom);
    //    int q = 0;
    int currLwb = nestedbnb_search_4dof_withCorr_IS(matrixS, matrixT, corr, corrNOS, corrNOT,
                                                    result, inlTh, 0,
                                                    0,0,trSSR,false);

    cout<<"final consensus = "<<currLwb<<endl;
}

int GLOBREG::FPA(const Matrix3X &matrixS, const Matrix3X &matrixT,
         vector<corrTab> &corr, vector<int> &corrNOS, vector<int> &corrNOT,
         double inlTh, Transform3 &result){
    int currLwb = 2;
    //1 compute the upper bound
    AxisAngle rsearchResult;

    Matrix3X matrixSTmp(matrixS.cols());
    Matrix3X matrixTTmp(matrixT.cols());
    for(size_t i=0;i<corr.size();i++){
        //translate the idx-th pair of corresponding points to the origin
        Vector3 CentS(-matrixS(0,corr[i].idxS),-matrixS(1,corr[i].idxS),-matrixS(2,corr[i].idxS));
        Vector3 CentT(-matrixT(0,corr[i].idxT),-matrixT(1,corr[i].idxT),-matrixT(2,corr[i].idxT));
        transPC(matrixSTmp, matrixS, CentS);
        transPC(matrixTTmp, matrixT, CentT);
        //evaluate upper bound
        corr[i].upBnd = rot1_withCorr_IntervalStab(matrixSTmp,matrixTTmp,corr,corrNOS,corrNOT,
                                      2*inlTh,0,currLwb-1,0,rsearchResult);

        if(corr[i].upBnd < currLwb){
            corrNOS[corr[i].idxS]--;
            corrNOT[corr[i].idxT]--;
            corr.erase(corr.begin()+i);
            i--;
        }
        else if(corr[i].upBnd > currLwb){
            Matrix3 R;
            fromAxisAngle(R,rsearchResult);
            //rotate matrixSTmp
            rotate1D_aroundZ(matrixSTmp,rsearchResult.w);
            //evaluate lower bound
            corr[i].lwBnd = evalObj(matrixSTmp,matrixTTmp,corr,inlTh);
            if(corr[i].lwBnd>currLwb){
                currLwb = corr[i].lwBnd;
                Transform3 guessAndResultTmp(R);
                result = TransMatCompute(guessAndResultTmp,CentS,CentT);
            }
        }
    }
    //remove outliers again using current best lwbnd
    for(size_t i=0;i<corr.size();i++){
        if(corr[i].upBnd < currLwb){
            corrNOS[corr[i].idxS]--;
            corrNOT[corr[i].idxT]--;
            corr.erase(corr.begin()+i);
            i--;
        }
    }
    return currLwb;
}


int GLOBREG::evalObj(Matrix3X matrixS, Matrix3X matrixT,
                       std::vector<corrTab> corr, double inlTh, bool enforceOnetoOne)
{
    int obj = 0;
    std::vector<bool> hasMatch(matrixS.n,false);
    double sqTh = inlTh*inlTh;
    for(int i=0;i<corr.size();i++){
        int idxS = corr.data()[i].idxS;
        int idxT = corr.data()[i].idxT;
        //evaluate only when current source point still has no match
        if(enforceOnetoOne){
            if(!hasMatch[idxS]){
                double disX = matrixS(0,idxS)-matrixT(0,idxT);
                double disY = matrixS(1,idxS)-matrixT(1,idxT);
                double disZ = matrixS(2,idxS)-matrixT(2,idxT);
                if(disX*disX+disY*disY+disZ*disZ<=sqTh){
                    hasMatch[idxS] = true;
                    obj++;
                }
            }
        }
        else{
            double disX = matrixS(0,idxS)-matrixT(0,idxT);
            double disY = matrixS(1,idxS)-matrixT(1,idxT);
            double disZ = matrixS(2,idxS)-matrixT(2,idxT);
            if(disX*disX+disY*disY+disZ*disZ<=sqTh){
                obj++;
            }
        }

    }
    return obj;
}


int GLOBREG::RANSAC(pcl::PointCloud<pcl::PointXYZ>::Ptr cloudS, pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT,
           vector<corrTab> &corr, double inlTh, Transform3 &result){
    double conf = 0.99;
    Matrix3X matrixS;
    DataTrans(cloudS, matrixS);
    Matrix3X matrixT;
    DataTrans(cloudT, matrixT);

    /* initialize random seed: */
    std::srand (std::time(NULL));
    int iter = 0;
    int maxIter;
    if (corr.size()<sqrt(1e7)){
        maxIter = corr.size()*(corr.size()-1)/2;
    }
    else{
        maxIter = 1e7;
    }


    //generate random series
    std::vector<std::vector<int>> pool(corr.size());
    for(size_t i=0;i<pool.size();i++){
        std::vector<int> tmp(corr.size()-i-1);
        std::iota( std::begin( tmp ), std::end( tmp ), i+1 );
        pool[i] = tmp;
    }
    //need only two point samples for 4DOF
    int idx[2];

    int objBest = 0;
    int skipNO = 0;
    while (iter<maxIter){
        iter++;

        idx[0] = rand() % (corr.size()-1);
        int repeatSample = 0;
        while(pool.data()[idx[0]].size()==0){
            idx[0] = rand() % (corr.size()-1);
            if(repeatSample>corr.size()){
                cout<<"not able to find a not repeated sample in "<<repeatSample<<" times"<<endl;
                return objBest;
            }
        }
        int tmp = rand() % pool.data()[idx[0]].size();
        idx[1] = pool.data()[idx[0]].data()[tmp];
        pool.data()[idx[0]].erase(pool.data()[idx[0]].begin()+tmp);
        //compute rigid transformation
        std::vector<Eigen::Vector3d> pcS(2), pcT(2);
        for(size_t i=0;i<2;i++){
            //source sample points
            pcS.data()[i](0) = cloudS->points[corr.data()[idx[i]].idxS].x;
            pcS.data()[i](1) = cloudS->points[corr.data()[idx[i]].idxS].y;
            pcS.data()[i](2) = cloudS->points[corr.data()[idx[i]].idxS].z;
            //target sample points
            pcT.data()[i](0) = cloudT->points[corr.data()[idx[i]].idxT].x;
            pcT.data()[i](1) = cloudT->points[corr.data()[idx[i]].idxT].y;
            pcT.data()[i](2) = cloudT->points[corr.data()[idx[i]].idxT].z;
        }

        //don't use sample points whose location is too close on the xy plane
        double diffxS = pcS.data()[0](0)-pcS.data()[1](0);
        double diffyS = pcS.data()[0](1)-pcS.data()[1](1);
        if(fabs(diffxS)<=1e-6&&fabs(diffyS)<=1e-6){
            skipNO++;
            continue;
        }
        double diffxT = pcT.data()[0](0)-pcT.data()[1](0);
        double diffyT = pcT.data()[0](1)-pcT.data()[1](1);
        if(fabs(diffxT)<=1e-6&&fabs(diffyT)<=1e-6){
            skipNO++;
            continue;
        }

        //compute rigid transformation
        Eigen::Vector3d aveS = (pcS.data()[0]+pcS.data()[1])/2;
        Eigen::Vector3d aveT = (pcT.data()[0]+pcT.data()[1])/2;

        //find 1d rotation
        double a = 0, b = 0;
        for(size_t i=0;i<2;i++){
            Eigen::Vector3d sTmp,tTmp;
            sTmp = pcS.data()[i]-aveS;
            tTmp = pcT.data()[i]-aveT;
            a += (sTmp(0)*tTmp(0)+sTmp(1)*tTmp(1));
            b += (sTmp(0)*tTmp(1)-sTmp(1)*tTmp(0));
        }

        double theta;
        if (fabs(a)<1e-9){
            if(a>=0){
                a=1e-9;
            }
            else{
                a=-1e-9;
            }
        }
        theta = atan(b/a);
        if(a<0){
            theta = theta+M_PI;
        }

        //compute optimal translation
        Vector3 optt(0,0,0);
        double cosTheta = cos(theta);
        double sinTheta = sin(theta);
        optt.x = aveT(0)-(cosTheta*aveS(0)-sinTheta*aveS(1));
        optt.y = aveT(1)-(sinTheta*aveS(0)+cosTheta*aveS(1));
        optt.z = aveT(2)-aveS(2);

        //transform point clouds
        Matrix3X matrixSTmp(matrixS);
        //rotate
        rotate1D_aroundZ(matrixSTmp, theta);
        //translate
        transPC(matrixSTmp, matrixSTmp,optt);
        //evaluate obj
        int obj = evalObj(matrixSTmp,matrixT,corr,inlTh,false);
        if (obj>objBest){
            //update solution
            objBest = obj;
            cout<<"updating solution, obj = "<<objBest<<endl;
            result = Transform3(optt);
            result.x[0] = cosTheta;
            result.x[4] = -sinTheta;
            result.x[1] = sinTheta;
            result.x[5] = cosTheta;
            //update max number of iterations
            double maxIterTmp = std::ceil(std::log(1-conf)/std::log(1-std::pow((objBest/(float)corr.size()),2)));
            if(maxIterTmp>0&&maxIterTmp<INT_MAX)
                maxIter = min(maxIter, (int)maxIterTmp);
            cout<<"maxIter = "<<maxIter<<endl;
        }
    }
    return objBest;
}
