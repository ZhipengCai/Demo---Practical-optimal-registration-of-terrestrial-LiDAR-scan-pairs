#include <iostream>
#include <time.h>
#include <fstream>
#include <string>
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand

//external libraries
#include <pcl/io/auto_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/registration/icp.h>
#include <pcl/common/transforms.h>
//source files
#include "../include/util/inputGen.h"
#include "../include/globReg4D/globReg4D.h"
#include "../include/globReg4D/LM.h"
#include "../include/globReg4D/GameTheoryAlbrarelli.h"
#include <pcl/registration/ia_kfpcs.h>
#include <chrono>
//!uncomment if you have installed and want to use S4PCS (please reset to the directory of the installed library in your computer)
//#include "../S4PCS/Super4PCS/build/install/include/pcl/registration/super4pcs.h"

using namespace std;
using namespace GenIn;
using namespace GLOBREG;

void saveResult(string fnameOut, Transform3 result){
    //create file
    ofstream myfile;
    myfile.open(fnameOut);
    for(size_t i=0;i<4;i++){
        for(size_t j=0;j<4;j++){
            myfile<<result.x[i+j*4]<<" ";
        }
        myfile<<endl;
    }
    myfile.close();
}

void optimization(pcl::PointCloud<pcl::PointXYZ>::Ptr cloudS, pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT,
                  string fnameOut, double inlTh, int maxCorr, int testMethod){
    if(inlTh<=0||maxCorr<=0){
        cout<<"inlier threshold or maximum correspondence number must > 0"<<endl;
        return;
    }
    if(testMethod<1||testMethod>7){
        cout<< "wrong testing method, use numbers from 1 to 7 to specify"<<endl;
        return;
    }


    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> colorS(cloudS, 255, 0, 0);
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> colorT(cloudT, 255, 255, 0);

    //voxel grid filter (VGF)
    cout<<"performing voxel grid sampling with grid size = "<<inlTh<<endl;
    VGF(cloudS, cloudS, inlTh);
    VGF(cloudT, cloudT, inlTh);

    //show point clouds after VGF
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewerVGF (new pcl::visualization::PCLVisualizer ("After Voxel Grid Downsampling"));
    viewerVGF->setBackgroundColor (0, 0, 0);
    viewerVGF->addPointCloud<pcl::PointXYZ> (cloudS, colorS, "source cloud");
    viewerVGF->addPointCloud<pcl::PointXYZ> (cloudT, colorT, "target cloud");
    viewerVGF->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "After Voxel Grid Downsampling");
    viewerVGF->spinOnce ();

    //extract ISS
    pcl::PointCloud<pcl::PointXYZ>::Ptr issS (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointIndicesPtr issIdxS(new pcl::PointIndices);
    pcl::PointCloud<pcl::PointXYZ>::Ptr issT (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointIndicesPtr issIdxT(new pcl::PointIndices);
    cout<<"extracting ISS keypoints..."<<inlTh<<endl;
    ISSExt(cloudS,issS,issIdxS,inlTh);
    ISSExt(cloudT,issT,issIdxT,inlTh);
    cout<<"size of issS = "<<issS->size()<<"; size of issT = "<<issT->size()<<endl;

    //translating the center of both point clouds to the origin
    Vector3 centS(0,0,0), centT(0,0,0);
    double rS, rT;
    CentAndRComp(issS,centS,rS);
    CentAndRComp(issT,centT,rT);
    Vector3 mCentS(-centS.x,-centS.y,-centS.z);
    Vector3 mCentT(-centT.x,-centT.y,-centT.z);
    transPC(issS,mCentS);
    transPC(issT,mCentT);

    //compute matches if needed
    //compute normal
    pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhS (new pcl::PointCloud<pcl::FPFHSignature33> ());
    pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhT (new pcl::PointCloud<pcl::FPFHSignature33> ());
    vector<corrTab> corr;
    vector<int> corrNOS,corrNOT;

    if(testMethod>=1 && testMethod<=5){
        //compute fpfh
        cout<<"computing fpfh..."<<endl;
        FPFHComp(cloudS,inlTh,issIdxS,fpfhS);
        FPFHComp(cloudT,inlTh,issIdxT,fpfhT);
        //match features
        cout<<"matching correspodences..."<<endl;
        corrComp(fpfhS,fpfhT,corr,maxCorr,corrNOS,corrNOT);
        cout<<"NO. corr = "<<corr.size()<<endl;
    }

    //start optimization
    Transform3 result;

    auto start = chrono::steady_clock::now();
    if(testMethod == 1){
        cout<<"running FMA+BnB..."<<endl;
        globReg4D(issS, issT, corr, corrNOS, corrNOT, inlTh,rS+rT,result);
    }
    else if(testMethod == 2){
        cout<<"running BnB..."<<endl;
        globReg4D(issS, issT, corr, corrNOS, corrNOT, inlTh,rS+rT,result,false);
    }
    else if(testMethod == 3){
        cout<<"running RANSAC..."<<endl;
        RANSAC(issS,issT,corr,inlTh,result);
    }
    else if(testMethod == 4){
        cout<<"running the Lifting Method (LM)..."<<endl;
        CApp regFun;
        regFun.LM(issS,issT,corr,inlTh,result);
    }
    else if(testMethod == 5){
        cout<<"running the Game Theory Approach (GTA)..."<<endl;
        GameTheoryAlbrarelli regFun;
        regFun.GTReg(issS,issT,corr,inlTh,result);
    }
    else if(testMethod == 6){
        cout<<"running K4PCS..."<<endl;
        pcl::registration::KFPCSInitialAlignment<pcl::PointXYZ, pcl::PointXYZ> kfpcs;
        pcl::PointCloud<pcl::PointXYZ> final;
        kfpcs.setInputSource(issS);
        kfpcs.setInputTarget(issT);
        kfpcs.setNumberOfThreads (1);
        kfpcs.setDelta (inlTh, false);
        kfpcs.setScoreThreshold (0.001);
        kfpcs.align(final);
        Eigen::Matrix<float, 4,4> a = kfpcs.getFinalTransformation();
        for(size_t i=0;i<4;i++){
            for(size_t j=0;j<4;j++){
                result.x[i+j*4] = (double)a(i,j);
            }
        }
    }
    //!uncomment if you have installed and want to use S4PCS
    //    else if(testMethod == 7){
    //        cout<<"running S4PCS..."<<endl;
    //        pcl::Super4PCS<pcl::PointXYZ, pcl::PointXYZ> s4pcs;
    //        pcl::PointCloud<pcl::PointXYZ> final;
    //        s4pcs.setInputSource(issS);
    //        s4pcs.setInputTarget(issT);
    //        s4pcs.options_.delta =  inlTh;
    //        s4pcs.options_.configureOverlap(0.5);
    //        //register
    //        s4pcs.align(final);
    //        Eigen::Matrix<float, 4,4> a = s4pcs.getFinalTransformation();
    //        for(size_t i=0;i<4;i++){
    //            for(size_t j=0;j<4;j++){
    //                result.x[i+j*4] = (double)a(i,j);
    //            }
    //        }
    //    }

    auto end = chrono::steady_clock::now();
    cout<<"runtime in ms = "<<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<endl;

    result = TransMatCompute(result,mCentS,mCentT);
    cout<<"saving result to: "<<fnameOut<<endl;
    saveResult(fnameOut, result);


    //after registration
    //transform point clouds

    Eigen::Matrix4f transform;
    for(size_t i=0;i<4;i++){
        for(size_t j=0;j<4;j++){
            transform(i,j) = result.x[i+4*j];
        }
    }

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloudSReg (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::transformPointCloud (*cloudS, *cloudSReg, transform);


    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer2 (new pcl::visualization::PCLVisualizer ("after registration"));
    viewer2->setBackgroundColor (0, 0, 0);
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> colorSReg(cloudSReg, 255, 0, 0);
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> colorTReg(cloudT, 255, 255, 0);
    viewer2->addPointCloud<pcl::PointXYZ> (cloudSReg, colorSReg, "source cloud");
    viewer2->addPointCloud<pcl::PointXYZ> (cloudT, colorTReg, "target cloud");
    viewer2->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "after registration");

    while (!viewer2->wasStopped ())
    {
        viewer2->spinOnce (100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }
}

void FPADependenceTest(pcl::PointCloud<pcl::PointXYZ>::Ptr cloudS, pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT,
                       double inlTh, int maxCorr, int testNO, vector<int> &sizeResult){
    if(inlTh<=0||maxCorr<=0){
        cout<<"inlier threshold or maximum correspondence number must > 0"<<endl;
        return;
    }

    //voxel grid filter (VGF)
    cout<<"performing voxel grid sampling with grid size = "<<inlTh<<endl;
    VGF(cloudS, cloudS, inlTh);
    VGF(cloudT, cloudT, inlTh);

    //extract ISS
    pcl::PointCloud<pcl::PointXYZ>::Ptr issS (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointIndicesPtr issIdxS(new pcl::PointIndices);
    pcl::PointCloud<pcl::PointXYZ>::Ptr issT (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointIndicesPtr issIdxT(new pcl::PointIndices);
    cout<<"extracting ISS keypoints..."<<inlTh<<endl;
    ISSExt(cloudS,issS,issIdxS,inlTh);
    ISSExt(cloudT,issT,issIdxT,inlTh);
    cout<<"size of issS = "<<issS->size()<<"; size of issT = "<<issT->size()<<endl;

    //translating the center of both point clouds to the origin
    Vector3 centS(0,0,0), centT(0,0,0);
    double rS, rT;
    CentAndRComp(issS,centS,rS);
    CentAndRComp(issT,centT,rT);
    Vector3 mCentS(-centS.x,-centS.y,-centS.z);
    Vector3 mCentT(-centT.x,-centT.y,-centT.z);
    transPC(issS,mCentS);
    transPC(issT,mCentT);

    //compute matches if needed
    //compute normal
    pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhS (new pcl::PointCloud<pcl::FPFHSignature33> ());
    pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhT (new pcl::PointCloud<pcl::FPFHSignature33> ());
    vector<corrTab> corr;
    vector<int> corrNOS,corrNOT;

    //compute fpfh
    cout<<"computing fpfh..."<<endl;
    FPFHComp(cloudS,inlTh,issIdxS,fpfhS);
    FPFHComp(cloudT,inlTh,issIdxT,fpfhT);
    //match features
    cout<<"matching correspodences..."<<endl;
    corrComp(fpfhS,fpfhT,corr,maxCorr,corrNOS,corrNOT);
    cout<<"NO. corr = "<<corr.size()<<endl;

    sizeResult[0] = corr.size();

    //convert data format
    Matrix3X matrixS;
    DataTrans(issS, matrixS);
    Matrix3X matrixT;
    DataTrans(issT, matrixT);

    Transform3 result;

    vector<int> sizeAfFMP(testNO);
    std::srand ( unsigned ( std::time(0) ) );

    //start testing
    for(size_t i=0;i<testNO;i++){
        //re-initialize and shuffle the correspondences
        vector<corrTab> corrCopy(corr);
        // using built-in random generator:
        std::random_shuffle ( corrCopy.begin(), corrCopy.end() );

        vector<int> corrNOSCopy(corrNOS),corrNOTCopy(corrNOT);

        cout<<"NO correspondences before FMP = "<<corrCopy.size()<<endl;
        FPA(matrixS,matrixT,corrCopy,corrNOSCopy,corrNOTCopy,inlTh,result);
        cout<<"NO correspondences after FMP = "<<corrCopy.size()<<endl;
        sizeResult[i+1] = corrCopy.size();
        sizeAfFMP[i] = corrCopy.size();
    }

    //compute average
    int sum = std::accumulate(sizeAfFMP.begin(), sizeAfFMP.end(), 0);
    float mean = sum/((float)testNO);
    float sd = 0;
    for(size_t i=0;i<testNO;i++){
        sd += (sizeAfFMP[i]-mean)*(sizeAfFMP[i]-mean);
    }
    int max = *max_element(sizeAfFMP.begin(), sizeAfFMP.end());
    int min = *min_element(sizeAfFMP.begin(), sizeAfFMP.end());
    sd = sqrt(sd)/(float)testNO;
    cout<<"mean = "<<mean<<"; sd = "<<sd<<"; max = "<<max<<"; min = "<<min;
}

void saveFMPDependence(string fname, std::vector<int> outVec){
    std::ofstream output_file(fname);
    for (const auto &e : outVec) output_file << e << "\n";
}

int main(int argc, char *argv[])
{
    if(argc == 4){
        //INPUT:
        // 1. path to the source point cloud
        string fnameS = argv[1];
        // 2. path to the target point cloud
        string fnameT = argv[2];
        // 3. output transformation file
        string fnameOut = argv[3];
        // 4. inlier threshold (default 0.05)
        double inlTh = 0.1;
        // 5. maximum number of correspondences per source point (default 10)
        int maxCorr = 10;
        // 6. method u want to run (from 1 to 6, default 1)
        //    1. FMP+BnB; 2. BnB; 3.RANSAC 4. LM; 5. GTA; 6. K4PCS; 7. S4PCS
        int testMethod = 1;
        //read point clouds
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloudS (new pcl::PointCloud<pcl::PointXYZ>);
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT (new pcl::PointCloud<pcl::PointXYZ>);
        cout << "reading input point clouds ..."<< endl;
        if (pcl::io::load<pcl::PointXYZ> (fnameS, *cloudS) == -1) //* load the file
        {
            cout << "Couldn't read file: " << fnameS << endl;
            return (-1);
        }
        if (pcl::io::load<pcl::PointXYZ> (fnameT, *cloudT) == -1) //* load the file
        {
            cout << "Couldn't read file: " << fnameT << endl;
            return (-1);
        }
        optimization(cloudS,cloudT,fnameOut,inlTh,maxCorr,testMethod);
    }
    else if(argc == 6){
        //INPUT:
        // 1. path to the source point cloud
        string fnameS = argv[1];
        // 2. path to the target point cloud
        string fnameT = argv[2];
        // 3. output transformation file
        string fnameOut = argv[3];
        // 4. inlier threshold (default 0.1)
        double inlTh = atof( argv[4] );
        // 5. method u want to run (from 1 to 6, default 1)
        //    1. FMP+BnB; 2. BnB; 3. LM; 4. GTA; 5.RANSAC 6. K4PCS; 7. S4PCS
        int testMethod = atoi( argv[5] );

        int maxCorr = 10;
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloudS (new pcl::PointCloud<pcl::PointXYZ>);
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT (new pcl::PointCloud<pcl::PointXYZ>);

        //read point clouds
        cout << "reading input point clouds ..."<< endl;
        if (pcl::io::load<pcl::PointXYZ> (fnameS, *cloudS) == -1) //* load the file
        {
            cout << "Couldn't read file: " << fnameS << endl;
            return (-1);
        }
        if (pcl::io::load<pcl::PointXYZ> (fnameT, *cloudT) == -1) //* load the file
        {
            cout << "Couldn't read file: " << fnameT << endl;
            return (-1);
        }
        optimization(cloudS,cloudT,fnameOut,inlTh,maxCorr,testMethod);
    }
    //test for the dependence of FPA to the order of input correspondences
    else if(argc == 5){

        //INPUT:
        // 1. path to the config file
        string fpathConfig = argv[1];
        // 2. inlier threshold (default 0.05)
        double inlTh = atof( argv[2] );
        // 3. number of test times
        int testNO = atoi( argv[3] );
        // 4. output file folder
        string fpathOut = argv[4];

        //read config file
        vector<int> vecOut(testNO+1);
        vector<pair<int, int>> data;

        int oneNO;
        string fnameConfig = fpathConfig+"RegConfig_Zhipeng.config";
        cout<<fnameConfig<<endl;
        ifstream myfile (fnameConfig);
        cout<<myfile.is_open()<<endl;
        if (myfile.is_open()){
            //read number of pairs
            myfile >> oneNO;
            cout<<oneNO<<endl;
            if(oneNO>0){
                data.resize(oneNO);
                for(size_t i=0;i<oneNO;i++){
                    myfile >> data[i].first;
                    myfile >> data[i].second;
                    cout<<"testing data["<<i<<"] = ("<<data[i].first<<", "<<data[i].second<<")"<<endl;
                    string fnameS, fnameT;
                    if(data[i].first < 10){
                        fnameS = fpathConfig+"s0"+std::to_string(data[i].first)+".ply";
                    }
                    else{
                        fnameS = fpathConfig+"s"+std::to_string(data[i].first)+".ply";
                    }

                    if(data[i].second < 10){
                        fnameT = fpathConfig+"s0"+std::to_string(data[i].second)+".ply";
                    }
                    else{
                        fnameT = fpathConfig+"s"+std::to_string(data[i].second)+".ply";
                    }


                    int maxCorr = 10;
                    pcl::PointCloud<pcl::PointXYZ>::Ptr cloudS (new pcl::PointCloud<pcl::PointXYZ>);
                    pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT (new pcl::PointCloud<pcl::PointXYZ>);

                    //read point clouds
                    cout << "reading input point clouds ..."<< endl;
                    if (pcl::io::load<pcl::PointXYZ> (fnameS, *cloudS) == -1) //* load the file
                    {
                        cout << "Couldn't read file: " << fnameS << endl;
                        return (-1);
                    }
                    if (pcl::io::load<pcl::PointXYZ> (fnameT, *cloudT) == -1) //* load the file
                    {
                        cout << "Couldn't read file: " << fnameT << endl;
                        return (-1);
                    }

                    //start the test
                    FPADependenceTest(cloudS,cloudT,inlTh,maxCorr,testNO,vecOut);

                    //save the result
                    string fOutName = fpathOut+"s"+std::to_string(data[i].first)+"-s"+std::to_string(data[i].second)+".txt";
                    cout<<"saving file to: "<<fOutName<<endl;
                    saveFMPDependence(fOutName, vecOut);
                }
            }
        }

        myfile.close();
    }
    //registering multiple point clouds using FMP+BnB
    else if(argc == 3){
        //Input:
        //1. path to the config file which stores the paths to all point clouds that we want to register
        string fnameC = argv[1];
        //2. the inlier threshold
        double inlTh = atof(argv[2]);
        int maxCorr = 10;
        //read config file
        string line;
        ifstream myfile (fnameC);


        if (myfile.is_open())
        {
            getline (myfile,line);
            int noPC = stoi(line);
            pcl::PointCloud<pcl::PointXYZ>::Ptr cloudS (new pcl::PointCloud<pcl::PointXYZ>);
            pcl::PointCloud<pcl::PointXYZ>::Ptr cloudSReg (new pcl::PointCloud<pcl::PointXYZ>);
            pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT (new pcl::PointCloud<pcl::PointXYZ>);
            pcl::PointCloud<pcl::PointXYZ>::Ptr cloudTReg (new pcl::PointCloud<pcl::PointXYZ>);
            pcl::PointCloud<pcl::PointXYZ>::Ptr issS (new pcl::PointCloud<pcl::PointXYZ>);
            pcl::PointIndicesPtr issIdxS(new pcl::PointIndices);
            pcl::PointCloud<pcl::PointXYZ>::Ptr issT (new pcl::PointCloud<pcl::PointXYZ>);
            pcl::PointIndicesPtr issIdxT(new pcl::PointIndices);
            pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhS (new pcl::PointCloud<pcl::FPFHSignature33> ());
            pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhT (new pcl::PointCloud<pcl::FPFHSignature33> ());
            vector<corrTab> corr;
            vector<int> corrNOS,corrNOT;
            Transform3 result;

            //show point clouds
            //before registration
            boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("registration UI"));
            viewer->setBackgroundColor (0, 0, 0);
            viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "registration UI");

            /* initialize random seed for generating colors */
            srand (time(NULL));

            int cR, cG, cB; //colors
            //read point clouds
            for(size_t i=0;i<noPC;i++)
            {
                cR = rand()%256;
                cG = rand()%256;
                cB = rand()%256;

                getline (myfile,line);
                cout << "reading: "<< line<<" ..."<< endl;
                char cloudName [50];
                sprintf (cloudName, "Point Cloud %d", (int)i+1);
                //read data into source pc, and register it to the target (ignore the first one)
                if(i%2 == 0){
                    if (pcl::io::load<pcl::PointXYZ> (line, *cloudS) == -1) //* load the file
                    {
                        cout << "Couldn't read file: " << line << endl;
                        return (-1);
                    }


                    //voxel grid filter
                    cout<<"performing voxel grid sampling with grid size = "<<inlTh<<endl;
                    VGF(cloudS, cloudS, inlTh);

                    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> colorSOri(cloudS, cR, cG, cB);
                    viewer->addPointCloud<pcl::PointXYZ> (cloudS, colorSOri, cloudName);
                    viewer->spinOnce ();

                    cout<<"extracting ISS keypoints..."<<endl;
                    ISSExt(cloudS,issS,issIdxS,inlTh);
                    cout<<"computing FPFH..."<<endl;
                    FPFHComp(cloudS,inlTh,issIdxS,fpfhS);

                    if(i==0){
                        continue;
                    }
                    else{

                        //match features (source to target)
                        cout<<"matching correspodences..."<<endl;
                        corrComp(fpfhS,fpfhT,corr,maxCorr,corrNOS,corrNOT);
                        cout<<"NO. corr = "<<corr.size()<<endl;

                        //translating the center of both point clouds to the origin
                        Vector3 centS(0,0,0), centT(0,0,0);
                        double rS, rT;
                        CentAndRComp(issS,centS,rS);
                        CentAndRComp(issT,centT,rT);
                        Vector3 mCentS(-centS.x,-centS.y,-centS.z);
                        Vector3 mCentT(-centT.x,-centT.y,-centT.z);
                        transPC(issS,mCentS);
                        transPC(issT,mCentT);

                        //optimization
                        cout<<"running FMP+BnB..."<<endl;
                        globReg4D(issS, issT, corr, corrNOS, corrNOT, inlTh,rS+rT,result);

                        //transform pc
                        result = TransMatCompute(result,mCentS,mCentT);

                        //after registration
                        //transform point clouds
                        Eigen::Matrix4f transform;
                        for(size_t i=0;i<4;i++){
                            for(size_t j=0;j<4;j++){
                                transform(i,j) = result.x[i+4*j];
                            }
                        }

                        pcl::transformPointCloud (*cloudS, *cloudSReg, transform);
                        transPC(issS,centS);
                        pcl::transformPointCloud (*issS, *issS, transform);

                        //add transformed pc into the viewer
                        pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> colorS(cloudSReg, cR, cG, cB);
                        viewer->removePointCloud(cloudName);
                        viewer->addPointCloud<pcl::PointXYZ> (cloudSReg, colorS, cloudName);
                        viewer->spinOnce ();
                    }

                }
                //read data into target pc and register it to the source
                else{
                    if (pcl::io::load<pcl::PointXYZ> (line, *cloudT) == -1) //* load the file
                    {
                        cout << "Couldn't read file: " << line << endl;
                        return (-1);
                    }

                    //voxel grid filter
                    cout<<"performing voxel grid sampling with grid size = "<<inlTh<<endl;
                    VGF(cloudT, cloudT, inlTh);

                    //add transformed pc into the viewer
                    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> colorTOri(cloudT, cR, cG, cB);
                    viewer->addPointCloud<pcl::PointXYZ> (cloudT, colorTOri, cloudName);
                    viewer->spinOnce ();

                    cout<<"extracting ISS keypoints..."<<endl;
                    ISSExt(cloudT,issT,issIdxT,inlTh);

                    cout<<"computing FPFH..."<<endl;
                    FPFHComp(cloudT,inlTh,issIdxT,fpfhT);

                    //match features (target to source)
                    cout<<"matching correspodences..."<<endl;
                    corrComp(fpfhT,fpfhS,corr,maxCorr,corrNOT,corrNOS);
                    cout<<"NO. corr = "<<corr.size()<<endl;

                    //translating the center of both point clouds to the origin
                    Vector3 centS(0,0,0), centT(0,0,0);
                    double rS, rT;
                    CentAndRComp(issS,centS,rS);
                    CentAndRComp(issT,centT,rT);
                    Vector3 mCentS(-centS.x,-centS.y,-centS.z);
                    Vector3 mCentT(-centT.x,-centT.y,-centT.z);
                    transPC(issS,mCentS);
                    transPC(issT,mCentT);

                    cout<<"running FMP+BnB..."<<endl;
                    globReg4D(issT, issS, corr, corrNOT, corrNOS, inlTh,rS+rT,result);
                    //transform pc
                    result = TransMatCompute(result,mCentT,mCentS);

                    //after registration
                    //transform point clouds
                    Eigen::Matrix4f transform;
                    for(size_t i=0;i<4;i++){
                        for(size_t j=0;j<4;j++){
                            transform(i,j) = result.x[i+4*j];
                        }
                    }

                    pcl::transformPointCloud (*cloudT, *cloudTReg, transform);
                    transPC(issT,centT);
                    pcl::transformPointCloud (*issT, *issT, transform);

                    //add transformed pc into the viewer
                    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> colorT(cloudTReg, cR, cG, cB);
                    viewer->removePointCloud(cloudName);
                    viewer->addPointCloud<pcl::PointXYZ> (cloudTReg, colorT, cloudName);
                    viewer->spinOnce ();
                }
            }
            myfile.close();
            while (!viewer->wasStopped())
            {
                viewer->spinOnce (100);
                boost::this_thread::sleep (boost::posix_time::microseconds (100000));
            }
        }
        else{
            cout << "Unable to open config file";
        }


    }
    return 0;
}
