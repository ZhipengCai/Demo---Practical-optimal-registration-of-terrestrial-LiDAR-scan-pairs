# Demo---Practical-optimal-registration-of-terrestrial-LiDAR-scan-pairs

About
=====
Efficient 4DOF (3DOF translation + 1D rotation angle, with known rotation axis) point cloud registration with guaranteed global optimality. 

![alt text](https://github.com/ZhipengCai/ZhipengCai.github.io/blob/master/papers/ISPRSJ18.jpg " ")

[[Paper](https://www.sciencedirect.com/science/article/pii/S0924271618303125?via%3Dihub)] 

[[Video demo](https://www.youtube.com/watch?v=MKzSN4bbs1o&feature=youtu.be)]

The demo is free for non-commercial academic use. Any commercial use is strictly 
prohibited without the authors' consent. Please acknowledge the authors by citing:

```
@article{cai2019practical,
  title={Practical optimal registration of terrestrial LiDAR scan pairs},
  author={Cai, Zhipeng and Chin, Tat-Jun and Bustos, Alvaro Parra and Schindler, Konrad},
  journal={ISPRS Journal of Photogrammetry and Remote Sensing},
  volume={147},
  pages={118--131},
  year={2019},
  publisher={Elsevier}
}
```
in any academic publications that have made use of this package or part of it.

Contact
=======
Homepage:[https://zhipengcai.github.io/](https://zhipengcai.github.io/) 

Email: zhipeng.cai@adelaide.edu.au

Do not hesitate to contact the authors if you have any question or find any bugs :)

Algorithm overview
==================

Given two input point clouds, and a given rotation axis (e.g., for terrestrial LiDAR scans, the rotation axis is the vertical/z-axis), our algorithm works by:

1. Compute candidate matches between two point clouds.

2. Prune the candidate matches via the proposed Fast Match Pruning (FMP, Algorithm 3 in the paper) algorithm. Besides being able to significantly reduce the input matches, FMP guarantees that **all removed matches are guaranteed to be outliers/ incorrect matches**, so that the global optimal solution remains the same before and after FMP.

3. Find the global optimal 4DOF pose by Branch-and-Bound (BnB, Algorithm 2 in the paper). The BnB works on the 3D translation space, which exhaustively searches for the best translation. A **global optimal** and **polynomial time** 1D rotation search algorithm is embedded inside the translation BnB to search for the best rotation given translation. This highly efficient rotation search algorithm makes our BnB very practical and have comparable/ faster speed to prevalent local methods.

The main contribution of this paper is in step 2 and 3.
 
Step 1 is achieved by first extracting [ISS](https://ieeexplore.ieee.org/document/5457637) keypoints between two point clouds, and then computing matches between the keypoints via [FPFH](https://ieeexplore.ieee.org/document/5152473). We use the implementation from [PCL](http://pointclouds.org/) for ISS and FPFH.

Please refer to the [paper](https://www.sciencedirect.com/science/article/pii/S0924271618303125?via%3Dihub) for more details.

Getting Started
===============
The demo has been tested on Linux (Ubuntu 14.04 LTS 64-bit). All methods reported in the experiments are included and are implemented in C++.

-------------
Pre-requisite
-------------
1. [CMake](https://cmake.org/) version 2.8 or above.

2. [PCL](http://pointclouds.org/) version 1.8 or above.

3. [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) version 3.2 or above.

4. [Cblas](http://www.netlib.org/blas/). There are also other implementations of Cblas, all should work.

5. (Optional) [Super4PCS library](http://nmellado.github.io/Super4PCS/a05034.html). Only for enabling the Super4PCS in the demo. If you don't want to test it, just ignore this item.

-------
Compile
-------
1. Clone this repository.

2. Open the terminal and enter the root directory of this repository.

3. Type "mkdir build && cd build && cmake .. && make"

4. Run the program using the command:

"./4DOFReg ../data/arch/s01.pcd ../data/arch/s02.pcd ../data/result1-2.txt 0.1 1"

+ Command explanation:

"./4DOFReg ../data/arch/s01.pcd (path to the source point cloud, **the point cloud file can be formatted in ".ply/.pcd"**) ../data/arch/s02.pcd (path to the target point cloud) ../data/result1-2.txt (path to the output file, need to mkdir if the folder does not exist) 0.1 (inlier threshold, any value between 0.05 to 0.2 should be fine for the reported real-world datasets) 1 (the method u want to run)"

--------------------------------
Corresponding ID for each method
--------------------------------

1. FMP+BnB

2. BnB

3. RANSAC

4. [LM](http://vladlen.info/papers/fast-global-registration.pdf)

5. [GTA](http://www.dsi.unive.it/~rodola/cvpr2010.pdf)

6. [K4PCS](https://www.ethz.ch/content/dam/ethz/special-interest/baug/igp/photogrammetry-remote-sensing-dam/documents/pdf/thei_weg_schind_IJPRS2013.pdf)

7. [S4PCS](http://geometry.cs.ucl.ac.uk/projects/2014/super4PCS/super4pcs_low.pdf) (i.e., Super4PCS). Disabled by default, see below for how to enable.

--------------------------------
Enable Super4PCS
--------------------------------

To enable S4PCS (i.e. Super4PCS), after installing its library, you need to uncomment the corresponding code in "./CMakeLists.txt" and "./src/main.cpp" and recompile, as follows.

+ The part to uncomment (Just remove one "#" sign for each line) in "CMakeLists.txt" file:

```
## Find super4pcs(used after installing super4pcs), install super4pcs library and uncomment if you want to use it.
##set this to the cmake sub-folder of your S4PCS installation dir. And change the filename of "Super4PCSConfig.cmake" to "FindSuper4PCS.cmake"
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/S4PCS/Super4PCS/build/install/lib/cmake")
#find_package(Super4PCS REQUIRED)
#include_directories(${Super4PCS_INCLUDE_DIR})
#link_directories(${Super4PCS_LIB_DIR})
#target_link_libraries(${PROJECT_NAME} ${Super4PCS_LIBRARIES})
```
  
+ The part to uncomment in "main.cpp" file (**remember to reset the link to the Super4PCS library in your computer**):
```
//!uncomment if you have installed and want to use S4PCS (please reset to the directory of the installed library in your computer)
//#include "../S4PCS/Super4PCS/build/install/include/pcl/registration/super4pcs.h"
```
and
 
```
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
```

