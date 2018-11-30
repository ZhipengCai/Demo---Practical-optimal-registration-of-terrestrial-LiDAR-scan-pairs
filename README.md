# Demo---Practical-optimal-registration-of-terrestrial-LiDAR-scan-pairs

Efficient 4DOF (3DOF translation + 1D rotation angle, with known rotation axis) point cloud registration with guaranteed global optimality. 

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

Do not hesitate to contact the authors if you have any question or find any bugs.

Algorithm overview
==================

Given two input point clouds, and a given rotation axis (e.g., for terrestrial LiDAR scans, the rotation axis is the vertical/z-axis), our algorithm works by:

1. Compute candidate matches between two point clouds.

2. Prune the candidate matches via the proposed Fast Match Pruning (FMP, Algorithm 3 in the paper) algorithm. Besides being able to significantly reduce the input matches, FMP guarantees that **all removed matches are guaranteed to be outliers/ incorrect matches**, so that the global optimal solution remains the same before and after FMP.

3. Find the global optimal 4DOF pose by Branch-and-Bound (BnB, Algorithm 2 in the paper). The BnB works on the 3D translation space, which exhaustively searches for the best translation. A **global optimal** and **polynomial time** 1D rotation search algorithm is embedded inside the translation BnB to search for the best rotation given translation. This highly efficient rotation search algorithm makes our BnB very practical and have comparable/ faster speed to prevalent local methods.

The main contribution of this paper is in step 2 and 3.
 
Step 1 is achieved by first extracting [ISS](https://ieeexplore.ieee.org/document/5457637) keypoints between two point clouds, and then computing matches between the keypoints via [FPFH](https://ieeexplore.ieee.org/document/5152473). We use the implementation from [PCL](http://pointclouds.org/) for ISS and FPFH.

Please refer to the [paper](https://www.sciencedirect.com/science/article/pii/S0924271618303125?via%3Dihub) for more details.

Compile
=======


