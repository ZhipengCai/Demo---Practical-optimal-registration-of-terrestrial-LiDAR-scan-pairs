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

Overview
=========

Given two input point clouds, and a given rotation axis (e.g., for terrestrial LiDAR scans, the rotation axis is the vertical/z-axis), our algorithm works by:

1. Extracting keypoints for both point clouds, via [ISS](https://ieeexplore.ieee.org/document/5457637)

2. Finding candidate matches between keypoints, via [FPFH](https://ieeexplore.ieee.org/document/5152473)

3. Prune the candidate matches via the proposed Fast Match Pruning (FMP, Algorithm 3 in the paper) algorithm. Besides being able to significantly reduce the input matches, FMP guarantees that all the removed matches are guaranteed to be outliers/incorrect matches, so that the global optimal solution remains the same before and after FMP.

4. Finding the global optimal 4DOF pose via an efficient Branch-and-Bound (BnB, Algorithm 2 in the paper). The BnB works on the 3D translation space, which exhaustively searches for the best translation. And a global optimal and polynomial time 1D rotation search algorithm is embedded inside the translation BnB to search for the best rotation given translation. This polynomial time rotation search algorithm is the key to make our BnB algorithm efficient.

The main contribution of this paper lies in step 3 and 4, step 1 and 2 are implemented using [PCL](http://pointclouds.org/) and can be replaced by other methods.

Please refer to the paper for more details.



