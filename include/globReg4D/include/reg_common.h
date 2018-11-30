/////////////////////////////////////////////////////////////////////////////
//
//              R O T A T I O N   S E A R C H
//
// This package contains the source code which implements the
// BnB rotation search algorithm and the nested 6 DoF registration
// algorithm proposed in
//
// A. Parra Bustos, T.-J. Chin, A. Eriksson, H. Li and D. Suter
// Fast Rotation Search with Stereographic Projections for 3D Registration
// IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI)
//
// Copyright (c) 2016 Alvaro PARRA BUSTOS (aparra@cs.adelaide.edu.au.)
// School of Computer Science, The University of Adelaide, Australia
// The Australian Center for Visual Technologies
// http://cs.adelaide.edu.au/~aparra
// Please acknowledge the authors by citing the above paper in any academic
// publications that have made use of this package or part of it.
//
/////////////////////////////////////////////////////////////////////////////

#ifndef REG_COMMON_
#define REG_COMMON_

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <string.h>
#include <iostream>
//!this line tried to remove the assertion step
#define NDEBUG
#include <assert.h>
//#include "matrix.h"
//#include "mex.h"
//#include "blas.h"
#include "cblas.h"

#include "math.h"

namespace reg
{
    //#if !defined(_WIN32)
    //#define dgemm dgemm_
    //#endif
    
#define EPSILON std::numeric_limits<double>::epsilon()
#define DUMMY_PRECISION 1e-12
#define ROT_DIV_ACC 1e-9 //accuracy for dividing 1D rotation interval
#define TRANS_DIV_ACC 1e-6 //accuracy for dividing translation block (set for lidar data, should change for synthetic data)
#define NORM_PRECISION 1e-6

#define PI  M_PI
#define TWOPI 2*M_PI

#define MAX(a,b) ((a) > (b) ? a : b)
#define MIN(a,b) ((a) < (b) ? a : b)


struct NCAng1d{ // normal constrained feasible angle for 1D rotation search
    int n;
    double range[8]; // maximum 4 (one range can pass through 0 degree) ranges
    NCAng1d(){
        n = 0;
    }
    NCAng1d &operator = (const NCAng1d &other){
        n = other.n;
        std::copy(other.range, other.range + 6, range);
        return *this;
    }
};

struct corrTab{ //correspondence table, including the indices in corrIdx and the upperbound
    int idxS; //first idx in corrIdx
    int idxT; //second idx in corrIdx
    int upBnd; //upperbound for this correspondences
    int lwBnd;
    double disFeature; //distance between the feature descriptors
};

struct corrTabWithNC{ //union of correspondence table and normal constrains
    corrTab corr;
    NCAng1d NC;
};

    struct Matrix2X
    {
        double *x;//pointer
        size_t n; //number of points
        
        Matrix2X():x(NULL),n(0) {}
        
        Matrix2X(size_t n): n(n)
        {
            assert(n>0);
            x = (double *)malloc(n*2*sizeof(double));
        }
        
        Matrix2X(const Matrix2X& other) : n(other.n)
        {
            assert(n>0);
            x = (double *)malloc(n*2*sizeof(double));
            std::copy(other.x, other.x+2*n, x);
            //            memcpy(x,other.x,n*2*sizeof(double));
        }
        
        ~Matrix2X()
        {
            free(x);
        }
        
        void setSize(size_t newSize)
        {
            assert(newSize>=0 && "size must be >=0");
            n=newSize;
            
            if (n==0)
            {
                free(x);
                x=NULL;
            }
            else
            {
                if(x==NULL)
                {
                    x = (double *)malloc(n*2*sizeof(double));
                }
                else
                {
                    x = (double *)realloc(x,n*2*sizeof(double));
                }
            }
        }
        
        double *getPr() const
        {
            return x;
        }
        
        double &operator()(int i, int j)
        {
            return x[i  + 2*j];
        }
        
        double operator()(int i, int j) const
        {
            return x[i  + 2*j];
        }
        
        //Copy assignment
        Matrix2X &operator=(const Matrix2X &other)
        {
            n=other.n;
            x = (double *)malloc(n*2*sizeof(double));
            //memcpy(x,other.x,n*3*sizeof(double));
            std::copy(other.x, other.x + n*2, x);
            return *this;
        }
        
        // sort using index
        void sort(unsigned int* idx)
        {
            double *tmp = (double *)malloc(n*2*sizeof(double));
            
            size_t i;
            for(i=0; i<n; i++)
            {
                tmp[i*2]   = x[idx[i]*2];
                tmp[i*2+1] = x[idx[i]*2+1];
                tmp[i*2+2] = x[idx[i]*2+2];
            }
            std::copy(tmp, tmp+2*n, x);
            free(tmp);
        }
        
        size_t cols() const
        {
            assert(n>=0 && "invalid number of columns");
            return n;
        }
    };
    
    
    
    
    struct Matrix3X
    {
        double *x;//pointer
        size_t n; //number of points
        
        Matrix3X():x(NULL),n(0) {}
        
        Matrix3X(size_t n): n(n)
        {
            //std::cout<<" n = "<< n <<std::endl;
            
            assert(n>0);
            x = (double *)malloc(n*3*sizeof(double));
        }
        
        Matrix3X(const Matrix3X& other) : n(other.n)
        {
            assert(n>0);
            x = (double *)malloc(n*3*sizeof(double));
            //        memcpy(x,other.x,n*3*sizeof(double));
            std::copy(other.x, other.x+3*n, x);            
        }
        
        ~Matrix3X()
        {
            free(x);
        }
        
        void setSize(size_t newSize)
        {
            assert(newSize>=0 && "size must be >=0");
            n=newSize;
            
            if (n==0)
            {
                free(x);
                x=NULL;
            }
            else
            {
                if(x==NULL)
                {
                    x = (double *)malloc(n*3*sizeof(double));
                }
                else
                {
                    x = (double *)realloc(x,n*3*sizeof(double));
                }
            }
        }
        
        double *getPr() const
        {
            return x;
        }
        
        double &operator()(int i, int j)
        {
            return x[i  + 3*j];
        }
        
        double operator()(int i, int j) const
        {
            return x[i  + 3*j];
        }
        
        //Copy assignment
        Matrix3X &operator=(const Matrix3X &other)
        {
            n=other.n;
            x = (double *)malloc(n*3*sizeof(double));
            //memcpy(x,other.x,n*3*sizeof(double));
            std::copy(other.x, other.x + n*3, x);
            return *this;
        }
        
        // sort using index
        void sort(unsigned int* idx)
        {
            double *tmp = (double *)malloc(n*3*sizeof(double));
            
            size_t i;
            for(i=0; i<n; i++)
            {
                tmp[i*3]   = x[idx[i]*3];
                tmp[i*3+1] = x[idx[i]*3+1];
                tmp[i*3+2] = x[idx[i]*3+2];
            }
            //memcpy(x, tmp, n*3*sizeof(double));
            std::copy(tmp, tmp+3*n, x);
            free(tmp);
        }
        
        size_t cols() const
        {
            assert(n>=0 && "invalid number of columns");
            return n;
        }
    };
    
    
    struct Matrix4X
    {
        double *x;//pointer
        size_t n; //number of points
        
        Matrix4X():x(NULL),n(0) {}
        
        Matrix4X(int n): n(n)
        {
            assert(n>0 && "");
            x = (double *)malloc(n*4*sizeof(double));
        }
        
        Matrix4X(const Matrix4X& other) : n(other.n)
        {
            assert(n>0 && "");
            x = (double *)malloc(n*4*sizeof(double));
            //            memcpy(x,other.x,n*4*sizeof(double));
            std::copy(other.x, other.x+4*n, x);
        }
        
        ~Matrix4X()
        {
            free(x);
        }
        
        void setSize(size_t size)
        {
            assert(size>0 && "size must be >0");
            n=size;
            if (n==0)
            {
                free(x);
                x=NULL;
            }
            else
            {
                if(x==NULL)
                {
                    x = (double *)malloc(n*4*sizeof(double));
                }
                else
                {
                    x = (double *)realloc(x,n*4*sizeof(double));
                }
            }
        }
        
        double &operator()(int i, int j)
        {
            return x[i  + 4*j];
        }
        
        double operator()(int i, int j) const
        {
            return x[i  + 4*j];
        }
        
        //Copy assignment
        Matrix4X &operator=(const Matrix4X &other)
        {
            n=other.n;
            x = (double *)malloc(n*4*sizeof(double));
            memcpy(x,other.x,n*4*sizeof(double));
            return *this;
        }
        
        // sort using index
        void sort(unsigned int* idx)
        {
            double *tmp = (double *)malloc(n*4*sizeof(double));
            
            unsigned int i;
            for(i=0; i<n; i++)
            {
                tmp[i*4]   = x[idx[i]*4];
                tmp[i*4+1] = x[idx[i]*4+1];
                tmp[i*4+2] = x[idx[i]*4+2];
                tmp[i*4+3] = x[idx[i]*4+3];
                
            }
            memcpy(x, tmp, n*4*sizeof(double));
            free(tmp);
        }
        
        size_t cols() const
        {
            assert(n>=0 && "invalid number of columns");
            return n;
        }
    };
    
    
    
    struct Vector
    {
        double *x;//pointer
        size_t n; //size
        
        Vector():x(NULL),n(0) {}
        
        
        Vector(size_t n): n(n)
        {
            assert(n>0);
            x = (double *)malloc(n*sizeof(double));
        }
        
        Vector(const Vector& other): n(other.n)
        {
            x = (double *)malloc(n*sizeof(double));
            memcpy(x,other.x,n*sizeof(double));
        }
        
        Vector &operator=(const Vector &other)
        {
            n=other.n;
            x = (double *)malloc(n*sizeof(double));
            memcpy(x,other.x,n*sizeof(double));
            return *this;
        }
        
        ~Vector()
        {
            free(x);
        }
        
        double &operator()(int i) const
        {
            return x[i];
        }
        
        
        void setSize(size_t size)
        {
            assert(size>=0 && "size must be >0");
            if (n==0)
            {
                free(x);
                x=NULL;
            }
            else
            {
                if(x==NULL)
                {
                    x = (double *)malloc(n*sizeof(double));
                }
                else
                {
                    x = (double *)realloc(x,n*sizeof(double));
                }
            }
        }
        
    };
    
    struct Matrix3
    {
        double m[9];
        double *getPr()
        {
            return m;
        }
    };
    
    
    struct Matrix2
    {
        double m[4];
        double *getPr()
        {
            return m;
        }
    };
    
    
    struct Vector3
    {
        union
        {
            double m[3];
            struct
            {
                double x, y, z;
            };
        };
        
        Vector3(double x, double y, double z): x(x),y(y),z(z) {}
        
        double norm() const
        {
            return sqrt(x*x + y*y + z*z);
        }
        
        double sqrdNorm() const
        {
            return x*x + y*y + z*z;
        }
        
        double operator[](int i) const
        {
            return m[i];
        }
        
        double *getPr()
        {
            return m;
        }
        
        
        Vector3 operator+( const Vector3& other )
        {
            Vector3 resp(this->x + other.x, this->y + other.y, this->z + other.z);
            return resp;
        }
        
        Vector3 operator-( const Vector3& other )
        {
            Vector3 resp(this->x - other.x, this->y - other.y, this->z - other.z);
            return resp;
        }
    };
    
    struct Vector2
    {
        union
        {
            double m[2];
            struct
            {
                double x, y;
            };
        };
        
        Vector2(double x, double y): x(x),y(y){}
        
        double norm() const
        {
            return sqrt(x*x + y*y);
        }
        
        double sqrdNorm() const
        {
            return x*x + y*y ;
        }
        
        double operator[](int i) const
        {
            return m[i];
        }
        
        double *getPr()
        {
            return m;
        }
        
        Vector2 operator+( const Vector2& other )
        {
            Vector2 resp(this->x + other.x, this->y + other.y);
            return resp;
        }
        
        Vector2 operator-( const Vector2& other )
        {
            Vector2 resp(this->x - other.x, this->y - other.y);
            return resp;
        }
    };
    
    
    //Multiply a matrix and a 3-vector
    inline Vector3 multiply(const Matrix3 &M, const double *v)
    {
        double x,y,z;
        const double *m = M.m;
        x = v[0]*m[0] + v[1]*m[3] + v[2]*m[6];
        y = v[0]*m[1] + v[1]*m[4] + v[2]*m[7];
        z = v[0]*m[2] + v[1]*m[5] + v[2]*m[8];
        return Vector3(x,y,z);
    }
    
    //Multiply a matrix and a 2-vector
    inline Vector2 multiply(const Matrix2 &M, const double *v)
    {
        double x,y;
        const double *m = M.m;
        x = v[0]*m[0] + v[1]*m[2];
        y = v[0]*m[1] + v[1]*m[3];
        return Vector2(x,y);
    }
    
    
    struct Transform3
    {
        double x[16];
        
        Transform3()
        {
            x[0]=1; x[4]=0;  x[8]=0; x[12]=0;
            x[1]=0; x[5]=1;  x[9]=0; x[13]=0;
            x[2]=0; x[6]=0; x[10]=1; x[14]=0;
            x[3]=0; x[7]=0; x[11]=0; x[15]=1;
        }
        
        Transform3(Matrix3 &R, Vector3 &t)
        {
            const double *r = R.m;
            
            x[0]=r[0]; x[4]=r[3];  x[8]=r[6]; x[12]=t.x;
            x[1]=r[1]; x[5]=r[4];  x[9]=r[7]; x[13]=t.y;
            x[2]=r[2]; x[6]=r[5]; x[10]=r[8]; x[14]=t.z;
            x[3]=0;    x[7]=0;    x[11]=0;    x[15]=1;
        }
        
        Transform3(Matrix3 &R)
        {
            const double *r = R.m;
            
            x[0]=r[0]; x[4]=r[3];  x[8]=r[6]; x[12]=0;
            x[1]=r[1]; x[5]=r[4];  x[9]=r[7]; x[13]=0;
            x[2]=r[2]; x[6]=r[5]; x[10]=r[8]; x[14]=0;
            x[3]=0;    x[7]=0;    x[11]=0;    x[15]=1;
        }
        
        Transform3(Matrix2 &R)
        {
            const double *r = R.m;
            
            x[0]=r[0]; x[4]=r[2]; x[ 8]=0.0; x[12]=0.0;
            x[1]=r[1]; x[5]=r[3]; x[ 9]=0.0; x[13]=0.0;
            x[2]=0.0;  x[6]=0.0;  x[10]=1.0; x[14]=0.0;
            x[3]=0.0;  x[7]=0.0;  x[11]=0.0; x[15]=1.0;
        }
        
        Transform3(Vector3 &t)
        {
            x[0]=1; x[4]=0;  x[8]=0; x[12]=t.x;
            x[1]=0; x[5]=1;  x[9]=0; x[13]=t.y;
            x[2]=0; x[6]=0; x[10]=1; x[14]=t.z;
            x[3]=0; x[7]=0; x[11]=0; x[15]=1;
        }
        
        Transform3 operator*(Transform3& tform )
        {
            //char *chn = (char *)"N";
            double alpha = 1.0, beta = 0.0;
            int m = 4;
            
            Transform3 resp; //TODO: avoid to call default constructor
            
            //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, m, m, alpha, this->x, m, tform.x, m, beta, resp.x, m);
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, alpha, this->x, m, tform.x, m, beta, resp.x, m);

            return resp;
        }
        
        double *getPr()
        {
            return x;
        }
    };
    
    
    
    /**
     * @brief matrixMultipy Computes C = R*X.
     * @param R Rotation matrix.
     * @param X 3xn matrix where each column is a 3D point.
     * @param C 3xn output matrix. C = R*X.
     * @param n Number of points.
     */
    inline
    void matrixMultipy(double *R, double*X, double *C ,int n)
    {
        assert(n>0 && "");
        //char *chn = (char *)"N";
        double alpha = 1.0, beta = 0.0;
        int m = 3;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, m, alpha, R, m, X, m, beta, C, m);
    }
    
    
    //Norm of a vector
    inline
    double norm(const double *x)
    {
        return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    }
    
    
    struct AxisAngle
    {
        double x, y, z; //axis
        double w;       //angle
        AxisAngle(){}
        AxisAngle(Vector3 &axis, double angle): x(axis.x), y(axis.y), z(axis.z), w(angle) {}
        AxisAngle(const Vector3 &axis, const double angle): x(axis.x), y(axis.y), z(axis.z), w(angle) {}
    };
    
    
    inline
    void fromAxisAngle(Matrix3 &R, const AxisAngle &a)
    {
        assert(fabs(sqrt(a.x*a.x + a.y*a.y + a.z*a.z)-1.0)<DUMMY_PRECISION &&
                 "axis must be an unit-vector");
        
        // Matrix3 R;
        
        //Convert to quaternion
        const double s = sin(0.5*a.w);
        const double w = cos(0.5*a.w);
        const double x = s*a.x;
        const double y = s*a.y;
        const double z = s*a.z;
        
        //Convert to matrix
        const double tx  = 2.0*x;
        const double ty  = 2.0*y;
        const double tz  = 2.0*z;
        const double twx = tx*w;
        const double twy = ty*w;
        const double twz = tz*w;
        const double txx = tx*x;
        const double txy = ty*x;
        const double txz = tz*x;
        const double tyy = ty*y;
        const double tyz = tz*y;
        const double tzz = tz*z;
        
        R.m[0] = 1.0-(tyy+tzz);
        R.m[1] = txy+twz;
        R.m[2] = txz-twy;
        R.m[3] = txy-twz;
        R.m[4] = 1.0-(txx+tzz);
        R.m[5] = tyz+twx;
        R.m[6] = txz+twy;
        R.m[7] = tyz-twx;
        R.m[8] = 1.0-(txx+tyy);
        
        // return R;
    }
    
    
    struct PointCloud
    {
        const Matrix3X &pts;
        
        PointCloud(const Matrix3X &pts): pts(pts) {}
        
        // Must return the number of data points
        inline size_t kdtree_get_point_count() const { return pts.cols(); }
        
        // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
        inline double kdtree_distance(const double *p1, const size_t idx_p2,size_t /*size*/) const
        {
            const double d0=p1[0]-pts.x[3*idx_p2];
            const double d1=p1[1]-pts.x[3*idx_p2+1];
            const double d2=p1[2]-pts.x[3*idx_p2+2];
            return d0*d0+d1*d1+d2*d2;
        }
        
        // Returns the dim'th component of the idx'th point in the class:
        // Since this is inlined and the "dim" argument is typically an immediate value, the
        //  "if/else's" are actually solved at compile time.
        inline double kdtree_get_pt(const size_t idx, int dim) const
        {
            if (dim==0) return pts.x[3*idx];
            if (dim==1) return pts.x[3*idx+1];
            return pts.x[3*idx+2];
        }
        
        // Optional bounding-box computation: return false to default to a standard bbox computation loop.
        //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
        //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
        template <class BBOX>
        bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
        
    };

    
} // End namespace reg

#endif // REG_COMMON_
