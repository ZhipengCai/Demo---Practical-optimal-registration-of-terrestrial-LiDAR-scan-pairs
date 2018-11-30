#include "../../include/globReg4D/GameTheoryAlbrarelli.h"

GameTheoryAlbrarelli::GameTheoryAlbrarelli()
{

}

double dist(const double *x, const double *y)
{
    double r, d = 0;

    for (int i=0; i<3; i++)
    {
        r = x[i]-y[i];
        d += r*r;
    }

    d = sqrt(d);
    return d;
}

void payoff(const pcl::PointCloud<pcl::PointXYZ>::Ptr X,
            const pcl::PointCloud<pcl::PointXYZ>::Ptr Y,
            const float &th, const bool &oto,
            const std::vector<reg::corrTab> &corr,
            double *C){
    int n = corr.size();

    double dx, dy;
    int i,j;

    for(i = 0;i<n*n;i++){
        C[i] = 0;
    }
    for (i=0; i<n; i++)
    {
        reg::corrTab p = corr.data()[i];

        for (j = i+1; j<n; j++)
        {

            reg::corrTab q = corr.data()[j];

            if  (oto && ( (p.idxS == q.idxS) || (p.idxT == q.idxT) )){
                continue;
            }

            double xp[3] = {X->points[p.idxS].x, X->points[p.idxS].y, X->points[p.idxS].z};
            double xq[3] = {X->points[q.idxS].x, X->points[q.idxS].y, X->points[q.idxS].z};


            double yp[3] = {Y->points[p.idxT].x, Y->points[p.idxT].y, Y->points[p.idxT].z};
            double yq[3] = {Y->points[q.idxT].x, Y->points[q.idxT].y, Y->points[q.idxT].z};


            dx = dist(xp,xq);
            dy = dist(yp,yq);

            if (th > 0)
            {
                if (fabs(dx-dy) <= th)
                    C[n*j+i] = C[n*i+j] = 1;
            }
            else
            {
                if (dx>dy)
                    C[n*j+i] = C[n*i+j] = dy/dx;
                else
                    C[n*j+i] = C[n*i+j] = dx/dy;
            }
        }
    }
}

template <typename T>
inline void mult(const T *A,const T *x, int size, T *y)
{
    const T *ptr;
    double yv;
    for(int j=0;j<size;++j,++y){
        ptr=x;
        //*y=0;
        yv = 0;
        for(int i=0;i<size;++i,++ptr,++A){
            yv += (*A)*(*ptr);
        }
        *y = static_cast<T>(yv);
    }
}

template <typename T>
inline void simplexify(T *x, int size){
    T sum=0;
    T *ptr=x;
    for(int i=0;i<size;++i,++ptr)
        if(*ptr>=0)
            sum+=*ptr;
        else
            *ptr=0;
    for(int i=0;i<size;++i,++x)
        *x/=sum;
}

template <typename T>
inline double dot(const T *x, const T *y, int size){
  double sum=0;
  for(int i=0;i<size;++i,++x,++y)
     sum+=*x**y;
  return sum;
}

template <typename T>
inline void scale(T *x, int size, double c){
  for(int i=0;i<size;++i,++x)
     *x*=c;
}

template <typename T>
inline void linear_comb(const T *x, T *y, int size,double alfa){
  for(int i=0;i<size;++i,++x,++y)
     *y=alfa*(*x-*y)+*y;
}

template <typename T>
inline double nash_error(const T *x, const T *Ax, double xAx, int size)
{
  double sum=0;
  const T *p_x=x;
  const T *p_Ax=Ax;
  double tmp;
  for(int i=0;i<size;++i,++p_x,++p_Ax){
     tmp=xAx-*p_Ax;
     if(tmp>*p_x)
        tmp=*p_x;
     sum+=tmp*tmp;
  }
  return sum;
}

template <typename T>
inline double selectStrategy(const T *Ax, const T *x, int size, double &delta, int &idx)
{
  const T *ptr1,*ptr2;
  double maxv = -std::numeric_limits<double>::infinity();
  double minv = std::numeric_limits<double>::infinity();

  int max_idx=-1,min_idx=-1;
  ptr1=Ax;
  ptr2=x;

  for(int i=0;i<size;++i,++ptr1,++ptr2){
     if(*ptr1>maxv){
        maxv=*ptr1;
        max_idx=i;
     }
     if(*ptr2>0&&*ptr1<minv){
        minv=*ptr1;
        min_idx=i;
     }
  }
  double xAx=dot(Ax,x,size);

  maxv-=xAx;
  minv=xAx-minv;

  idx=max_idx;
  delta=maxv;
  if(maxv<minv){
     idx=min_idx;
     delta=-minv;
  }

  return nash_error(x,Ax,xAx,size);
}

template <typename T>
void iidyn(const T* A, T* x, int size, T toll, int max_iters, double* out_iters, double* out_err)
{
  int niter=0;
  /* Calculate Ax */
  T *Ax=new T[size];
  simplexify(x,size);
  mult(A,x,size,Ax);

  /*** OCCHIO **/
  toll*=toll;

  double delta = 0.;
  double err = std::numeric_limits<double>::max();

  while (niter < max_iters)
  {
     int idx=-1;

     err = selectStrategy(Ax,x,size,delta,idx);
     if (err < toll)
        break;

     double den=A[idx*(size+1)]-Ax[idx]-delta;
     bool do_remove=false;
     double mu,tmp;
     if(delta>=0)
     {
        mu=1;
        if(den<0)
        {
           tmp=-delta/den;
           if(mu>tmp)mu=tmp;
           if(mu<0)mu=0;
        }
     }
     else
     {
        //mu=1.0-1.0/(1-x[idx]);
        mu=x[idx]/(x[idx]-1);
        do_remove=true;
        if(den<0)
        {
           tmp=-delta/den;
           if(mu<tmp)
           {
              mu=tmp;
              do_remove=false;
           }
           if(mu>0)mu=0;
        }
     }
     scale(x,size,1-mu);
     x[idx]=do_remove?0:x[idx]+mu;

     simplexify(x,size);

     linear_comb(A+idx*size, Ax, size, mu);

     ++niter;
  }

  delete[] Ax;

  *out_iters = double(niter);
  *out_err = err;
}

void GameTheoryAlbrarelli::GTReg(pcl::PointCloud<pcl::PointXYZ>::Ptr cloudS, pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT,
                                 vector<corrTab> &corr,
                                 double inlTh, reg::Transform3 &result){

    int n = corr.size();
    double *payoffMat = new double[n*n];
    //compute the payoff matrix
    payoff(cloudS, cloudT, inlTh, true, corr, payoffMat);
    //initialize the mixed strategies
    double *mixedS = new double[n];
    double sumx = 0;
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,0.03);
    for(size_t i = 0;i<n;i++){
        mixedS[i] = std::fabs((double)1.0/n+distribution(generator));
        sumx = sumx+mixedS[i];
    }

    for(size_t i=0;i<n;i++){
        mixedS[i]/=sumx;
    }
    int maxIters = 10000;
    double toll = 1e-8;
    double outIters, outErr;
    iidyn<double>(payoffMat, mixedS, n, toll, maxIters, &outIters, &outErr);

    int numCon = 0;
    for(size_t i=0;i<n;i++){
        if(mixedS[i]>1e-9){
            numCon++;
        }
    }

    // Check if the estimation with correspondences gives the same results
    Eigen::Matrix4f T_SVD;
    pcl::Correspondences corrPCL;
    corrPCL.reserve (numCon);
    for (size_t i=0; i<n; i++){
        if(mixedS[i]>1e-9)
        corrPCL.push_back (pcl::Correspondence (corr.data()[i].idxS, corr.data()[i].idxT, 0.f));
    }
    const pcl::registration::TransformationEstimationSVD<pcl::PointXYZ, pcl::PointXYZ> trans_est_svd;
    trans_est_svd.estimateRigidTransformation(*cloudS, *cloudT, corrPCL, T_SVD);


    for(size_t i=0;i<4;i++){
        for(size_t j=0;j<4;j++){
            result.x[j+i*4] = T_SVD(j,i);
        }
    }

    delete[] mixedS;
    delete[] payoffMat;
}
