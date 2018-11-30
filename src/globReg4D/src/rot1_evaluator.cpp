#include "../../../include/globReg4D/include/rot1_evaluator.h"
#include "../../../include/globReg4D/include/util_sort.h"
#include "../../../include/globReg4D/include/reg_binary_tree.h"
#include "../../../include/globReg4D/include/reg_common.h"
#include "../../../include/globReg4D/include/geometry.h"


using namespace reg::search;


//!ArcEvaluatorWithCorr_IS
//used to sort interval_array in ascending order
inline bool compareIntervalEnd(const intervalEnd &intA, const intervalEnd &intB)
{
    return (intA.location<intB.location/*||(intA.location==intB.location&&intA.isStart)*/);
}

inline void insertInterval(std::vector<intervalEnd> &intervalArr, const double &startPt, const double &endPt, const int &corrIdx){
intervalEnd intEndTmp;
intEndTmp.formIntervalEnd(startPt,true,corrIdx);
intervalArr.push_back(intEndTmp);
intEndTmp.formIntervalEnd(endPt,false,corrIdx);
intervalArr.push_back(intEndTmp);
}


ArcEvaluatorWithCorr_IS::ArcEvaluatorWithCorr_IS(const Matrix3X &X, const Matrix3X &Y, const Vector &TH,
                                                 const std::vector<corrTab> &corr, const std::vector<int> &corrNOS,
                                                 const std::vector<int> &corrNOT):
    M_in(X), B_in(Y), TH_in(TH), corr_in(corr), _size(0)
{
    int i;

    for(i=0;i<corrNOS.size();i++){
        if(corrNOS[i]>0){
            matchListS_in.push_back(i);
        }
    }
    for(i=0;i<corrNOT.size();i++){
        if(corrNOT[i]>0){
            matchListT_in.push_back(i);
        }
    }
}

ArcEvaluatorWithCorr_IS::~ArcEvaluatorWithCorr_IS()
{
}

int ArcEvaluatorWithCorr_IS::size() const
{
    return _size;
}



void ArcEvaluatorWithCorr_IS::sweep(double &outAngle, int &outUpbnd, bool oneToOne){

    interval_array.clear();
    //save the size for interval stabbing
    _size = M_in.cols();
    //-----------------------------
    //        precomputation
    //-----------------------------

    const int msize_total = (int) M_in.cols();
    const int bsize_total = (int) B_in.cols();
    const int msize = (int)matchListS_in.size();
    const int bsize = (int)matchListT_in.size();

    // Norm of points on XY-plane.
    Vector M_len(msize_total);
    int i,idx;

    for(i=0; i<msize; i++)
    {
        idx = matchListS_in[i];
        double &x = M_in.x[3*idx];
        double &y = M_in.x[3*idx+1];
        M_len(idx)=sqrt(x*x + y*y);
    }

    Vector B_len(bsize_total);
    for(i=0; i<bsize; i++)
    {
        idx = matchListT_in[i];
        double &x = B_in.x[3*idx];
        double &y = B_in.x[3*idx+1];

        B_len(idx)=sqrt(x*x + y*y);
    }

    // Z-coordinate of points.

    Vector M_z(msize_total);
    Vector B_z(bsize_total);

    for(i=0; i<msize; i++)
    {
        M_z(matchListS_in[i])=M_in(2,matchListS_in[i]);
    }

    for(i=0; i<bsize; i++)
    {
        B_z(matchListT_in[i])=B_in(2,matchListT_in[i]);
    }

    // Azimuth of points.
    Vector M_azi(msize_total);
    Vector B_azi(bsize_total);

    for(i=0; i<msize; i++)
    {
        M_azi(matchListS_in[i]) = atan2(M_in(1,matchListS_in[i]), M_in(0,matchListS_in[i]));
    }

    for(i=0; i<bsize; i++)
    {
        B_azi(matchListT_in[i]) = atan2(B_in(1,matchListT_in[i]), B_in(0,matchListT_in[i]));
    }

    //map
    Vector TH = TH_in;

    //-----------------------------
    //        Sweep
    //-----------------------------

    double dz,d,thMz,rth,dev,beg,end;
    for (i=0; i<corr_in.size(); i++)
    {
        int idxS = corr_in[i].idxS;
        int idxT = corr_in[i].idxT;
        dz = B_z(idxT)-M_z(idxS);
        d = B_len(idxT)-M_len(idxS);
        thMz = TH(idxS)*TH(idxS)-dz*dz;
        if(d*d<=thMz){
            rth = sqrt(thMz);
            //insert the intersection interval to int_idxS

            //insert [0,2pi] if M is too short
            if (M_len(idxS)<=DUMMY_PRECISION)
            {
                insertInterval(interval_array,0,TWOPI,idxS);
            }
            else
            {
                dev = reg::geometry::circleintersection(M_len(idxS), B_len(idxT), rth);

                if (fabs(dev-PI) <= DUMMY_PRECISION )
                {
                    /*That could be improved by instead of insert adding 1 to the final quality*/
                    insertInterval(interval_array,0,TWOPI,i);
                }
                else
                {
                    beg = fmod(B_azi(idxT)-dev-M_azi(idxS), TWOPI);
                    if (beg<0)
                    {
                        beg += TWOPI;
                    }
                    end = fmod(B_azi(idxT)+dev-M_azi(idxS), TWOPI);
                    if (end<0)
                    {
                        end += TWOPI;
                    }
                    if (end>=beg)
                    {
                        insertInterval(interval_array,beg,end,idxS);
                    }
                    else
                    {
                        insertInterval(interval_array,beg,TWOPI,idxS);
                        insertInterval(interval_array,0,end,idxS);
                    }
                }
            }
        }
    }

    intervalStab(outAngle,outUpbnd,oneToOne);
}

void ArcEvaluatorWithCorr_IS::intervalStab(double &outAngle, int &outUpbnd, bool oneToOne){
    std::vector<int> ACTab(_size,0);
    int currUpbnd = 0;
    outUpbnd = 0;
    //1. sort interval_array
    std::sort(interval_array.begin(),interval_array.end(),compareIntervalEnd);
    double currLoc = 0;
    int NOEnd = 0;
    if (!oneToOne){
        for(size_t i=0;i<interval_array.size();i++){
            //is a starting point
            if(interval_array[i].isStart){
                ACTab[interval_array[i].corrIdx]++;
                if(ACTab[interval_array[i].corrIdx] == 1){
                    currUpbnd++;
                    if(currUpbnd>outUpbnd){
                        outUpbnd = currUpbnd;
                        outAngle = interval_array[i].location;
                    }
                }
            }
            else{
                ACTab[interval_array[i].corrIdx]--;
                if(ACTab[interval_array[i].corrIdx] == 0){
                    NOEnd++;
                }
            }
            if(interval_array[i].location>currLoc){
                currUpbnd-=NOEnd;
                NOEnd = 0;
                if(currLoc == outAngle){
                    outAngle = (currLoc+interval_array[i].location)/2;
                }
                currLoc = interval_array[i].location;
            }
        }
        currUpbnd-=NOEnd;
    }
    else{
        for(size_t i=0;i<interval_array.size();i++){
            //is a starting point
            if(interval_array[i].isStart){
                ACTab[interval_array[i].corrIdx]++;
                currUpbnd++;
                if(currUpbnd>outUpbnd){
                    outUpbnd = currUpbnd;
                    outAngle = interval_array[i].location;
                }
            }
            else{
                ACTab[interval_array[i].corrIdx]--;
                NOEnd++;
            }
            if(interval_array[i].location>currLoc){
                currUpbnd-=NOEnd;
                NOEnd = 0;
                currLoc = interval_array[i].location;
            }
        }
        currUpbnd-=NOEnd;
    }
}
