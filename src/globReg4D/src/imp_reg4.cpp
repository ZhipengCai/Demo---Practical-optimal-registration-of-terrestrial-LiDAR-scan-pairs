//4Dof rotation search implemented by Zhipeng

#include <iostream>
#include "../../../include/globReg4D/include/registration.h"
#include "../../../include/globReg4D/include/state_priority_hashtable.h"
#include "../../../include/globReg4D/include/search.h"
#include "../../../include/globReg4D/include/rot1_evaluator.h"
#include <fstream>

//qt libs

using namespace reg;
using namespace reg::search;
typedef TranslationSearchSpaceRegion3DOF TSSR;

#define BRANCHING_FACTOR 8

inline
void translate(Matrix3X &TX, const Matrix3X &X, const Vector3 &t)
{
    const size_t n = X.cols();
    double *tx = TX.x;
    double *x = X.x;

    for (size_t i=0;i<n; i++)
    {
        tx[i*3]   = x[i*3]  +t.x;
        tx[i*3+1] = x[i*3+1]+t.y;
        tx[i*3+2] = x[i*3+2]+t.z;
    }
}

//translate only the points within the matchlist
inline
void translateML(Matrix3X &TX, const Matrix3X &X, const Vector3 &t, int *matchList, int matchListSize){
    double *tx = TX.x;
    double *x = X.x;
    int idx;
    for(size_t i = 0;i<matchListSize;i++)
    {
        idx = matchList[i];
        tx[idx*3] = x[idx*3] + t.x;
        tx[idx*3+1] = x[idx*3+1] + t.y;
        tx[idx*3+2] = x[idx*3+2] + t.z;
    }
}



double computeDis(const Matrix3X &M, const Matrix3X &B, const int idxM, const int idxB){
    double disX = M.x[3*idxM]-B.x[3*idxB];
    double disY = M.x[3*idxM+1]-B.x[3*idxB+1];
    double disZ = M.x[3*idxM+2]-B.x[3*idxB+2];
    return sqrt(disX*disX+disY*disY+disZ*disZ);
}


Matrix4X reg::search::transform( Matrix4X &x, Transform3 &tform)
{
    //TM = guessAndResult * TM;
    //y = tform*x;
    //char *chn = (char *)"N";
    double alpha = 1.0, beta = 0.0;
    int m = 4;
    int n =x.cols();

    assert(n>0);

    Matrix4X y(n);

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, m, alpha, tform.x, m, x.x, m, beta, y.x, m);

    return y;
}

int reg::search::nestedbnb_search_4dof_withCorr_IS(const Matrix3X &M, const Matrix3X &B, std::vector<corrTab> &corr,
                                                     std::vector<int> &corrNOS, std::vector<int> &corrNOT,
                                                     Transform3 &guessAndResult, double th, double tolZ,
                                                     int gap, int lwbnd, TranslationSearchSpaceRegion3DOF &tr_region, bool enforceOneToOne)
{
    int i,j,k, np;
    int state_lwbnd, upbnd;

    // Counters of bounds' evaluation calls
    int countUpbndEval, countLwbndEval;

    Vector3 trCentre(0,0,0);

    const size_t msize = M.cols();
    const int buckets = MAX((int)msize/10,10);

    StatePriorityHashtable<TSSR, int, SearchState> table(buckets); //init priority queue

    // Initialise
    countUpbndEval = countLwbndEval= 0;

    trCentre = ssrCentre(tr_region);
    Matrix3X TM(msize);
    translate(TM, M, trCentre);

    AxisAngle rsearchResult;
    upbnd = rot1_withCorr_IntervalStab(TM,B,corr,corrNOS,corrNOT,
                          th+ssrUncertainty(tr_region),tolZ,lwbnd,0,rsearchResult,enforceOneToOne);

    countUpbndEval++;
    assert( upbnd>=lwbnd && "Bound error");

    if (upbnd-lwbnd <= gap)
    {
        return lwbnd;
    }

    TSSR **ssrArray = new TSSR*[BRANCHING_FACTOR];

    SearchState<TSSR, int> *state;
    state = new SearchState<TSSR>(tr_region, upbnd);
    table.push(state);

    int iter = 0;

    while (table.size()){
        iter++;
        state = table.pop();
        // Evaluate lower boud
        // Obtain lower bound calling the rotsearch with tr_c and 0 uncertainty trans.
        trCentre = ssrCentre(state->ssr);
        translate(TM, M, trCentre);

        state_lwbnd = rot1_withCorr_IntervalStab(TM, B, corr, corrNOS, corrNOT, th,
                                                 tolZ, lwbnd, 0, rsearchResult,enforceOneToOne);

        countLwbndEval++;

        // Update solution
        if (state_lwbnd > lwbnd)
        {
            lwbnd = state_lwbnd;

            Matrix3 R;
            fromAxisAngle(R,rsearchResult);
            Transform3 Rtform(R);
            Transform3 Ttform(trCentre);

            guessAndResult = Rtform*Ttform;

            table.prune(lwbnd);
        } // End update solution

        // Stopping criterion
        if (state->bnd - lwbnd <= gap)
        {
            delete state;
            break;
        }

        // Branch
        np = reg::search::split(state->ssr, ssrArray);
        delete state;
        for(i=0; i<np; i++)
        {
            // Compute upper bound by calling rotsearch with tc and t uncertainty delta
            trCentre = ssrCentre( *(ssrArray[i]) );

            translate (TM, M, trCentre);

            upbnd = rot1_withCorr_IntervalStab(TM, B, corr, corrNOS, corrNOT, th+ssrUncertainty(*(ssrArray[i])), tolZ,
                                  lwbnd, 0, rsearchResult,enforceOneToOne);

            if ( upbnd > lwbnd )
            {
                state = new SearchState<TSSR>(*(ssrArray[i]), upbnd);
                table.push(state);
            }
            delete ssrArray[i];
        }

        countUpbndEval += np;
    } //bnb loop

    delete []ssrArray;
    return lwbnd;
}

int reg::search::computeConsensus_withCorr(const Matrix3X &M, const Matrix3X &B, std::vector<corrTab> &corr, Transform3 &guessAndResult, double th, std::vector<int> &consensusSet)
{
    consensusSet.clear();

    //transform M
    Matrix4X TM(M.cols());
    for(size_t i=0; i<M.cols(); i++)
    {
        TM.x[4*i]   = M.x[3*i];
        TM.x[4*i+1] = M.x[3*i+1];
        TM.x[4*i+2] = M.x[3*i+2];
        TM.x[4*i+3] = 1;
    }

    //TM = guessAndResult * TM;
    TM = transform(TM,guessAndResult);

    for(size_t i=0;i<corr.size();i++){
        int idxS = corr[i].idxS;
        int idxT = corr[i].idxT;
        Vector3 pM(TM(0,idxS),TM(1,idxS),TM(2,idxS));
        Vector3 pB(B(0,idxT),B(1,idxT),B(2,idxT));
        Vector3 pDiff = pM-pB;
        if(pDiff.norm()<=th){
            consensusSet.push_back(i);
            //qDebug()<<"("<<corr.data()[i].idxS<<", "<<corr.data()[i].idxT<<")";
        }
    }

    return consensusSet.size();
}

