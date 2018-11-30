#include "../../../include/globReg4D/include/reg_common.h"
#include "../../../include/globReg4D/include/search.h"
#include "../../../include/globReg4D/include/registration.h"

using namespace reg::search;

int reg::search::rot1_withCorr_IntervalStab(const Matrix3X &X, const Matrix3X &Y,
                                            std::vector<corrTab> &corr,
                                            std::vector<int> &corrNOS, std::vector<int> &corrNOT,
                                            const double &th, const double &tolZ, const int &lwbnd,
                                            const int &gap, AxisAngle &rsearchResult, bool enforceOneToOne)
{
    assert(th>0 && gap>=0 && lwbnd>=0 && "");
    typedef RotationSearchSpaceRegion1DOF RSSR;
    int upbnd;
    const size_t numOfPts = X.cols();
    assert(numOfPts>0 && "");
    Vector TH(numOfPts);
    for (int i=0;i<numOfPts;i++)
    {
        double *m = X.x + 3*i; //.col(matchList[i])
        double sqNormM = m[0]*m[0]+m[1]*m[1]+m[2]*m[2];
        TH.x[i]=th+sqrt(2*(sqNormM-sqNormM*cos(tolZ)));
    }


    ArcEvaluatorWithCorr_IS dsi(X,Y,TH,corr,corrNOS,corrNOT);

    double outAngle;
    dsi.sweep(outAngle,upbnd,enforceOneToOne);
    const Vector3 axis(0,0,1);
    rsearchResult=AxisAngle(axis, outAngle);

    return upbnd;
}
