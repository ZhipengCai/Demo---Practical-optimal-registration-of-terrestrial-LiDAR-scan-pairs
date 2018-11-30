#ifndef REG_ARC_EVALUATOR_
#define REG_ARC_EVALUATOR_

#include "reg_common.h"
#include "data_indexation.h"
#include "state.h"
#include "interval_tree.h"
#include <algorithm>

namespace reg {
    namespace search {

        //!intersection interval of a circular arc and an epsilon-ball, used for doing sweep algorithm
        struct intervalEnd{
            double location; //location of the end point
            bool isStart; //is this end point a starting point of an interval
            int corrIdx;
            void formIntervalEnd(const double &location_in, const bool &isStart_in, const int &corrIdx_in){
                location = location_in;
                isStart = isStart_in;
                corrIdx = corrIdx_in;
            }
        };

        typedef RotationSearchSpaceRegion1DOF RSSR1;
        //!1D rotation evaluator with correspondences via interval stabbing
        class ArcEvaluatorWithCorr_IS
        {
        public:

            ArcEvaluatorWithCorr_IS(const Matrix3X &X, const Matrix3X &Y, const Vector &TH,
                                 const std::vector<corrTab> &corr, const std::vector<int> &corrNOS,
                                 const std::vector<int> &corrNOT);
            ~ArcEvaluatorWithCorr_IS();

            void sweep(double &outAngle, int &outUpbnd, bool oneToOne = false);

            void intervalStab(double &outAngle, int &outUpbnd, bool oneToOne);

            int size() const;

        private:

            const Matrix3X &M_in;
            const Matrix3X &B_in;
            const Vector &TH_in;
            const std::vector<corrTab> &corr_in; //correspondence set
            std::vector<int> matchListS_in;
            std::vector<int> matchListT_in;
            std::vector<intervalEnd> interval_array;

            int _size;
        };

    } // End namespace sarch
} // End namespace reg

#endif
