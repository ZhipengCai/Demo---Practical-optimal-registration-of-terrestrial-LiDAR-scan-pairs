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

#ifndef REG_REG_SEARCH_
#define REG_REG_SEARCH_

#include "reg_common.h"
#include "state.h"
#include "rot1_evaluator.h"

namespace reg {
namespace search {

Matrix4X transform( Matrix4X &x, Transform3 &tform);

//--------------------------------------------------------------------------
//     1D Rotation search
//--------------------------------------------------------------------------

int rot1_withCorr_IntervalStab(const Matrix3X &X, const Matrix3X &Y, std::vector<corrTab> &corr,
                               std::vector<int> &corrNOS, std::vector<int> &corrNOT,
                               const double &th, const double &tolZ,
                               const int &lwbnd, const int &gap, AxisAngle &rsearchResult,
                               bool enforceOneToOne = false);

//--------------------------------------------------------------------------
//     Registration
//--------------------------------------------------------------------------
int nestedbnb_search_4dof_withCorr_IS(const Matrix3X &M, const Matrix3X &B,
                                        std::vector<corrTab> &corr, std::vector<int> &corrNOS,
                                        std::vector<int> &corrNOT, Transform3 &guessAndResult,
                                        double th, double tolZ, int gap, int lwbnd,
                                        TranslationSearchSpaceRegion3DOF &tr_region, bool enforceOneToOne = false);

int computeConsensus_withCorr(const Matrix3X &M, const Matrix3X &B,
                              std::vector<corrTab> &corr, Transform3 &guessAndResult,
                              double th, std::vector<int> &consensusSet);
} // End namespace sarch
} // End namespace reg

#endif

