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


#ifndef REG_SEARCH_
#define REG_SEARCH_

#include "reg_common.h"
#include "data_indexation.h"

namespace reg {
    namespace search {
        
        /**
         * @brief Find optimal quality of a transformation
         * @param psi Indexation of point sets.
         * @param knownSolutionQual Known solution quality.
         * @param gap BnB stop gap.
         * @param guessAndResult Guess search region and final region such that
         *        upbnd-lwbnd<=gap
         * @return Quality of the central transform of the optimal region.
         */
        template <class SSR, unsigned int BRANCHING_FACTOR>
        int bnb_search_queue(const DataIndexation<SSR> &psi,
                             int lwbnd, int gap, SSR &guessAndResult);
        
        
        template <class SSR, unsigned int BRANCHING_FACTOR>
        int bnb_search_table(const DataIndexation<SSR> &dsi,
                             int lwbnd, int gap, int buckets,
                             SSR &guessAndResult);
        
        template <class SSR, unsigned int BRANCHING_FACTOR>
        int searchTableDF(const DataIndexation<SSR> &dsi,
                          int lwbnd, int gap, int buckets,
                          SSR &guessAndResult);
        
        
        template <class SSR, unsigned int BRANCHING_FACTOR>
        int bnb_search_ml_table(const DataIndexation<SSR> &dsi,
                                int lwbnd, int gap, int buckets,
                                SSR &guessAndResult);
        
        
    } // End namespace sarch
} // End namespace reg

#include "search.hpp"

#endif

