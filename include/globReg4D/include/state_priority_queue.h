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


#ifndef REG_STATE_PRIORITY_QUEUE_H_
#define REG_STATE_PRIORITY_QUEUE_H_

#include "state.h"
#include <cstddef> //NULL

namespace reg {
namespace search {


template <class SSR, typename Scalar=int >
class StatePriorityQueue
{
public:
    enum OptProblem{MINIMISATION, MAXIMISATION};

private:

    class Node
    {
    public:
        SearchState<SSR, Scalar> *state;
        Node *left, *right;

        Node(): state(NULL), left(NULL), right(NULL) {}
        Node(SearchState<SSR, Scalar> *state):state(state), left(NULL), right(NULL){}
        ~Node() {if (state!=NULL) delete state;}
    };

    const OptProblem optProblem;
    Node *head, *tail;
    unsigned int m_size;

public:
    StatePriorityQueue(OptProblem op=MAXIMISATION);
    ~StatePriorityQueue();

    SearchState<SSR, Scalar> *pop();
    void push(SearchState<SSR, Scalar> *state);

    /**
     * @brief Remove and free states with upper bound lower or equal to lwbnd.
     * @param lwbnd Known lower bound.
     */
    void prune(int curbest);

    unsigned int size() const;
};


} // End namespace search
} // End namespace reg

#include "state_priority_queue.hpp"

#endif
