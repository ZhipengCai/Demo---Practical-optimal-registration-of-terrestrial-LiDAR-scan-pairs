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

#include <stdlib.h>
#include <string.h>
#include "reg_common.h"


#ifndef REG_SORT_H
#define REG_SORT_H


namespace reg
{
namespace util
{

//TODO: use a template
struct value_index
{
    double value;
    unsigned int index;
};

struct value_index_uint
{
    unsigned int value;
    unsigned int index;
};


int compare_value_index(const void * a, const void * b);
int compare_value_index_uint(const void * a, const void * b);

unsigned int* sort_index(double *a, unsigned int len);
unsigned int* sort_index(unsigned int *a, unsigned int len);

void sorted_by_index(double *a, unsigned int* idx, unsigned int len);
void sorted_by_index(int *a, unsigned int* idx, unsigned int len);


/*
 * Sort an array according to indexes in idx. An aux array is given.
 */
void sorted_by_index2(double *a, unsigned int* idx, unsigned int len, double *tmp);


} // End namespace util
} // End namespace reg

#endif
