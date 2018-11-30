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

#ifndef REG_BINARYTREE_H
#define REG_BINARYTREE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "math.h"


namespace reg
{
namespace binarytree {


/* Definition for interval.*/
typedef struct interval
{
    double lw;
    double up;
} interval;

/* Definitions for binary search tree.*/
typedef struct payload
{
    double val;
    int order;
} payload;

typedef struct treeNode
{
    payload data;
    struct treeNode *left;
    struct treeNode *right;
} treeNode;

treeNode *Insert(treeNode*, payload, treeNode*);
void free_Binarytree(treeNode *node);
int queryLower(treeNode*,double,treeNode*);
int queryUpper(treeNode*,double,treeNode*);
double queryMiddle(treeNode *,double,treeNode *);
void PrintInorder(treeNode*);
int size_Binarytree(treeNode*);
int count_pointers_Binarytree(treeNode*);

} // End namespace binarytree
} // End namespace reg

#endif
