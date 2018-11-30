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

#include "../../../include/globReg4D/include/interval_tree.h"
#include <stdlib.h>    // qsort
#include <stdio.h>     // printf
#include "time.h"


//--------------------------------------------
//            Operators
//--------------------------------------------

std::ostream &operator<<(std::ostream& os, const reg::IntervalNode *node)
{
    if (node==NULL)
    {
        return os;
    }

    os << node->left;
    os << "["<< node->lw << " "<< node->up <<"] ";
    os <<node->right;
    return os;
}

std::ostream& operator<<(std::ostream &os, const reg::IntervalTree &t)
{
    os << t.root;
    os << std::endl;
    return os;
}


//--------------------------------------------
//            Interval Node
//--------------------------------------------

using namespace reg;


IntervalNode::~IntervalNode()
{
    delete left;
    delete right;
}

std::ostream& operator<<(std::ostream &os, const IntervalNode& p)
{
    return os<< "[" << p.lw << ", " << p.up << "]\n";
}



//--------------------------------------------
//            Interval Tree
//--------------------------------------------

IntervalTree::IntervalTree(): root(NULL) {}

IntervalTree::~IntervalTree()
{
    delete root;
}


void IntervalTree::insert(IntervalNode *node, float lw, float up)
{
    IntervalNode *parent, *temp;

    // Interval inside of node
    if (node->lw <= lw && up <= node->up)
    {
        return;
    }

    // Insert to the left
    if (up < node->lw )
    {
        if (node->left==NULL)
        {
            node->left = new IntervalNode(lw, up);
        }
        else
        {
            insert(node->left, lw, up);
        }
    }

    // Insert to the right
    else if( lw > node->up )
    {
        if (node->right==NULL)
        {
            node->right = new IntervalNode(lw, up);
        }
        else
        {
            insert(node->right, lw, up);
        }
    }

    // Intersection
    else
    {
        // If [lw up] is equal to [node.lw node.up] no insertion is needed, then
        // only inequalities on boundaries are tested.
        if (lw < node->lw)
        {
            // Update node
            node->lw = lw;

            // Update node->left
            parent = node;
            temp = parent->left;
            while (temp)
            {
                // Remove temp if it is contained into node
                if (temp->lw > lw)
                {
                    if (parent!=node)
                    {
                        parent->right = temp->left;
                        temp->left = NULL;
                        delete temp;
                        temp = parent->right;
                    }
                    else
                    {
                        parent->left = temp->left;
                        temp->left = NULL;
                        delete temp;
                        temp = parent->left;
                    }
                    continue;
                }

                // Stop condition: lw contained in temp
                if ( (temp->lw <= lw) && (lw <= temp->up) )
                {
                    break;
                }

                parent = temp;
                temp = (lw < temp->lw ) ? temp->left : temp->right;
            }// End while

            // temp points to the intersection node
            if (temp != NULL)
            {
                // Update node lw value
                node->lw = temp->lw;
                if (parent != node)
                {
                    parent->right = temp->left;
                }
                else //parent == node
                {
                    node->left = temp->left;
                }
                /* free intersected node */
                temp->left = NULL;
                delete temp;
            }
        }

        if (up > node->up)
        {
            node->up = up;
            /* Find intersection right node*/
            parent = node;
            temp = parent->right;
            while (temp)
            {
                /* remove node */
                if (temp->up < up)
                {
                    if (parent!=node)
                    {
                        parent->left = temp->right;
                        temp->right = NULL;
                        delete temp;
                        temp = parent->left;
                    }
                    else
                    {
                        parent->right = temp->right;
                        temp->right = NULL;
                        delete temp;
                        temp = parent->right;
                    }
                    continue;
                }
                if ( (temp->lw <= up) && (up <= temp->up) )
                {
                    break;
                }

                parent = temp;
                temp = (up > temp->up) ? temp->right : temp->left;
            }

            if (temp != NULL)
            {
                node->up = temp->up;
                if (parent != node)
                {
                    parent->left = temp->right;
                }
                else
                {
                    node->right = temp->right;
                }
                temp->right = NULL;
                delete temp;

            }
        }
    } // Intersection

}



void IntervalTree::insert(float lw, float up)
{
    if (root==NULL)
    {
        root = new IntervalNode(lw, up);
    }
    else
    {
        insert(root, lw, up);
    }

}

bool IntervalTree::matchValue(IntervalNode *node, float val)
{
    if (node==NULL)
    {
        return false;
    }

    if ( (node->lw <= val) && (val <= node->up) )
    {
        return true;
    }

    if (val < node->lw)
    {
        return matchValue(node->left, val);
    }

    return matchValue(node->right, val);
}


bool IntervalTree::matchValue(float val)
{
    return matchValue(root, val);
}

bool IntervalTree::matchInter(IntervalNode *node, float lw, float up)
{
    if(node==NULL)
    {
        return false;
    }

    // Interval inside of node
    if (node->lw <= lw && up <= node->up)
    {
        return true;
    }

    // Match at the left
    if (up < node->lw )
    {
        return matchInter(node->left, lw, up);
    }

    // Match at the right
    else if( lw > node->up)
    {
        return matchInter(node->right, lw, up);
    }

    // Intersection
    return true;
}

bool IntervalTree::matchInter(float lw, float up)
{
    return matchInter(root, lw, up);
}

int IntervalTree::size(IntervalNode *node)
{
    if (node==NULL)
    {
        return 0;
    }
    return size(node->left) + size(node->right) + 1;
}

int IntervalTree::size()
{
    return size(root);
}


int IntervalTree::countPointers(IntervalNode *node)
{
    if (node==NULL)
    {
        return 0;
    }
    return countPointers(node->left)  +
            countPointers(node->right) + 2;
}

int IntervalTree::countPointers()
{
    return countPointers(root);
}








//---------------------------
//        Testing
//---------------------------

typedef IntervalNode Interval;

bool reg::inspect_tree(Interval *node, IntervalNode *interval_list, int *idx)
{
    if (node==NULL)
    {
        return true;
    }

    if(!inspect_tree(node->left, interval_list, idx))
    {
        return false;
    }

    if( node->lw != interval_list[*idx].lw ||
        node->up != interval_list[*idx].up)
    {
        return false;
    }
    (*idx)++;

    if(!inspect_tree(node->right, interval_list, idx))
    {
        return false;
    }

    return true;
}


bool reg::inspect_tree(IntervalTree *tree, IntervalNode *interval_list, int len_list)
{
    bool resp;
    int idx = 0;
    resp =  inspect_tree(tree->root, interval_list, &idx) ;
    resp = (idx==len_list) && resp;
    return resp;
}



int compareInterval(const void * i1, const void * i2)
{
    return ((Interval*)i1)->lw - ((Interval*)i2)->lw;
}

/* return: merged_list size */
int merge(float *a, float *b, int len, Interval **merged_list)
{
    int i, m;
    float top_a, top_b, l, r;
    Interval *interval_list;

    interval_list = new Interval[len]; //(Interval*)malloc(len*sizeof(Interval));
    *merged_list = new Interval[len]; //(Interval*)malloc(len*sizeof(Interval));
    m=0;
    for (i=0;i<len;i++)
    {
        interval_list[i].lw = a[i];
        interval_list[i].up = b[i];
    }

    qsort((void *)(interval_list), len, sizeof(Interval), compareInterval);

    (*merged_list)[0].lw = interval_list[0].lw;
    (*merged_list)[0].up = interval_list[0].up;
    m++;

    for (i=1; i<len; i++)
    {
        l = interval_list[i].lw;
        r = interval_list[i].up;
        top_a = (*merged_list)[m-1].lw;
        top_b = (*merged_list)[m-1].up;

        if (top_b < l)
        {
            (*merged_list)[m].lw = l;
            (*merged_list)[m].up = r;
            m++;
        }
        else if ( top_b < r)
        {
            (*merged_list)[m-1].lw = top_a;
            (*merged_list)[m-1].up = r;
        }
    }
    delete []interval_list;

    return m;
}




int reg::test_IntervalTree()
{
    int i,k;
    int number_of_intervals, m;
    int resp;
    const int TEST_REPETITIONS = 100000;

    float *a, *b;

    IntervalTree *tree;
    clock_t begin, end;
    clock_t dt_tree, dt_ka;
    Interval *merged_list;

    dt_tree = dt_ka = 0;
    resp = 0;

    // Number of random intervals
    number_of_intervals=100;

    // Repite test
    srand (time(NULL));
    for (k=0; k<TEST_REPETITIONS; k++)
    {
        /* generate random intervals */

        a = new float[number_of_intervals];
        b = new float[number_of_intervals];
        for (i=0;i<number_of_intervals;i++)
        {
            a[i] = rand() % 1000 + 1;
            b[i] = a[i] + rand() % 10 + 1;
        }

        /* merge intergals using known algorithm */
        begin = clock();

        m = merge(a, b, number_of_intervals, &merged_list);
        end = clock();
        dt_ka += end-begin;

        /*
        printf("result:\n\n");
        for(i=0; i<m; i++)
        {
            printf("[%f, %f]\n", merged_list[i].lw, merged_list[i].up);
        }
        */

        /* Build and populate interval tree */
        begin = clock();
        tree = new IntervalTree();
        for(i=0; i<number_of_intervals; i++)
        {
           /* printf(" inserting [%f, %f]\n", a[i], b[i]);*/
            tree->insert(a[i], b[i]);
        }
        end = clock();
        dt_tree += end-begin;


        /* check if tree has the same intervals than the known algorithm */
        if (!inspect_tree(tree, merged_list, m))
        {
            printf("ERROR\n");
            resp++;
        }

//        if (k%10==0)
//        {
//            std::cout << *tree <<std::endl;
//        }

        delete[] a;
        delete[] b;
        delete[] merged_list;
        delete  tree;
    }

    end = clock();

    printf("dt known algorithm = %f s\n", (dt_ka)/(double) CLOCKS_PER_SEC);
    printf("dt tree            = %f s\n", (dt_tree)/(double) CLOCKS_PER_SEC);
    if(resp>0)
    {
        printf("Numer of errors = %d\n", resp);
    }
    else
    {
        printf("Test passed without errors\n");
    }

    return resp;
}



