#ifndef REG_INTERVAL_TREE_
#define REG_INTERVAL_TREE_

#include <iostream>


namespace reg
{
class IntervalTree;

typedef struct IntervalNode
{
    float lw;
    float up;
    struct IntervalNode *left;
    struct IntervalNode *right;

    IntervalNode(): lw(0), up(0), left(NULL), right(NULL) {}
    IntervalNode(float lw, float up): lw(lw), up(up), left(NULL), right(NULL) {}
    ~IntervalNode();
} IntervalNode;

} // End namespace reg


//-------------------------------------
//        Operators
//-------------------------------------

std::ostream& operator<<(std::ostream &, const reg::IntervalNode *);
std::ostream& operator<<(std::ostream &, const reg::IntervalTree &);



namespace reg
{
class IntervalTree
{
private:
    IntervalNode *root;
    void insert(IntervalNode *node, float lw, float up);

    bool matchValue(IntervalNode *node, float val);
    bool matchInter(IntervalNode *node, float lw, float up);
    int size(IntervalNode *node);
    int countPointers(IntervalNode *node);


public:
    IntervalTree();
    ~IntervalTree();
    void insert(float lw, float up);

    bool matchValue(float val);
    bool matchInter(float lw, float up);
    int size();
    int countPointers();


    //friend std::ostream &printIntervalTree(std::ostream &os, const IntervalNode &node);
    friend std::ostream& (::operator <<) (std::ostream &os, const reg::IntervalTree &);

    friend bool inspect_tree(IntervalNode *node, IntervalNode *interval_list, int *idx);
    friend bool inspect_tree(IntervalTree *tree, IntervalNode *interval_list, int len_list);

 };



int test_IntervalTree();

} // End namespace reg
#endif
