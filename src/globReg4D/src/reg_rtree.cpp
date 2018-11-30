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

#include "../../../include/globReg4D/include/reg_rtree.h"

reg::RTree::Node::Node(bool leaf): leaf(leaf), n(0)
{
    if (leaf)
    {
        circle = new Circle[RTREE_M];
    }
}

reg::RTree::Node::~Node()
{
    if (leaf)
    {
        delete []circle;
    }
}

reg::RTree::RTree() : _depth(0), _full(false), _node_path(NULL), _idx_path(NULL)
{
    root = new Node(true);
}

std::ostream& operator<<(std::ostream &os, const reg::RTree *tree)
{
    if (tree!=NULL)
    {
        os<<"no implemented!";
    }
    return os;
}

void reg::RTree::freeRTree(Node *node)
{
    if(!node->leaf)
    {
        for(int i=0; i<(node->n); i++)
        {
            freeRTree(node->child[i]);
        }
    }
    delete node;
}

reg::RTree::~RTree()
{
    freeRTree(root);
    if(_node_path!=NULL)
    {
        delete []_node_path;
        delete []_idx_path;
    }
}


void reg::RTree::addPatch(const Vector3 &centre, double angle)
{
    assert(angle>EPSILON && angle<PI && "inconsistent angle");
    assert((sqrt(centre.x*centre.x+centre.y*centre.y+centre.z*centre.z)-1.0)<=DUMMY_PRECISION &&
             "norm of centre must be = 1");

    const double psi_p = acos(centre.z);
    const double psi_l = psi_p - angle;

    assert(fabs(psi_p + angle)>DUMMY_PRECISION && "inconsistent angle");

    // Line projection
    if ( fabs(psi_l) < DUMMY_PRECISION )
    {
        const bool hp_sign = angle<=.5*PI;
        //std::cout<<"line projection!"<<std::endl;
        if (hp_sign) // Pos
        {
            insertPos( halfPlaneSterProj(centre, angle) );
        }
        else // Neg
        {
            Vector3 negCentre(-centre.x, -centre.y, -centre.z);
            const double negAngle = PI-angle;
            insertNeg( halfPlaneSterProj(negCentre, negAngle) );
        }
    }
    // Circle projection
    else
    {
        bool c_sign;
        Circle c = patchSterProj(centre, angle, c_sign);
        if(c_sign) insertPos(c);
        else       insertNeg(c);
    }
}

void reg::RTree::printNegCircles() const
{
    std::list<Circle>::const_iterator it;
    for (it=negCirclesList.begin(); it!=negCirclesList.end(); ++it)
    {
        Circle c = *it;
        //std::cout<< "circle ("<<c.x<<", "<<c.y<<", "<<c.r<<")"<<std::endl;
    }
}


void reg::RTree::updateMBR(Node **path, int *idx_path)
{
    int idx, i;
    Node *node, *parent;

    for (int k=_depth-1; k>=1; k--)
    {
        node = path[k];
        parent = path[k-1];

        Rectangle r = node->MBR[0];
        for (i=1; i<node->n; i++)
        {
            r = mbr(r, node->MBR[i]);
        }

        // if parent's mbr has not changed, return
        idx = idx_path[k-1];
        Rectangle pmbr = parent->MBR[idx];
        if (r != pmbr)
        {
            parent->MBR[idx]=r;
        }
        else
        {
            return;
        }
    }
}


void reg::RTree::Node::add(Node *node, Rectangle r)
{
    assert(n<RTREE_M && "inconsistent state");
    child[n] = node;
    MBR[n] = r;
    n++;
    assert(n<=RTREE_M && "inconsistent state");
}

void reg::RTree::Node::add(Circle c, Rectangle r)
{
    assert(leaf && "Trying to insert a circle into a no leaf node");
    assert(n<RTREE_M && "inconsistent state");

    MBR[n] = r;
    circle[n] = c;
    n++;

    assert(n<=RTREE_M && "inconsistent state");
}

void reg::RTree::splitInter(Node **path, int *idx_path)
{
    int i,j;
    int idx;
    int seed1, seed2;
    int rem;
    double d, maxd;
    double enl1, enl2;
    Node *node, *parent;

    for (int k=_depth; k-- > 0; )
    {
        node = path[k];

        if (node->n < RTREE_M)
        {
            return;
        }

        if (k==0)
        {
            // create new root
            Node *nr = new Node(false);
            nr->child[0] = node;
            nr->n = 1;
            root = nr;
            _depth++;

            parent = nr;
            idx = 0;
        }
        else
        {
            parent = path[k-1];
            idx = idx_path[k-1];
        }

        maxd = -1;
        seed1 = seed2 = 0;
        for(i=0; i<node->n; i++)
        {
            for (j=i+1; j<node->n; j++)
            {
                d = dist(node->MBR[i], node->MBR[j]);
                if (d>maxd)
                {
                    maxd = d;
                    seed1=i; seed2=j;
                }
            }
        }

        Node *nod1 = new Node(false);
        Node *nod2 = new Node(false);

        Rectangle mbr_nod1 = node->MBR[seed1];
        Rectangle mbr_nod2 = node->MBR[seed2];
        nod1->MBR[0] = mbr_nod1;
        nod2->MBR[0] = mbr_nod2;
        nod1->child[0] = node->child[seed1];
        nod2->child[0] = node->child[seed2];
        nod1->n = nod2->n = 1;

        /*move*/
        node->MBR[seed1] = node->MBR[0];
        node->MBR[seed2] = node->MBR[1];
        node->child[seed1] = node->child[0];
        node->child[seed2] = node->child[1];

        rem = node->n-2;

        for (i=2;i<node->n;i++)
        {
            if (nod1->n+rem == RTREE_m)
            {
                for (; i < node->n; i++)
                {
                    Rectangle r = node->MBR[i];
                    nod1->add(node->child[i], r);
                    mbr_nod1 = mbr(mbr_nod1, r);
                }
                break;
            }

            if (nod2->n+rem == RTREE_m)
            {
                for (; i < node->n; i++)
                {
                    Rectangle r = node->MBR[i];
                    nod2->add(node->child[i], r);
                    mbr_nod2 = mbr(mbr_nod2, r);
                }
                break;
            }

            Rectangle r = node->MBR[i];
            Rectangle mbr_nod1_post = mbr(mbr_nod1, r);
            Rectangle mbr_nod2_post = mbr(mbr_nod2, r);

            enl1 = area(mbr_nod1_post) - area(mbr_nod1);
            enl2 = area(mbr_nod2_post) - area(mbr_nod2);

            if (enl1 < enl2)
            {
                nod1->add(node->child[i], r);
                mbr_nod1 = mbr_nod1_post;
            }
            else if (enl2 < enl1)
            {
                nod2->add(node->child[i], r);
                mbr_nod2 = mbr_nod2_post;
            }
            // equal enlargement area
            else if ( mbr_nod1 < mbr_nod2 ) //compare areas
            {
                nod1->add(node->child[i], r);
                mbr_nod1 = mbr_nod1_post;
            }
            else if ( mbr_nod2 < mbr_nod1 )
            {
                nod2->add(node->child[i], r);
                mbr_nod2 = mbr_nod2_post;
            }
            // equal mbr area
            else if (nod1->n < nod2->n)
            {
                nod1->add(node->child[i], r);
                mbr_nod1 = mbr_nod1_post;
            }
            else
            {
                nod2->add(node->child[i], r);
                mbr_nod2 = mbr_nod2_post;
            }

            rem--;
        }
        delete node;

        parent->child[idx]       = nod1;
        parent->child[parent->n] = nod2;
        parent->MBR[idx]       = mbr_nod1;
        parent->MBR[parent->n] = mbr_nod2;
        parent->n++;
    }

}

void reg::RTree::splitLeaf(Node *leaf, Node **node_path, int *idx_path)
{
    assert(leaf->leaf && "inconsistent state");
    assert(leaf->n == RTREE_M && "inconsistent state");

    int j;
    int seed1, seed2;
    double d, maxd;
    double enl1, enl2;

    maxd = -1;
    seed1 = seed2 = 0;
    for (int i=0; i<leaf->n; i++)
    {
        for (j=i+1; j<leaf->n; j++)
        {
            d = dist(leaf->circle[i], leaf->circle[j]);
            if (d>maxd)
            {
                maxd=d;
                seed1=i; seed2=j;
            }
        }
    }

    Rectangle mbr_l1 = leaf->MBR[seed1];
    Rectangle mbr_l2 = leaf->MBR[seed2];

    Node *l1 = new Node(true);
    Node *l2 = new Node(true);

    l1->add(leaf->circle[seed1], mbr_l1);
    l2->add(leaf->circle[seed2], mbr_l2);

    // move
    leaf->circle[seed1] = leaf->circle[0];
    leaf->circle[seed2] = leaf->circle[1];

    leaf->MBR[seed1] = leaf->MBR[0];
    leaf->MBR[seed2] = leaf->MBR[1];

    leaf->child[seed1]= leaf->child[0];
    leaf->child[seed2]= leaf->child[1];

    int rem = leaf->n-2;
    for (int k=2; k<leaf->n; k++)
    {
        if (l1->n + rem == RTREE_m)
        {
            for (; k<leaf->n; k++)
            {
                Rectangle r = leaf->MBR[k];
                l1->add(leaf->circle[k], r);
                mbr_l1 = mbr(mbr_l1, r);
            }
            break;
        }

        if (l2->n + rem == RTREE_m)
        {
            for (; k<leaf->n; k++)
            {
                Rectangle r = leaf->MBR[k];
                l2->add(leaf->circle[k], r);
                mbr_l2 = mbr(mbr_l2, r);
            }
            break;
        }

        Circle c = leaf->circle[k];
        Rectangle r = leaf->MBR[k];

        Rectangle mbr_l1_post = mbr(mbr_l1, r);
        Rectangle mbr_l2_post = mbr(mbr_l2, r);

        enl1 = area(mbr_l1_post) - area(mbr_l1);
        assert(enl1>=0);

        enl2 = area(mbr_l2_post) - area(mbr_l2);
        assert(enl2>=0);

        if (enl1 < enl2)
        {
            l1->add(c,r);
            mbr_l1 = mbr_l1_post;
        }
        else if (enl2 < enl1)
        {
            l2->add(c,r);
            mbr_l2 = mbr_l2_post;
        }
        // Equal enlargement area
        else if (mbr_l1 < mbr_l2)
        {
            l1->add(c,r);
            mbr_l1 = mbr_l1_post;
        }
        else if (mbr_l2 < mbr_l1)
        {
            l2->add(c,r);
            mbr_l2 = mbr_l2_post;
        }
        // Equal mbr area
        else if (l1->n < l2->n)
        {
            l1->add(c,r);
            mbr_l1 = mbr_l1_post;
        }
        else
        {
            l2->add(c,r);
            mbr_l2 = mbr_l2_post;
        }

        rem--;
    }


    /* add leaf nodes */
    if (_depth==0) /*leaf is the root*/
    {
        // Update root
        root = new Node(false);
        root->add(l1, mbr_l1);
        root->add(l2, mbr_l2);
        _depth=1;
    }
    else
    {
        Node *parent = node_path[_depth-1];
        int idx = idx_path[_depth-1];
        parent->child[idx]      = l1;
        parent->child[parent->n]= l2;
        parent->MBR[idx]        = mbr_l1;
        parent->MBR[parent->n]  = mbr_l2;
        parent->n++;

        // update path nodes
        updateMBR(node_path, idx_path);
        splitInter(node_path, idx_path);
    }
    delete leaf;
}


void reg::RTree::printPosCircles(Node *node) const
{
    if (node->leaf)
    {
        for (int i=0; i<node->n; i++)
        {
            Circle c = node->circle[i];
            //std::cout<< "circle ("<<c.x<<", "<<c.y<<", "<<c.r<<")\n";
        }
    }
    else
    {
        for (int i=0; i<node->n; i++)
        {
            printPosCircles(node->child[i]);
        }
    }
}

void reg::RTree::printPosCircles() const
{
    printPosCircles(root);
}


bool reg::RTree::intersectsPosHalfPlanes(double x, double y) const
{
    std::list<HalfPlane>::const_iterator it;
    for(it=posHalfPlanesList.begin(); it!=posHalfPlanesList.end(); ++it)
    {
        if (contains(*it, x, y))
        {
            return true;
        }
    }
    return false;
}

bool reg::RTree::intersectsNegHalfPlanes(double x, double y) const
{
    std::list<HalfPlane>::const_iterator it;
    for(it=negHalfPlanesList.begin(); it!=negHalfPlanesList.end(); ++it)
    {
        if (!contains(*it, x, y))
        {
            return true;
        }
    }
    return false;
}


bool reg::RTree::intersectsNegHalfPlanes(Circle c) const
{
    std::list<HalfPlane>::const_iterator it;
    for(it=negHalfPlanesList.begin(); it!=negHalfPlanesList.end(); ++it)
    {
        if (!contains(*it, c))
        {
            return true;
        }
    }
    return false;
}


bool reg::RTree::intersectsNegCircles(double x, double y) const
{
    std::list<Circle>::const_iterator it;
    for(it=negCirclesList.begin(); it!=negCirclesList.end(); ++it)
    {
        if (!contains(*it, x, y))
        {
            return true;
        }
    }
    return false;
}

bool reg::RTree::intersectsNegCircles(Circle c) const
{
    std::list<Circle>::const_iterator it;
    for(it=negCirclesList.begin(); it!=negCirclesList.end(); ++it)
    {
        // true if c included in *it
        if (!contains(*it, c))
        {
            return true;
        }
    }
    return false;
}


bool reg::RTree::intersectsPosCircles(Node *node, double x, double y) const
{
    int i;

    if (node->leaf)
    {
        for (i=0; i<node->n; i++)
        {
            if ( contains (node->circle[i], x, y) )
            {
                return true;
            }
        }
        return false;
    }

    for (i=0; i<node->n; i++)
    {
        if ( contains(node->MBR[i], x, y) )
        {
            if ( intersectsPosCircles(node->child[i], x, y) )
            {
                return true;
            }
        }
    }
    return false;
}



bool reg::RTree::intersectsPosCircles(Node *node, Circle c, bool c_sign) const
{
    int i;

    if (node->leaf)
    {
        if(c_sign) //Pos c
        {
            for (i=0; i<node->n; i++)
            {
                if ( intersects(node->circle[i], c ) )
                {
                    return true;
                }
            }
        }
        else // Neg c
        {
            for (i=0; i<node->n; i++)
            {
                if ( !contains(c, node->circle[i] ) )
                {
                    return true;
                }
            }
        }
        return false;
    }

    Rectangle mbrc = mbr(c);
    //Recursion
    if(c_sign) //Pos c
    {

        for (i=0; i<node->n; i++)
        {
            if ( intersects(node->MBR[i], mbrc) )
            {
                if ( intersectsPosCircles(node->child[i], c, c_sign) )
                {
                    return true;
                }
            }
        }
    }
    else // Neg c
    {
        for (i=0; i<node->n; i++)
        {
            if( !intersects(mbrc, node->MBR[i]) )
            {
                return true;
            }

            if( !contains(c, node->MBR[i]) )
            {
                if(intersectsPosCircles(node->child[i], c, c_sign) )
                {
                    return true;
                }
            }
        } //for
    }
    return false;
}



bool reg::RTree::matchPoint(const Vector3 &p) const
{
    if(_full)
    {
        return true;
    }

    if ( fabs( 1. -p.z ) < DUMMY_PRECISION )
    {
        return negCirclesList.size()>0 || posHalfPlanesList.size()>0 || negHalfPlanesList.size()>0;
    }

    // stereographic projection
    const double x = p.x/(1.-p.z);
    const double y = p.y/(1.-p.z);

    if (intersectsPosHalfPlanes(x, y))
    {
        return true;
    }

    if (intersectsNegHalfPlanes(x, y))
    {
        return true;
    }

    if (intersectsNegCircles(x, y))
    {
        return true;
    }

    return intersectsPosCircles(root, x, y);
}


bool reg::RTree::intersectsPosCircles(Node *node, HalfPlane hp, bool pos) const
{
    if (node->leaf)
    {
        if(pos)
        {
            for (int i=0; i<node->n; i++)
            {
                if ( intersects (hp, node->circle[i]) )
                {
                    return true;
                }
            }
        }
        else //Neg HP
        {
            for (int i=0; i<node->n; i++)
            {
                if ( !contains (hp, node->circle[i]) )
                {
                    return true;
                }
            }
        }
        return false;
    }

    //Recursion
    if(pos)
    {
        for (int i=0; i<node->n; i++)
        {
            if( contains(hp, node->MBR[i]))
            {
                return true;
            }
            if ( intersects(node->MBR[i], hp) )
            {
                if ( intersectsPosCircles(node->child[i], hp, pos) )
                {
                    return true;
                }
            }
        }
    }
    else //Neg half plane
    {
        for (int i=0; i<node->n; i++)
        {
            if( !intersects(hp, node->MBR[i]))
            {
                return true;
            }
            if ( !contains(hp, node->MBR[i]) )
            {
                if ( intersectsPosCircles(node->child[i], hp, pos) )
                {
                    return true;
                }
            }
        }
    }
    return false;
}

bool reg::RTree::matchPatch(const Vector3 &patchCentre, double patchAngle) const
{
        if(_full)
        {
            return true;
        }
        if( fabs(patchAngle-PI)<DUMMY_PRECISION)
        {
            return true;
        }

        const double psi_p = acos(patchCentre.z);

        // Half-plane projection
        if ( fabs(psi_p-patchAngle) <= DUMMY_PRECISION )
        {
            //std::cout<<"eval line projection"<<std::endl;

            // If intersect at the pole
            if ( negCirclesList.size()>0 ||
                 posHalfPlanesList.size()>0 ||
                 negHalfPlanesList.size()>0)
            {
                return true;
            }

            //Negative halfPlane
            if(patchAngle<DUMMY_PRECISION)
            {
                // Intersection at the inf. already checked!
                return false;
            }

            if(patchAngle>0.5*PI)
            {
                Vector3 posPatchCentre(
                            -patchCentre.x,
                            -patchCentre.y,
                            -patchCentre.z);

                const double posPatchAngle = PI-patchAngle;

                HalfPlane negHP = halfPlaneSterProj(posPatchCentre, posPatchAngle);
                return intersectsPosCircles(root, negHP, false);
            }
            else //Positive half plane
            {
                HalfPlane posHP = halfPlaneSterProj(patchCentre, patchAngle);
                return intersectsPosCircles(root, posHP, true);
            }
        }

        // Circle projection


        bool circleSign;
        Circle projCircle = patchSterProj(patchCentre, patchAngle, circleSign);

        if (circleSign) //Positive
        {
            if (intersectsNegHalfPlanes(projCircle))
            {
                return true;
            }

            //should intersect pos half planes!

            if (intersectsNegCircles(projCircle))
            {
                return true;
            }

            return intersectsPosCircles(root, projCircle, circleSign);
        }

        // If N is indexed (as inf)
        if(  negCirclesList.size()>0 ||
             posHalfPlanesList.size()>0 ||
             negHalfPlanesList.size()>0 )
        {
            return true;
        }

        // no negative circle in tree
        //return _match_ext_circle(tree->root, c, circle_mbr(c));
        return intersectsPosCircles(root, projCircle, circleSign);
}



void reg::RTree::insertPos(HalfPlane hp)
{
    if(_full) return;
    posHalfPlanesList.push_back(hp);
}

void reg::RTree::insertNeg(HalfPlane hp)
{
    if(_full) return;
    negHalfPlanesList.push_back(hp);
}


void reg::RTree::insertPos(Circle c)
{
    assert(_depth>=0 && "_depth must be >=0");

    if(_full) return;

    int pre_depth=_depth;
    insertPos(root, c, _node_path, _idx_path, 0);
    if (_depth>pre_depth)
    {
        if(_node_path!=NULL)
        {
            delete []_node_path;
            delete []_idx_path;
        }
        _node_path = new Node*[_depth];
        _idx_path = new int[_depth];
    }
}


void reg::RTree::insertPos(Node *node, Circle c,
                           Node **node_path, int *idx_path, int lev)
{
    if (node->leaf)
    {
        // No action if circle is contained in other circle
        for (int i=0; i<node->n; i++)
        {
            if ( contains(node->circle[i], c) )
            {
                return;
            }
        }

        node->circle[node->n] = c;
        node->MBR[node->n] = mbr(c);
        node->n++;

        if (node->n == RTREE_M)
        {
            splitLeaf(node, node_path, idx_path);
        }

    }
    else // Internal node
    {
        // Find child with min enlargement
        double enl, min_enl;
        int idx = 0;
        Rectangle cmbr = mbr(c);
        Rectangle min_mbr = mbr(node->MBR[0], cmbr);
        enl = min_enl = area(min_mbr) - area(node->MBR[0]);

        for (int i=1; enl && i<node->n; i++)
        {
            Rectangle mbr_aux = mbr(node->MBR[i], cmbr);
            enl = area(mbr_aux) - area(node->MBR[i]);
            if (enl<min_enl)
            {
                min_enl = enl;
                min_mbr = mbr_aux;
                idx = i;
            }
        }

        // Update mbr of the selected child
        node->MBR[idx] = min_mbr;

        // Update path
        node_path[lev] = node;
        idx_path[lev] = idx;

        // Traverse
        Node *child = node->child[idx];
        insertPos(child, c, node_path, idx_path, lev+1);
    }
}

void reg::RTree::insertNeg(Circle c)
{
    assert(c.r>0 && "radius must be >0");

    // empty circle
    if (c.r<DUMMY_PRECISION)
    {
        _full = true;
        //depth = 0; //why?
        return;
    }

    // First insertion
    if (negCirclesList.size()==0)
    {
        negCirclesList.push_back(c);
        return;
    }

    // Merge with existent circles
    bool ignore;
    //cn = tree->negc_list;
    ignore = false;

    std::list<Circle>::iterator itr;
    for(itr = negCirclesList.begin(); itr != negCirclesList.end();)
    {
        Circle cn = *itr;
        if(contains(cn, c))   /* if c interior -> delete circle cn->c */
        {
            itr=negCirclesList.erase(itr);
        }
        else
        {
            if(contains(c, cn))
            {
                /*ignore circle if it contains some node circle*/
                ignore=true;
                break;
            }
          ++itr;
        }
    }


    if(!ignore) // Add new circle
    {
        negCirclesList.push_back(c);
    }

}
