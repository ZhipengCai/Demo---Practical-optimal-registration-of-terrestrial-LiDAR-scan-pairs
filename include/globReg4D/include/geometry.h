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

#ifndef REG_GEOMETRY_
#define REG_GEOMETRY_

#include "reg_common.h"
#include <limits>       // std::numeric_limits
#include <iostream>

namespace reg
{
namespace geometry
{

//------------------
//----- Circle -----
//------------------
struct Circle
{
    double x, y, r;

    Circle(){}

    Circle(double x, double y, double r):
        x(x), y(y), r(r) {}
};

inline
std::ostream& operator<<(std::ostream &os, const Circle &c)
{
    os << "[("<<c.x<<", "<<c.y<<"), "<<c.r<<"])"<<std::endl;
    return os;
}


inline bool operator==(const Circle& lhs, const Circle& rhs)
{
    return  lhs.x==rhs.x && lhs.y==rhs.y && lhs.r==rhs.r;
}
inline bool operator!=(const Circle& lhs, const Circle& rhs){return !operator==(lhs,rhs);}
inline bool operator< (const Circle& lhs, const Circle& rhs)
{
    return lhs.r < rhs.r;
}
inline bool operator> (const Circle& lhs, const Circle& rhs){return  operator< (rhs,lhs);}
inline bool operator<=(const Circle& lhs, const Circle& rhs){return !operator> (lhs,rhs);}
inline bool operator>=(const Circle& lhs, const Circle& rhs){return !operator< (lhs,rhs);}


//designed for 1dof
double circleintersection(double R, double d, double r);

//designed for 3dof
double circleIntersectionAngle(double R, double d, double r);

Circle circleSterProj(const Vector3 &p, double alpha);
Circle patchSterProj(const Vector3 &p, double alpha, bool &circleSign);



//-------------------
//---- HalfPlane ----
//-------------------

struct HalfPlane
{
    double x, y;
    HalfPlane(double x, double y):x(x), y(y){}
};


inline bool operator==(const HalfPlane& lhs, const HalfPlane& rhs)
{
    return  lhs.x==rhs.x && lhs.y==rhs.y;
}

    inline bool operator!=(const HalfPlane& lhs, const HalfPlane& rhs){return !operator==(lhs,rhs);}


inline
std::ostream& operator<<(std::ostream &os, const HalfPlane &h)
{
    os << "[("<<h.x<<", "<<h.y<<")])"<<std::endl;
    return os;
}

HalfPlane halfPlaneSterProj(const Vector3 &p, double alpha);


//-------------------
//---- Rectangle ----
//-------------------

struct Rectangle
{
    double ax, ay, bx, by;
    Rectangle(double ax, double ay, double bx, double by):
        ax(ax), ay(ay), bx(bx), by(by) {}
    Rectangle(){}
    ~Rectangle(){}
};

inline
std::ostream& operator<<(std::ostream &os, const Rectangle &r)
{
    os <<"[("<<r.ax<<", "<<r.ay<<"), ("<<r.bx<<","<<r.by<<")])"<<std::endl;
    return os;
}


inline bool operator==(const Rectangle& lhs, const Rectangle& rhs)
{
    return  lhs.ax == rhs.ax && lhs.ay == rhs.ay &&
            lhs.bx == rhs.bx && lhs.by == rhs.by;
}
inline bool operator!=(const Rectangle& lhs, const Rectangle& rhs){return !operator==(lhs,rhs);}
inline bool operator< (const Rectangle& lhs, const Rectangle& rhs)
{
    return (lhs.bx-lhs.ax)*(lhs.by-lhs.ay) < (rhs.bx-rhs.ax)*(rhs.by-rhs.ay);
}
inline bool operator> (const Rectangle& lhs, const Rectangle& rhs){return  operator< (rhs,lhs);}
inline bool operator<=(const Rectangle& lhs, const Rectangle& rhs){return !operator> (lhs,rhs);}
inline bool operator>=(const Rectangle& lhs, const Rectangle& rhs){return !operator< (lhs,rhs);}


//-------------------
//------ Cube -------
//-------------------

typedef struct Cube
{
    double ax, ay, az, bx, by, bz;
    Cube(double ax, double ay, double az, double bx, double by, double bz):
        ax(ax), ay(ay), az(az), bx(bx), by(by), bz(bz)
    {
        assert(ax<=bx && ay<=by && az<=bz);
    }
    Cube() {}
    ~Cube(){}
} Cube;

inline
std::ostream& operator<<(std::ostream &os, const Cube &c)
{
    os<<"[("<<c.ax<<", "<<c.ay<<", "<<c.az<<"), ("
     <<c.bx<<", "<<c.by<<", "<<c.bz<<")])"<<std::endl;
    return os;
}

inline bool operator==(const Cube& lhs, const Cube& rhs)
{
    return  lhs.ax == rhs.ax && lhs.ay == rhs.ay && lhs.az == rhs.az &&
            lhs.bx == rhs.bx && lhs.by == rhs.by && lhs.bz == rhs.bz;
}
inline bool operator!=(const Cube& lhs, const Cube& rhs){return !operator==(lhs,rhs);}
inline bool operator< (const Cube& lhs, const Cube& rhs)
{
    return (lhs.bx-lhs.ax) < (rhs.bx-rhs.ax);
}
inline bool operator> (const Cube& lhs, const Cube& rhs){return  operator< (rhs,lhs);}
inline bool operator<=(const Cube& lhs, const Cube& rhs){return !operator> (lhs,rhs);}
inline bool operator>=(const Cube& lhs, const Cube& rhs){return !operator< (lhs,rhs);}



//--------------
//  area
//--------------

inline double area(Rectangle r)
{
    return (r.bx-r.ax)*(r.by-r.ay);
}

inline double area(Circle c)
{
    return PI*c.r*c.r;
}


//--------------
//  volume
//--------------

inline double volume(Cube c)
{
    const double x = (c.bx-c.ax);
    return x*x*x;
}



//--------------
//  dist
//--------------
inline
double dist(Rectangle r1, Rectangle r2)
{
    /* find if rectangles intersec */
    if(     (r1.ax <= r2.bx) &&
            (r1.bx >= r2.ax) &&
            (r1.by >= r2.ay) &&
            (r1.ay <= r2.by ) )
    {
        return 0;
    }

    double mindist;
    double d, dx, dy;

    // corners
    double c1_x[4] = {r1.ax, r1.bx, r1.bx, r1.ax};
    double c2_x[4] = {r2.ax, r2.bx, r2.bx, r2.ax};
    double c1_y[4] = {r1.ay, r1.ay, r1.by, r1.by};
    double c2_y[4] = {r2.ay, r2.ay, r2.by, r2.by};

    mindist = std::numeric_limits<double>::max();

    int i, j;
    for(i=0;i<4;i++)
    {
        for(j=0;j<4;j++)
        {
            dx = c1_x[i]-c2_x[j];
            dy = c1_y[i]-c2_y[j];
            d = dx*dx+dy*dy;
            if (d<mindist)
            {
                mindist=d;
            }
        }
    }

    return mindist;
}

inline
double dist(Circle c1, Circle c2)
{
    const double dx = c1.x - c2.x;
    const double dy = c1.y - c2.y;
    const double d = sqrt(dx*dx+dy*dy)-c1.r-c2.r;
    if (d<0)
    {
        return 0;
    }
    return d;
}


//--------------
//     mbr
//--------------

inline
Rectangle mbr(Rectangle r1, Rectangle r2)
{
    Rectangle r( MIN(r1.ax, r2.ax),
                 MIN(r1.ay, r2.ay),
                 MAX(r1.bx, r2.bx),
                 MAX(r1.by, r2.by) );
    return r;
}

inline
Rectangle mbr(Circle c)
{
    Rectangle mbr( c.x-c.r, c.y-c.r,
                   c.x+c.r, c.y+c.r);
    return mbr;
}


//-------------------
//     Contains
//-------------------

inline
bool contains(HalfPlane hp, double x, double y)
{
    return hp.x * x + hp.y * y >= hp.x*hp.x + hp.y*hp.y;
}

inline
bool contains(HalfPlane hp, Rectangle r)
{
    return contains(hp, r.ax, r.ay) &&
           contains(hp, r.bx, r.by) &&
           contains(hp, r.ax, r.by) &&
           contains(hp, r.bx, r.ay);
}

inline
bool contains(HalfPlane hp, Circle c)
{
    //if contains => no intersection between hp line and circle
    //if there is no intersection check if centre point belong to hp positive side
    double m2 = -hp.x/hp.y; m2 *= m2;
    double n2 = (hp.x/hp.y)*(hp.x-c.x) + hp.y-c.y; n2 *= n2;
    double r2 = c.r*c.r;
    if ( m2*n2 <= (m2+1)*(n2-r2)  ) //intersection between hp line and circle
    {
        return false;
    }
    return contains(hp, c.x, c.y);
}

inline
bool contains(Rectangle r1, Rectangle r2)
{
    return r1.ax<=r2.ax && r1.ay<=r2.ay && r1.bx>=r2.bx && r1.by>=r2.by;
}

inline
bool contains(Rectangle r, double x, double y)
{
    return x>=r.ax && r.bx>=x && y>=r.ay && r.by>=y;
}

inline
bool contains(Circle c, double x, double y)
{
    const double dx = c.x - x;
    const double dy = c.y - y;
    if (dx*dx + dy*dy <= c.r * c.r )
    {
        return true;
    }
    return false;
}

inline
bool contains(Circle c1, Circle c2)
{
    const double dx = c1.x - c2.x;
    const double dy = c1.y - c2.y;
    const double d = sqrt(dx*dx + dy*dy);
    return d+c2.r<=c1.r;
}

inline
bool contains(Circle c, Rectangle r)
{
    // true if all r's corners are inside of c
    double dx, dy;
    const double r2 = c.r*c.r;

    dx = c.x - r.ax;
    dy = c.y - r.ay;
    if(dx*dx+dy*dy > r2)
    {
        return false;
    }

    dx = c.x - r.bx;
    dy = c.y - r.by;
    if(dx*dx+dy*dy>r2)
    {
        return false;
    }

    dx = c.x - r.ax;
    dy = c.y - r.by;
    if(dx*dx+dy*dy>r2)
    {
        return false;
    }

    dx = c.x - r.bx;
    dy = c.y - r.ay;
    if(dx*dx+dy*dy>r2)
    {
        return false;
    }

    return true;
}


//-------------------
//     intersects
//-------------------

inline
bool intersects(Rectangle r1, Rectangle r2)
{
    /*
     *A's Left Edge to left of B's right edge  And
    A's right edge to right of B's left edge And
    A's top above B's bottom And
    A's bottom below B's Top*/
    return  (r1.ax <= r2.bx) && (r1.bx >= r2.ax) &&
            (r1.by >= r2.ay) && (r1.ay <= r2.by);

}

inline
bool intersects(HalfPlane hp, Circle c)
{
    //intersection between hp line and circle => intersects
    //if there is no intersection check if centre point belong to hp positive side
    double m2 = -hp.x/hp.y; m2 *= m2;
    double n2 = (hp.x/hp.y)*(hp.x-c.x) + hp.y-c.y; n2 *= n2;
    double r2 = c.r*c.r;
    if ( m2*n2 <= (m2+1)*(n2-r2)  ) //intersection between hp line and circle
    {
        return true;
    }
    return contains(hp, c.x, c.y);

}

inline
bool intersects(Circle c, HalfPlane hp){ return intersects(hp, c);}

inline bool intersects(HalfPlane hp, Rectangle r)
{
    return contains(hp, r.ax, r.ay) ||
           contains(hp, r.bx, r.by) ||
           contains(hp, r.ax, r.by) ||
           contains(hp, r.bx, r.ay);
}

inline bool intersects(Rectangle r, HalfPlane hp){
    return intersects(hp,r);
}

inline
bool intersects(Circle c1, Circle c2)
{
    const double dx = c1.x - c2.x;
    const double dy = c1.y - c2.y;
    const double d2 = dx*dx + dy*dy;
    const double rsum = c1.r+c2.r;
    return d2 <= rsum*rsum;
}


} // End namespace geometry
} // End namespace reg
#endif
