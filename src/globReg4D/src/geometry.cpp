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

#include "../../../include/globReg4D/include/geometry.h"

double reg::geometry::circleintersection(double R, double d, double r)
{    
    assert(R>=0 && d>=0 && r>=0 && "parametres must be positive");


    //assert(d<(R+r));
    // Return value is between 0 and pi.
    
    double rat, x, angle;
            
    if (d<=DUMMY_PRECISION)
    {
        return PI;
    }

//    if( fabs(d-(R+r))<DUMMY_PRECISION )
//    {
//        return 0;
//    }

    x = (d*d-r*r+R*R)/(2*d);

    rat = x/R;
    if (rat<=-1.0)
    {
        return PI;
    }

    angle= acos(rat);
    assert(angle<=PI && "angle must be < PI");
    return angle;
}



double reg::geometry::circleIntersectionAngle(double R, double d, double r)
{
    assert(R>DUMMY_PRECISION && d>DUMMY_PRECISION && r>DUMMY_PRECISION &&
             "parametres must be > 0");
    // Return value is between 0 and pi.
    // Only intersection if circles touches
    // Return negative if no intersection

    const double x = (d*d-r*r+R*R)/(2.0*d);
    const double rat = x/R;
    if (rat<=-1.0)
    {
        return -1.0;
    }
    return acos(rat);
}


reg::geometry::Circle
reg::geometry::circleSterProj(const Vector3 &p, double alpha)
{
    assert(( sqrt(p.x*p.x+p.y*p.y+p.z*p.z)-1.0)<DUMMY_PRECISION && "norm of p must be =1");
    assert(alpha<.5*PI && "alpha <.5*pi");
    assert(alpha>DUMMY_PRECISION && "alpha must be >0");

    const double psi_p = acos(p.z);
    const double psi_u = psi_p + alpha;
    const double psi_l = psi_p - alpha;

    assert(psi_l<=PI && psi_u>=0 && "inconsistent inclination");

    const double x= p.x;
    const double y= p.y;

    assert( fabs(1-p.z) > DUMMY_PRECISION);

    const double ppnorm = sqrt(x*x+ y*y);
    const double Rl_d = 1. - cos(psi_l);
    const double Ru_d = 1. - cos(psi_u);
    assert(Rl_d>=0);
    assert(Ru_d>=0);

    double Rl = sin(psi_l)/Rl_d;
    double Ru = sin(psi_u)/Ru_d;

    assert((Rl<0 && psi_l<0) || (Rl==0 && fabs(psi_l)<DUMMY_PRECISION) || (Rl>0 && psi_l>0));
    assert((Ru<0 && psi_u>PI) || (Rl==0 && fabs(psi_u-PI)<DUMMY_PRECISION) || (Ru>0 && psi_u<PI));

    Circle c;
    if (ppnorm<DUMMY_PRECISION)
    {
        c.x = 0;
        c.y = 0;
    }
    else
    {
        const double Rd = 0.5*(Rl+Ru); //centre distance from origin
        c.x = (x/ppnorm) * Rd;
        c.y = (y/ppnorm) * Rd;
    }

    c.r = fabs(.5*(Rl-Ru));
    return c;
}


reg::geometry::Circle
reg::geometry::patchSterProj(const Vector3 &p, double angle, bool &pos)
{
    assert((sqrt(p.x*p.x+p.y*p.y+p.z*p.z)-1.0)<DUMMY_PRECISION && "norm of p must be 1");
    assert(angle>DUMMY_PRECISION && angle < PI && "inconsistent angle");

    const double psi_p = acos(p.z);
    const double psi_u = psi_p + angle;
    const double psi_l = psi_p - angle;

    assert(psi_l<=PI && psi_u>=0);
    pos = psi_l>=0;

    const double x = p.x;
    const double y = p.y;

    const double ppnorm = sqrt(x*x + y*y);
    const double Rl_d = 1. - cos(psi_l);
    const double Ru_d = 1. - cos(psi_u);
    assert(Rl_d>=0);
    assert(Ru_d>=0);

    const double Rl = sin(psi_l)/Rl_d;
    const double Ru = sin(psi_u)/Ru_d;

    assert((Rl<0 && psi_l<0)  || (Rl==0 && Rl_d==0) || (Rl>0 && psi_l>0));
    assert((Ru<0 && psi_u>PI) || (Rl==0 && Ru_d==0) || (Ru>0 && psi_u<PI));

    Circle c;
    if (ppnorm<DUMMY_PRECISION)
    {
        c.x = 0;
        c.y = 0;
    }
    else
    {
        const double Rd = 0.5*(Rl+Ru); //centre distance from origin. It can be negative!
        assert((Rd>0&&pos) || (Rd<0&&!pos) || (Rd==0));

        c.x = (x/ppnorm) * Rd;
        c.y = (y/ppnorm) * Rd;
    }

    c.r = fabs(0.5*(Rl-Ru));
    return c;
}


reg::geometry::HalfPlane
reg::geometry::halfPlaneSterProj(const Vector3 &p, double alpha)
{
    assert(alpha>DUMMY_PRECISION && alpha<PI/2. && "inconsistent angle");

    const double inv_d = tan(alpha); // 1/distance
    const double norm_xy = sqrt(p.x*p.x + p.y*p.y); // xy-norm

    assert(inv_d>DUMMY_PRECISION && norm_xy>DUMMY_PRECISION && inv_d*norm_xy>DUMMY_PRECISION);
    HalfPlane hp( p.x/(norm_xy*inv_d),
                  p.y/(norm_xy*inv_d));

    return hp;
}

