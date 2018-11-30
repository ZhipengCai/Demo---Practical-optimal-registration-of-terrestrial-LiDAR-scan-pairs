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


#ifndef REG_STATE_
#define REG_STATE_

#include <iostream>
#include "reg_common.h"

namespace reg {
    namespace search {
        
        //----------------------------------------------------
        //       Search Space Regions (SSR)
        //----------------------------------------------------
        
        struct RotationSearchSpaceRegion1DOF
        {
            double lw; // Lower angle
            double up; // Upper angle
            
            RotationSearchSpaceRegion1DOF() {lw=0;up=TWOPI;}
            RotationSearchSpaceRegion1DOF(double lw, double up):
            lw(lw), up(up) {}
        };
        
        inline
        std::ostream &operator<<(std::ostream& os, const RotationSearchSpaceRegion1DOF &reg)
        {
            os << "["<< reg.lw <<" "<< reg.up <<"]";
            return os;
        }
    
        struct TranslationSearchSpaceRegion3DOF
        {
            double ax, ay, az;
            double bx, by, bz;
            
            TranslationSearchSpaceRegion3DOF() {}
            TranslationSearchSpaceRegion3DOF(double ax, double ay, double az,
                                             double bx, double by, double bz):
            ax(ax), ay(ay), az(az),
            bx(bx), by(by), bz(bz) {}
            
            TranslationSearchSpaceRegion3DOF(double *x):
            ax(x[0]), ay(x[1]), az(x[2]),
            bx(x[3]), by(x[4]), bz(x[5]) {}
        };
        
        
        inline
        std::ostream &operator<<(std::ostream& os, const TranslationSearchSpaceRegion3DOF &reg)
        {
            os << "[("<< reg.ax <<" "<< reg.ay <<" "<< reg.az <<")"
            << " ("<< reg.bx <<" "<< reg.by <<" "<< reg.bz <<")]";
            return os;
        }

        inline
        double centre(RotationSearchSpaceRegion1DOF ssr)
        {
            return .5*(ssr.lw+ssr.up);
        }
        
        
        inline
        double ssrMaxAngle(RotationSearchSpaceRegion1DOF ssr)
        {
            // Compute max. rotation angle
            double r = .5*(ssr.up-ssr.lw);
            return  r;
        }
        
        
        inline
        int split(RotationSearchSpaceRegion1DOF ssr,
                  RotationSearchSpaceRegion1DOF **partition)
        {


            if (fabs(ssr.lw - ssr.up)<ROT_DIV_ACC)
                return 0;

            const double centre = .5*(ssr.lw + ssr.up);

            partition[0] = new RotationSearchSpaceRegion1DOF(ssr.lw, centre);
            partition[1] = new RotationSearchSpaceRegion1DOF(centre, ssr.up);

            return 2;
        }
        

        inline
        int split(TranslationSearchSpaceRegion3DOF ssr,
                  TranslationSearchSpaceRegion3DOF **partition)
        {
            const double ax = ssr.ax;
            const double ay = ssr.ay;
            const double az = ssr.az;

            const double bx = ssr.bx;
            const double by = ssr.by;
            const double bz = ssr.bz;

            if(bx-ax < TRANS_DIV_ACC)
            {
                //std::cout<<"warning! "<<std::endl;
                return 0;
            }

            const double mx = .5*(ax + bx);
            const double my = .5*(ay + by);
            const double mz = .5*(az + bz);

            partition[0] = new TranslationSearchSpaceRegion3DOF(ax, ay, az, mx, my, mz);
            partition[1] = new TranslationSearchSpaceRegion3DOF(mx, ay, az, bx, my, mz);
            partition[2] = new TranslationSearchSpaceRegion3DOF(ax, my, az, mx, by, mz);
            partition[3] = new TranslationSearchSpaceRegion3DOF(mx, my, az, bx, by, mz);
            partition[4] = new TranslationSearchSpaceRegion3DOF(ax, ay, mz, mx, my, bz);
            partition[5] = new TranslationSearchSpaceRegion3DOF(mx, ay, mz, bx, my, bz);
            partition[6] = new TranslationSearchSpaceRegion3DOF(ax, my, mz, mx, by, bz);
            partition[7] = new TranslationSearchSpaceRegion3DOF(mx, my, mz, bx, by, bz);

            return 8;
        }


        inline
        Vector3 ssrCentre(TranslationSearchSpaceRegion3DOF ssr)
        {
            Vector3 c (.5*(ssr.ax + ssr.bx),
                       .5*(ssr.ay + ssr.by),
                       .5*(ssr.az + ssr.bz));
            return c;
        }

        inline
        double ssrUncertainty(TranslationSearchSpaceRegion3DOF ssr)
        {
            /* compute max. translation uncertainty */
            const double dx = ssr.ax - ssr.bx;
            const double dy = ssr.ay - ssr.by;
            const double dz = ssr.az - ssr.bz;
            return .5*sqrt(dx*dx + dy*dy + dz*dz);
        }
        
        inline bool
        box_intersects_phi_sphere(double ax, double ay, double az,
                                  double bx, double by, double bz)
        {
            double cx, cy, cz;
            const double pi2 = PI*PI;
            
            cx = .5*(ax + bx); cx *= cx;
            cy = .5*(ay + by); cy *= cy;
            cz = .5*(az + bz); cz *= cz;
            
            if (cx + cy + cz <= pi2 ) // centre
            {
                return true;
            }
            
            ax *= ax; ay *= ay; az *= az;
            if (ax + ay + az <= pi2 ) // a
            {
                return true;
            }
            
            bx *= bx; by *= by; bz *= bz;
            
            //check remaining 7 corners
            return  bx + by + bz <= pi2 ||
            ax + by + az <= pi2 ||
            bx + by + az <= pi2 ||
            bx + ay + az <= pi2 ||
            bx + ay + bz <= pi2 ||
            ax + ay + bz <= pi2 ||
            ax + by + bz <= pi2;
        }
        
        
        
        //----------------------------------------------------
        //       Search State
        //----------------------------------------------------
        
        template<class SSR, typename Scalar=int>
        class SearchState
        {
        public:
            SSR ssr;
            Scalar bnd;
            SearchState(SSR ssr, Scalar bnd): ssr(ssr), bnd(bnd) {}
            ~SearchState(){}
            
            friend std::ostream &operator<<(std::ostream& os, const SearchState &ss)
            {
                os<< "ssr "<<ss.ssr<<" "<<" bnd "<<ss.bnd ;
                return os;
            }

        };
        
        
        template<class SSR, typename Scalar=int>
        class SearchStateMatches
        {
        public:
            SSR ssr;
            Scalar bnd;
            Vector3 tx;
            Vector3 ty;
            SearchStateMatches(SSR ssr, Scalar bnd, Vector3 tx, Vector3 ty):
            ssr(ssr), bnd(bnd), tx(tx), ty(ty) {}
            ~SearchStateMatches(){}
            
            friend std::ostream &operator<<(std::ostream& os, const SearchStateMatches &ss)
            {
                os<< "ssr "<<ss.ssr<<" "<<" bnd "<<ss.bnd ;
                return os;
            }

        };
        
        
        template<class SSR, typename Scalar=int>
        class SearchStateML
        {
        public:
            SSR ssr;
            Scalar bnd;
            int *matchList;
            int matchListSize;
            
            SearchStateML(SSR ssr, Scalar bnd,
                          int *matchList, int matchListSize):
            ssr(ssr), bnd(bnd), matchList(matchList),
            matchListSize(matchListSize)
            {
                assert(bnd>=0 && matchListSize>=0 /*&& matchListSize<=bnd*/ && "invalid bound value");
            }
            
            friend std::ostream &operator<<(std::ostream& os, const SearchStateML &ss)
            {
                os<< "ssr "<<ss.ssr<<" "<<" bnd "<<ss.bnd <<" ml size "<<ss.matchListSize;
                return os;
            }


            ~SearchStateML()
            {
                delete []matchList;
                matchList=NULL;
            }
        };
        
        
    } // End namespace sarch
} // End namespace reg

#endif
