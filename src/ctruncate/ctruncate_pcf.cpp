//
//  pcf.cpp
//  
//
/*  Compute parabolic cylinder function Dv(x)
  
 Translation from routines in program mpbdv.for by Shanjie Zhang and Jianming Jin
 distributed "Computation of Special Functions".
 Copyright 1996 by John Wiley & Sons, Inc.
 Permission has been granted for purchasers of the book to incorporate these
 programs as long as the copyright is acknowledged. */

#include <cmath>
#include <limits>
#include <vector>
#include "clipper/clipper.h"


namespace ctruncate {
                    
    double dvsa( double v, double x)
        //Compute parabolic cylinder function for small argument (<5.8)
        /* v is order
         x is argument
         return Dv(x) */
        {
            const double EPS = 1.0e-15;
            double pd;
            double ep(std::exp(-0.25*x*x) );
            double v0(0.5*(1.0-v) );
            
            if ( v == 0.0 ) {
                pd = ep;
            } else {
                if ( x == 0.0 ) {
                    pd = ( v0 <= 0.0 && std::abs(v0 - round(v0) ) < std::numeric_limits<float>::epsilon() ) ? 0.0 :
                    std::sqrt(clipper::Util::pi() ) * std::pow(2.0,v/2.0) /tgamma(v0);
                } else {
                    double r(1.0), r1(1.0);
                    pd = tgamma(-0.5*v);
                    int m(0);
                    do {
                        r *= -std::sqrt(2.0)*x/double(++m);
                        r1 = r*tgamma(0.5*(double(m)-v) );;
                        pd += r1;
                    } while ( m != 251 && std::abs(r1) >= std::abs(pd)*EPS );
                    pd *= std::pow(2.0,-0.5*v-1.0) * ep/tgamma(-v);
                }
            }
            return pd;
        }
        
        double vvla(double v, double x);
        
        double dvla(double v, double x)
        //compute parabolic cylinder function Dv(x) for large x (>5.8)
        {
            const double EPS = 1.0e-12;
            double ep(std::exp(-0.25*x*x) );
            double a0(std::pow(std::abs(x),v)*ep);
            double r(1.0), pd(1.0);
            int k(1);
            do {
                r *= -0.5*(2.0*double(k)-v-1.0)*(2.0*double(k)-v-2.0)/(double(k)*x*x);
                pd += r;
            } while ( ++k != 17 && std::abs(r/pd) > std::numeric_limits<float>::epsilon() );
            pd *= a0;
            if (x < 0.0 ) {
                pd *= cos(clipper::Util::pi()*v);
                pd += clipper::Util::pi()*vvla(v,-x)/tgamma(-v);
            }
            return pd;
        }

        double vvla(double v, double x)
        // compute parabolic cylinder function Vv(x) for large x
        {
            const double EPS = 1.0e-12;
            double pd(1.0);
            double ep(std::exp(0.25*x*x) );
            double a0(std::pow(std::abs(x),-v-1.0) * std::sqrt(2.0/clipper::Util::pi() ) * ep );
            double r(1.0);
            int k(1);
            do {
                r *= 0.5*(2.0*double(k)+v-1.0)*(2.0*double(k)+v)/(double(k)*x*x);
                pd += r;
            } while ( ++k != 19 && std::abs(r/pd) > std::numeric_limits<float>::epsilon() );
            pd *= a0;
            
            if (x < 0.0 ) {
                pd *= -cos(clipper::Util::pi()*v);
                pd += std::pow(sin(clipper::Util::pi()*v),2)*tgamma(-v)/(clipper::Util::pi()*dvla(v,-x));
            }
            return pd;
        }
                

#ifdef HAVE_LIBGSL
#include "gsl/gsl_sf_hyperg.h"
	
	double Utils::pbdv( double v, double x )
	/*!Compute parabolic cylinder function using libgsl functions
	 Note overflow for x > 20
	 \param v order
	 \param x argument
	 \return Dv(x) */
	{
		double Dv(0.0);
		if (x < 0.0) {
			double F_0 = gsl_sf_hyperg_1F1((v+1.0)/2.0,0.5,-x*x/2.0);
			double F_01 = gsl_sf_hyperg_1F1(1.0+v/2.0,1.5, -x*x/2.0);
			double G_0 = tgamma((1.0-v)/2.0);
			double G_01 = tgamma(-v/2.0);
			Dv = std::exp(0.25*x*x)*std::sqrt(clipper::Util::pi() ) * std::pow(2.0,v/2.0) * (F_0/G_0-std::sqrt(2.0)*x*F_01/G_01);
		} else {
			double U = gsl_sf_hyperg_U(-v/2.0,0.5, x*x/2.0);
			Dv = std::exp(-0.25*x*x)*std::pow(2.0,v/2.0)*U;
		}
		return Dv;
	}

#else
	
    double Utils::pbdv( double v, double x )
        /*!Compute parabolic cylinder function.
		Note overflow for x > 22 on 64-bit
         \param v order
         \param x argument
         \return Dv(x) */
        {
            double pd0,pd1;
            //double v1(v+std::copysign(1.0,v) );
            double v1 = v + ((v >= 0.0) ? 1.0 : -1.0 );
            int nv(v1);
            double v0(v1-double(nv) );
            int na(std::abs(nv) );
            double ep(std::exp(-0.25*x*x) );
            
            std::vector<double> dv(na+1); //
            
            if (v1 >= 0.0) {
                if (v0 == 0.0 ) {
                    pd0 = ep;
                    pd1 = x*ep;
                } else {
                    pd0 = ( std::abs(x) > 5.8 ) ? dvla(v0,x) : dvsa(v0,x);
                    pd1 = ( na >= 1 ) ? ( ( std::abs(x) > 5.8 ) ? dvla(v0+1.0,x) : dvsa(v0+1.0,x) ) : pd0;
                }
                dv[0] = pd0; //
                dv[1] = pd1; //
                for (int k = 2 ; k != na+1 ; ++k ) {
                    double pdf = x*pd1 - (k+v0-1.0)*pd0;
                    dv[k]= pdf; //
                    pd0 = pd1;
                    pd1 = pdf;
                }
            } else {
                if (x <= 0.0 ) {
                    pd0 = (std::abs(x) > 5.8 ) ? dvla(v0,x) : dvsa(v0,x);
                    pd1 = (std::abs(x) > 5.8 ) ? dvla(v0-1.0,x) : dvsa(v0-1.0,x);
                    dv[0] = pd0; //
                    dv[1] = pd1; //
                    for (int k = 2 ; k != na+1 ; ++k ) {
                        double pdf = (-x*pd1+pd0)/(double(k)-v0-1.0);
                        dv[k]= pdf; //
                        pd0 = pd1;
                        pd1 = pdf;
                    }
                } else if ( x <= 2.0 ) {
                    double v2 = (nv == 0 ) ? v1 - 1.0 : v1;
                    int nk(int(-v2) );
                    pd0 = dvsa(v2+1.0,x);
                    pd1 = dvsa(v2,x);
                    dv[nk--] = pd1;
                    dv[nk--] = pd0;
                    for ( int k = nk ; k != -1 ; --k) {
                        double pdf = x*pd0 + (double(k) -v0 +1.0)*pd1;
                        dv[k] = pdf;
                        pd1 = pd0;
                        pd0 = pdf;
                    }
                } else {
                    //double dv0 = (std::abs(x) > 5.8 ) ? dvla(v0,x) : dvsa(v0,x);
                    pd1 = 0.0; pd0 = 1.0e-30;
                    double pdf;
                    for ( int k = (na+100) ; k != -1 ; --k) {
                        pdf = x*pd0 + (double(k)-v0 +1.0)*pd1;
                        if (k <= na ) dv[k] = pdf;
                        pd1 = pd0;
                        pd0 = pdf;
                    }
                    double s0 = ((std::abs(x) > 5.8 ) ? dvla(v0,x) : dvsa(v0,x) )/pdf;
                    for (int k = 0 ; k != na+1 ; ++k ) dv[k] *= s0;
                }
            }
            /*
             std::vector<double> dp(na);
             for (int k = 0 ; k != na ; ++k ) {
             dp(k) = (v1 >= 0.0 ) ? 0.5*x*dv[k]-dv[k+1] : -0.5*x*dv[k] - (std::abs(v0)+double(k))*dv[k+1];
             }
             double pdd(dp[na-1]; */
            return dv[na-1];
        }
#endif //HAVE_LIBGSL
}

