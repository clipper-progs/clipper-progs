//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#include "ctruncate_twin.h"
#include "ctruncate_utils.h"
#include "ccp4/csymlib.h"
#include "twinlaws.h"
#include <string.h>
#include <cmath>
#include <iomanip>
#include <cstdlib>

namespace ctruncate {
	
    /*! Summary comments using alpha from H-test and L statistic from L-test.
     \param alpha Alpha from the H-test
     \param lval L statistic from L-test. */
    void twin_summary(clipper::ftype alpha, clipper::ftype lval) {
        printf("TWINNING SUMMARY\n\n");
        printf("Twinning fraction from H-test: %6.2f\n",alpha);
        printf("L-statistic from L-Test:       %6.2f\n\n",lval);
		printf("   Relation between L statistics and twinning fraction:\n");
		printf("      Twinning fraction = 0.000  L statistics = 0.500:\n");
		printf("      Twinning fraction = 0.100  L statistics = 0.440:\n");
		printf("      Twinning fraction = 0.500  L statistics = 0.375:\n\n");
		
        if ( alpha <= 0.0001 ) {
            if ( 0.440 <= lval && lval <= 0.500 ) {
                printf("NO Twinning detected\n\n");
            } else if ( 0.375 <= lval && lval < 0.440 ) {
                printf("L-test statistics indicate partial twinning\n");
                printf("   It is quite likely that your data were merged into a HIGHER symmetry space group than the\n   true space group.\n");
                printf("   Please revise the space group assignment if there are problems with model building or refinement.\n");
                printf("   (Very week data in higher resolution shell may be a reason of this L-value. Run twinning tests\n  with resolution cut off 3A.)\n\n");
			}
        } else if ( 0.0001 < alpha && alpha <= 0.1 ) {
            printf("No twinning or very low twinning fraction.\n\n");
            printf("   Twinning, if any, can be safely be ignored. However, twin refinement may be attempted, but not before the\n model is completely build.\n\n");
        } else if ( 0.1 < alpha && alpha < 0.4 ) {
            printf("It is highly probable that your crystal is TWINNED.\n\n");
            printf("   Please use twin refinement after your model is almost completed and R-free is below 40%%.\n\n");
        } else {
            if ( 0.440 <= lval && lval <= 0.500 ) {
                printf("Your data might have been scaled in a LOWER symmetry space group than the true space group.\n\n");
                printf("   Please run pointless, aimless (scala) and ctruncate to revise space group assignment.\n\n");
            } else if ( 0.375 <= lval && lval < 0.440 ) {
                printf("It is highly probable that your crystal is TWINNED.\n\n");
                printf("   Please use twin refinement after your model is almost completed and R-free is below 40%%\n\n");
            }
        }
    }
    
    
    //-------Twin Isymops-------------------------------------------------------
    
	/*! class to hold twinning test results
	 */
	
	/*! output name of operator
	 */
	std::string TwinSymops::Twinop::description() const {
		clipper::String s;
		clipper::Mat33<int> mat33=this->rot();
		MatrixToString(mat33,s);
		s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
		return s;
	 }
	 
    /*! constructor for list of twinning operators.
     This uses either first principles.
     \param cell Cell of system
     \param spg Spacegroup of system
     \return type TwinSymops */
    TwinSymops::TwinSymops(const clipper::Cell& cell, const clipper::Spacegroup& spg) {
		_cell = cell;
		_spg = spg;
        principles(cell, spg);
    }
    
	/*! operator() to perform calculation
	 \param cell Cell of system
     \param spg Spacegroup of system
	 \return in number of symops */
	int TwinSymops::operator()(const clipper::Cell& cell, const clipper::Spacegroup& spg) {
		_cell = cell;
		_spg = spg;
        principles(cell, spg);
		return size();
    }
		
    /*! output name of operator
     */
    std::string TwinSymops::description(const int& i) const {
        clipper::String s;
        clipper::Mat33<int> mat33=_twinops[i].rot();
        MatrixToString(mat33,s);
		s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
        return s;
    }
	
	/*! output type of operator
	 */
	std::string TwinSymops::Twinop::type(const clipper::Spacegroup& spg) const {
		bool m(false);
		//std::cout << "type: " << _sc << " " << spg.symbol_laue();
		if (_sc <= 0.0001 ) { // low obligueness
			std::string s(spg.symbol_laue() );
			if ( s[0] == '3' || s[0] == '4' || s[0] == '6' || ( s[0] == 'm' && s[1] == '3' ) ) { // trigonal, tetragonal, hexagonal or cubic
				clipper::Mat33<int> tmp = rot(); 
				clipper::Rotation rot(clipper::Mat33<clipper::ftype>(
																	 clipper::ftype(tmp(0,0) )/12.0,
																	 clipper::ftype(tmp(0,1) )/12.0,
																	 clipper::ftype(tmp(0,2) )/12.0,
																	 clipper::ftype(tmp(1,0) )/12.0,
																	 clipper::ftype(tmp(1,1) )/12.0,
																	 clipper::ftype(tmp(1,2) )/12.0,
																	 clipper::ftype(tmp(2,0) )/12.0,
																	 clipper::ftype(tmp(2,1) )/12.0,
																	 clipper::ftype(tmp(2,2) )/12.0));
				if (std::fabs(rot.w() - clipper::Util::pi() ) <= 0.0001 ) m = true; // two-fold
				//std::cout << " " << std::fabs(rot.w() - clipper::Util::pi() );
			}
		}
		//std::cout << std::endl;
		return ( m ? "m" : "pm" );
	}
	
    /*! private function to generate tabulated twinops
    */
    /*void TwinSymops::table( const clipper::Cell& cell, const clipper::Spacegroup& spg) {
        
        clipper::Spacegroup pg(clipper::Spgr_descr(spg.generator_ops().pgrp_ops() ) );
        int sg = spg.spacegroup_number();
        
		clipper::Mat33<clipper::ftype> twinop(0,0,0,0,0,0,0,0,0);
		int scalefac = 1;
		clipper::String s;
		if ( (sg >= 75 && sg <= 80) || sg == 146 || (sg >= 168 && sg <= 173) || (sg >= 195 && sg <= 199) ) { 
			s = "k, h, -l";
			twinop(0,1) = 12;
			twinop(1,0) = 12;
			twinop(2,2) = -12;
            //_twinops.push_back(clipper::Isymop(clipper::RTop<int>(twinop)) );
            _twinops.push_back(clipper::Isymop(clipper::Symop(clipper::RTop_frac(twinop)), clipper::Grid(12,12,12) ) );
		}
		else if( sg >= 149 && sg <= 154 ) {
			s = "-h, -k, l";
			twinop(0,0) = -12;
			twinop(1,1) = -12;
			twinop(2,2) = 12;
            //_twinops.push_back(clipper::Isymop(clipper::RTop<int>(twinop)) );
            _twinops.push_back(clipper::Isymop(clipper::Symop(clipper::RTop_frac(twinop)), clipper::Grid(12,12,12) ) );
		}
		else if( sg >= 143 && sg <= 145 ) {
			s = "k, h, -l";
			twinop(0,1) = 12;
			twinop(1,0) = 12;
			twinop(2,2) = -12;
            //_twinops.push_back(clipper::Isymop(clipper::RTop<int>(twinop)) );
            _twinops.push_back(clipper::Isymop(clipper::Symop(clipper::RTop_frac(twinop)), clipper::Grid(12,12,12) ) );
			
			s = "-k, -h, -l";
			twinop(0,1) = -12;
			twinop(1,0) = -12;
            //_twinops.push_back(clipper::Isymop(clipper::RTop<int>(twinop)) );
            _twinops.push_back(clipper::Isymop(clipper::Symop(clipper::RTop_frac(twinop)), clipper::Grid(12,12,12) ) );
			
			s = "-h, -k, l";
			twinop(0,1) = 0;
			twinop(1,0) = 0;
			twinop(0,0) = -12;
			twinop(1,1) = -12;
			twinop(2,2) = 12;
            //_twinops.push_back(clipper::Isymop(clipper::RTop<int>(twinop)) );
            _twinops.push_back(clipper::Isymop(clipper::Symop(clipper::RTop_frac(twinop)), clipper::Grid(12,12,12) ) );
            
        }
		else if( pg.spacegroup_number() == 16 ) { //PG222 or "P 2 2"
			// Can have pseudo-merohedral twinning in PG222 (orthorhombic) if a=b, b=c or c=a
			if ( fabs( 1.0 - cell.b()/cell.a() ) < 0.02 ) { 
				s = "k, h, -l";
				twinop(0,1) = 12;
				twinop(1,0) = 12;
				twinop(2,2) = -12;
                //_twinops.push_back(clipper::Isymop(clipper::RTop<int>(twinop)) );
                _twinops.push_back(clipper::Isymop(clipper::Symop(clipper::RTop_frac(twinop)), clipper::Grid(12,12,12) ) );
 
			}
			if ( fabs( 1.0 - cell.c()/cell.b() ) < 0.02 ) { 
				s = "-h, l, k";
				twinop(0,1) = 0;
				twinop(1,0) = 0;
				twinop(2,2) = 0;
				twinop(0,0) = -12;
				twinop(1,2) = 12;
				twinop(2,1) = 12;
                //_twinops.push_back(clipper::Isymop(clipper::RTop<int>(twinop)) );
                _twinops.push_back(clipper::Isymop(clipper::Symop(clipper::RTop_frac(twinop)), clipper::Grid(12,12,12) ) );
 
			}
			if ( fabs( 1.0 - cell.a()/cell.c() ) < 0.02 ) {
				s = "l, -k, h";
				twinop(0,0) = 0;
				twinop(1,2) = 0;
				twinop(2,1) = 0;
				twinop(1,1) = -12;
				twinop(0,2) = 12;
				twinop(2,0) = 12;
                //_twinops.push_back(clipper::Isymop(clipper::RTop<int>(twinop)) );
                _twinops.push_back(clipper::Isymop(clipper::Symop(clipper::RTop_frac(twinop)), clipper::Grid(12,12,12) ) );

			}
		}
		else if( pg.spacegroup_number() == 3 ) { //PG2 or "P 2y"
			// can have pseudomerohedral twinning in PG2 (monoclinic) if
			// beta = 90
			// a=c
			// cos(beta) = -a/2c, -c/a, -c/2a
			// sin(beta) = a/c
			if ( fabs( 1.0 - cell.a()/cell.c() ) < 0.02  || fabs( sin(cell.beta()) - cell.a()/cell.c() ) < 0.02 ) {
				s = "l, -k, h";
				twinop(1,1) = -12;
				twinop(0,2) = 12;
				twinop(2,0) = 12;
                //_twinops.push_back(clipper::Isymop(clipper::RTop<int>(twinop)) );
                _twinops.push_back(clipper::Isymop(clipper::Symop(clipper::RTop_frac(twinop)), clipper::Grid(12,12,12) ) );

			}
			
			if ( cell.beta() < clipper::Util::d2rad(93.0) ) {
				s = "-h, -k, l";
				twinop(0,0) = -12;
				twinop(0,2) = 0;
				twinop(1,1) = -12;
				twinop(2,0) = 0;
				twinop(2,2) = 12;
                //_twinops.push_back(clipper::Isymop(clipper::RTop<int>(twinop)) );
                _twinops.push_back(clipper::Isymop(clipper::Symop(clipper::RTop_frac(twinop)), clipper::Grid(12,12,12) ) );

			}	 
			
			if ( fabs( cos(cell.beta()) + 0.5*cell.a()/cell.c() ) < 0.02 ) {
				s = "-h, -k, h+l";
				twinop(0,0) = -12;
				twinop(0,2) = 0;
				twinop(1,1) = -12;
				twinop(2,0) = 12;
				twinop(2,2) = 12;
                //_twinops.push_back(clipper::Isymop(clipper::RTop<int>(twinop)) );
                _twinops.push_back(clipper::Isymop(clipper::Symop(clipper::RTop_frac(twinop)), clipper::Grid(12,12,12) ) );

			}	  
			
			if ( fabs( cos(cell.beta()) + 0.5*cell.c()/cell.a() ) < 0.02 ) { 
				s = "h+l, -k, -l";
				twinop(0,0) = 12;
				twinop(0,2) = 12;
				twinop(1,1) = -12;
				twinop(2,0) = 0;
				twinop(2,2) = -12;
                //_twinops.push_back(clipper::Isymop(clipper::RTop<int>(twinop)) );
                _twinops.push_back(clipper::Isymop(clipper::Symop(clipper::RTop_frac(twinop)), clipper::Grid(12,12,12) ) );

			}	  
			if ( fabs( cos(cell.beta()) + cell.c()/cell.a() ) < 0.02 ) { 
				s = "h+2l, -k, -l";
				twinop(0,0) = 12;
				twinop(0,2) = 24;
				twinop(1,1) = -12;
				twinop(2,0) = 0;
				twinop(2,2) = -12;
                //_twinops.push_back(clipper::Isymop(clipper::RTop<int>(twinop)) );
                _twinops.push_back(clipper::Isymop(clipper::Symop(clipper::RTop_frac(twinop)), clipper::Grid(12,12,12) ) );
  
			}
        }
        //_scalefac = scalefac;
        return ;
    } */

    /*! compute twinning operators using fp
     */
    void TwinSymops::principles(const clipper::Cell& cell, const clipper::Spacegroup& spgr) {		
        const int scalefac = 12;           // scale factor for integer symops
        double sc_tol = 3.5;    // tolerance for determining candidate twinops, obliquity in degrees
        int lc = 48;             // maximum number of twinops that can be stored
        int nc = 0;
        int ivb = 0;
        int ierr = 0;
        int nc2;
        
        std::vector< clipper::Vec3<int> > trans;
        std::vector< clipper::Mat33<int> > rot;
        //std::vector< clipper::Mat33<int> > twin(lc);   // NB need to dimension these o/w output from fortran corrupted
        std::vector<double> score(lc);
        
        int nsymops = spgr.num_symops();
        //printf("nsymops = %d\n",nsymops);
        clipper::Grid g( 12, 12, 12 );
        for (int i=0; i!=nsymops; ++i) {
            clipper::Isymop isymop( spgr.symop(i), g );
            /*for (int j=0; j<3; j++) {
             printf("%6d %6d %6d   %6d\n", isymop.rot()(j,0), isymop.rot()(j,1), isymop.rot()(j,2), isymop.trn()[j] );
             }
             printf("\n");*/
            trans.push_back( isymop.trn() );
            rot.push_back( isymop.rot() );
        }
        
        int vv[nsymops][3][3];
        int ww[nsymops][3];
        int uu_c[lc][3][3];
        double *sc_c = &score[0];
        double twin_cell[6];
        { 
            twin_cell[0] = cell.a();
            twin_cell[1] = cell.b();
            twin_cell[2] = cell.c();
            twin_cell[3] = cell.alpha();
            twin_cell[4] = cell.beta();
            twin_cell[5] = cell.gamma();
            for (int i=0; i!=nsymops; ++i) for (int j=0 ; j != 3; ++j ) for (int k = 0 ; k != 3 ; ++k ) vv[i][j][k] = rot[i](k,j);
            for (int i=0; i!=nsymops; ++i) for (int j=0 ; j != 3; ++j ) ww[i][j] = trans[i][j];
        }
        
        yyy_cell2tg(twin_cell, sc_tol, nsymops, vv, ww, lc, nc, nc2, uu_c, sc_c, ivb, ierr);
                    
        if (nc > 1) {
            for (int k=1; k!=nc; ++k) {
                // Transpose matrix once because passed from Fortran to C, then a second time because Andrey's
                // convention is that (h,k,l) is postmultiplied by the twinning operator, whereas mine is that
                // (h,k,l) is premultiplied by the the twinning operator. Net result: don't do anything!
                clipper::Mat33<int> twinop = clipper::Mat33<int>(uu_c[k][0][0],uu_c[k][0][1],uu_c[k][0][2],uu_c[k][1][0],uu_c[k][1][1]
                                                                 ,uu_c[k][1][2],uu_c[k][2][0],uu_c[k][2][1],uu_c[k][2][2]);
                _twinops.push_back(Twinop(clipper::RTop<int>(twinop),sc_c[k]) );
                //_scalefac = scalefac;
            }
        }
    }
    
    /*void TwinSymops::summary() const {
        if (_mode == FP) {
        printf("First principles calculation of potential twinning operators using code by Andrey Lebedev:\n");
        if (!size() ) 
            printf("First principles calculation has found no potential twinning operators\n\n");
        else
            printf("First principles calculation has found %d potential twinning operators\n\n", size() );
            /*clipper::String s;
            MatrixToString(twinoper,s);
            std::cout << "Twinning operator: " << s << std::endl;
            std::cout << "      Obliquity of the twin operator is " << sc_c[k] << " degrees" << std::endl;
            std::cout << "      (above 2 degrees probably be ignored)" << std::endl;*/

        /*} else {
            printf("The tablulated operators give %d potential twinning operators\n\n", size() );
        }
    }*/
	
	void TwinSymops::output() const {
		if (!size() ) 
			printf("First principles calculation has found no potential twinning operators\n\n");
		else {
			printf("First principles calculation has found %d potential twinning operators\n\n", size() );
			std::cout << "   # twinning operator   score  type" << std::endl;
			for (int i= 0 ; i != size() ; ++i)
				std::cout << std::right << std::setw(4) << i << " " << std::left << std::setw(17) << (_twinops[i]).description() << "   " << std::fixed << std::setw(4) << std::setprecision(2) << (_twinops[i]).score() << "   " << (_twinops[i]).type(_spg) << std::endl;
			std::cout << " m  merohedral\n pm pseudo-merohedral" << std::endl;
			std::cout << "The score gives an indication of the closure of the twinning operation.  The lower the values\n the more higher the overlap." << std::endl;
		}
	}
	
	
	std::stringstream& TwinSymops::xml_output(std::stringstream& ss)
	{
        ss << "<TwinOps>" << std::endl;
        if (size() != 0 ) {
            for (int i= 0 ; i != size() ; ++i) {
                ss << "<Operator number=\"" << i << "\">"  << std::endl;
                ss << "  <Description>" << (_twinops[i]).description() << "</Description>" << std::endl;
                ss << "  <Score>" << std::fixed << std::setw(4) << std::setprecision(2) << (_twinops[i]).score() << "</Score>" << std::endl;
                ss << "  <Type>" << (_twinops[i]).type(_spg) << "</Type>" << std::endl;
                ss << "</Operator>" << std::endl;
            }
            ss << "<Comment id='TwinOps'> m  merohedral\n pm pseudo-merohedral" << std::endl;
            ss << "The score gives an indication of the closure of the twinning operation.  The lower the values\n the more higher the overlap.</Comment>" << std::endl;
        } else {
            ss << "<Comment id='TwinOps'>No potential twinning operators found.</Comment>" << std::endl;
        }
		ss << "</TwinOps>" << std::endl;
		return ss;
	}

    //-------L-test-------------------------------------------------------------
    
    /*! constructor for the L-test.
     \param isig Experimental intensities
     \param reso Resolution range for calculation
     \param nbins Number of bins for the CDF (default 20)
     \return type L_test */
    /*template<class T, template<class> class D> L_test::L_test( clipper::HKL_data< D<T> >& isig, clipper::Range<clipper::ftype> reso, int nbins) {
        _reso = reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
		_alpha = -1.0;
        _cdf.resize(nbins);
		std::vector<clipper::Symop> tncs;
        calc(isig,tncs);
        return;
    }

	template<class T, template<class> class D> L_test::L_test( clipper::HKL_data< D<T> >& isig, std::vector<clipper::Symop>& tncs, clipper::Range<clipper::ftype> reso, int nbins) {
        _reso = reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
		_alpha = -1.0;
        _cdf.resize(nbins);
		
        calc(isig,tncs);
        return;
    } */
	

    /*! recalculate L-test.
     \param reso Resolution limit for calculation (default 3.0A)
     \return L-statistic */
    /*template<class T, template<class> class D> clipper::ftype L_test::operator() (clipper::HKL_data< D<T> >& isig, clipper::Range<clipper::ftype> reso ) {
        _reso = reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
	_alpha = -1.0;
		std::vector<clipper::Symop> tncs;
		
        return calc(isig,tncs);
    }
        
	template<class T, template<class> class D> clipper::ftype L_test::operator() (clipper::HKL_data< D<T> >& isig, std::vector<clipper::Symop>& tncs, clipper::Range<clipper::ftype> reso ) {
        _reso = reso;
		if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
		_alpha = -1.0;
		
        return calc(isig,tncs);
    }*/

	
	/*! return alpha estimate from <|L|>
	 \return alpha */
	clipper::ftype L_test::estimateAlphafromL() const {
		const clipper::ftype PERFECT_TWIN = 0.375;
		const clipper::ftype NO_TWIN = 0.5;
		clipper::ftype alpha;
		
		if (_Lav >= NO_TWIN ) {
			alpha = 0.0;
		} else if ( _Lav <= PERFECT_TWIN ) {
			alpha = 0.5;
		} else {
			clipper::ftype x1(0.01), x2(0.49), x(0.50), y(0.0);
			clipper::ftype y1 = meanL(x1);
			clipper::ftype y2 = meanL(x2);
			
			do {
				x = interpolate(_Lav, y1, y2, x1, x2);
				y = meanL(x);
				if ( (_Lav - y)/(y2-y1) > 0.0 ) {
					y1 = y;
					x1 = x;
					x2 = (x2+x)*0.5;
					y2 = meanL(x2);
				} else {
					y2 = y;
					x2 = x;
					x1 = (x1+x)*0.5;
					y1 = meanL(x1);
				}
			} while ( std::abs(y - _Lav) > 0.001 );
			alpha = x;
		}
		return alpha;
	}
	
	/*! return alpha estimate from <|L*L|>
	 \return alpha */
	clipper::ftype L_test::estimateAlphafromL2() const {
		const clipper::ftype PERFECT_TWIN = 0.2;
		const clipper::ftype NO_TWIN = 0.333333;
		clipper::ftype alpha;

		if (_L2av >= NO_TWIN ) {
			alpha = 0.0;
		} else if ( _L2av <= PERFECT_TWIN ) {
			alpha = 0.5;
		} else {
			clipper::ftype x1(0.01), x2(0.49), x(0.50), y(0.0);
			clipper::ftype y1 = meanL2(x1);
			clipper::ftype y2 = meanL2(x2);
			
			do {
				x = interpolate(_L2av, y1, y2, x1, x2);
				y = meanL2(x);
				if ( (_L2av - y)/(y2-y1) > 0.0 ) {
					y1 = y;
					x1 = x;
					x2 = (x2+x)*0.5;
					y2 = meanL2(x2);
				} else {
					y2 = y;
					x2 = x;
					x1 = (x1+x)*0.5;
					y1 = meanL2(x1);
				}
			} while ( std::abs(y - _L2av) > 0.001 );
			alpha = x;
		}
		return alpha;		
	}
	
	/*! return alpha estimate from N(L)
	 \return alpha */
	clipper::ftype L_test::estimateAlphafromNL() const {
		const clipper::ftype UPPER = 1.0;
		const clipper::ftype LOWER = 0.0;
		clipper::ftype alpha(0.0);
		clipper::ftype na(0.0);
		
		for (int i=0;i != _cdf.size() ;++i) {
			clipper::ftype l = (double(i+1))/double(_cdf.size());
			clipper::ftype nl = double(_cdf[i])/_NLT;
						
			clipper::ftype x1(0.01), x2(0.49), x(0.50), y(0.0);
			clipper::ftype y1 = meanNL(x1,l);
			clipper::ftype y2 = meanNL(x2,l);
			if (nl < y1) {
				alpha += 0.0;
			} else if ( nl > y2 ) {
				alpha += 0.5;
			} else {
				do {
					x = interpolate(nl, y1, y2, x1, x2);
					y = meanNL(x,l);
					if ( (nl - y)/(y2-y1) > 0.0 ) {
						y1 = y;
						x1 = x;
					} else {
						y2 = y;
						x2 = x;
					}
				} while ( std::abs(y - nl) > 0.0001 );
				alpha += x;
			}
			na += 1.0;
		}
		return alpha/na;
	}
	
    /*! print summary from L-test
     */
    void  L_test::summary() const {
        printf("\nL test for twinning: (Padilla and Yeates Acta Cryst. D59 1124 (2003))\n");
        printf("L statistic = %6.3f  (untwinned 0.5 perfect twin 0.375)\n", _Lav);
        printf("Data has used to %6.2f - %6.2f A resolution\n",1.0/std::sqrt(_reso.min()), 1.0/std::sqrt(_reso.max()) );
		printf("   Relation between L statistics and twinning fraction:\n");
        printf("      Twinning fraction = 0.000  L statistics = 0.500:\n");
        printf("      Twinning fraction = 0.100  L statistics = 0.440:\n");
        printf("      Twinning fraction = 0.500  L statistics = 0.375:\n");
		
        
        printf("\n");
		if (_Lav < 0.44) {
			printf("The L test suggests data is twinned\n");			
            printf("All data regardless of I/sigma(I) has been included in the L test\n");
		}
        printf("\n\n");
    }
    
    /*! print loggraph of cdf
     */
    void L_test::loggraph() const {
        printf("$TABLE: L test for twinning:\n");
		printf("$GRAPHS");
		printf(": cumulative distribution function for |L|, twin fraction of %4.2f:0|1x0|1:1,2,3,4:\n",fraction() );
		printf("$$ |L|   N(L) Untwinned Twinned $$\n$$\n");
		printf("0.0000 0.0000 0.0000 0.0000\n");
		
		for (int i=0;i != _cdf.size() ;++i) {
			double x = (double(i+1))/double(_cdf.size());
			printf("%6.4f %6.4f %6.4f %6.4f\n", x, double(_cdf[i])/_NLT, x, 0.5*x*(3.0-x*x)  );
		}
		printf("$$\n\n");
    }
	
	
	std::stringstream& L_test::xml_output(std::stringstream& ss)  {
		ss << "  <L-test>" << std::endl;
		ss << "     <Twinned>" << ((_Lav < 0.44) ? "Yes" : "No") << "</Twinned>" << std::endl;
		ss << "     <Result>" << _Lav << "</Result>" << std::endl;
		ss << "     <TwinFraction>" << fraction() << "</TwinFraction>" << std::endl;
		ss << "     <Comment id=\"L-test\">Relation between L statistics and twinning fraction:\n  Twinning fraction = 0.000  L statistics = 0.500:\n  Twinning fraction = 0.100  L statistics = 0.440:\n  Twinning fraction = 0.500  L statistics = 0.375:</Comment>" << std::endl;
		ss << "  </L-test>" << std::endl;
		return ss;
	}
	
	std::stringstream& L_test::xml_graph(std::stringstream& ss) const {
		ss << "<CCP4Table groupID=\"graph\" id=\"L-test\" title=\"L-test for twinning\">" << std::endl;
		ss << "<plot>" << std::endl;
		ss << "<title>cumulative distribution function for |L|, twin fraction of "<< std::setprecision(2) << fraction() <<"</title>" << std::endl;
		ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"4\" >" <<  std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>red</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"3\" >" << std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>green</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"  1\" ycol=\"  2\" >" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>black</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "</plot>" << std::endl;
		ss << "<headers separator=\" \">\n |L| N(L) Untwinned Twinned \n </headers>" << std::endl;
		ss << "<data>" << std::endl;;
		ss << " 0.0000 0.0000 0.0000 0.0000" << std::endl;
		for (int i=0;i != _cdf.size() ;++i) {
			double x = (double(i+1))/double(_cdf.size());
			ss << std::fixed << std::setprecision(4) << " " << x << " " << double(_cdf[i])/_NLT << " " <<  x << " " <<  0.5*x*(3.0-x*x) << std::endl;
		}		
		ss << "</data>" << std::endl;
		ss << "</CCP4Table>" << std::endl;
		
		return ss;
	}
    
    template<class T, template<class> class D> clipper::ftype L_test::calc(clipper::HKL_data< D<T> >& isig, std::vector<clipper::Symop>& tncs) {
        typedef clipper::HKL_data_base::HKL_reference_index HRI;

		clipper::Spacegroup spgr = isig.hkl_info().spacegroup();
        clipper::Cell cell = isig.hkl_info().cell();
		std::vector<clipper::HKL> steps;
		
		for ( int g = 1 ; g != 7 ; g += 1) {
			for ( int delta1 = -g; delta1 <= g; delta1 += 1 ) {
				for ( int delta2 = -g; delta2 <= g; delta2 += 1 ) {
					for ( int delta3 = -g; delta3 <= g; delta3 += 1 ) {
						if ( !(delta1==0 && delta2==0 && delta3==0) && (std::abs(delta1)==g || std::abs(delta2)==g || std::abs(delta3)==g) ) {
							clipper::HKL hkl2,hkl3;
							hkl2.h() = delta1;
							hkl2.k() = delta2;
							hkl2.l() = delta3;
							bool sys_abs(false);
							for ( int i = 1; i != spgr.num_symops(); ++i ) {
								hkl3 = hkl2.transform(spgr.symop(i));
								clipper::ftype shift = hkl2.sym_phase_shift(spgr.symop(i));
								if ( hkl3 == hkl2 ) {
									if ( cos(shift) < 0.999 ) {
										sys_abs = true; // flag sysabs
										break;
									}
								}
							}
							for (std::vector<clipper::Symop>::const_iterator i = tncs.begin() ; i != tncs.end() ; ++i ) {
								clipper::Coord_frac v1( i->trn());
								for ( int j = 0; j != spgr.num_symops(); ++j ) {
									clipper::Coord_frac v2 = v1.transform(spgr.symop(j) );
									//clipper::ftype shift = hkl2.sym_phase_shift(clipper::Symop(clipper::RTop_frac(clipper::Mat33<clipper::ftype>::identity(), v2) ));
									clipper::ftype shift = -clipper::Util::twopi()*(clipper::ftype(hkl2.h() )*v2[0]+clipper::ftype(hkl2.k() )*v2[1]+clipper::ftype(hkl2.l() )*v2[2]);
									if ( cos(shift) < 0.9 ) {
										sys_abs = true;
										break;
									}
								}
							}
							if (!sys_abs) {
								steps.push_back(hkl2);
							}
						}
					}
				}
			}
			if (steps.size() >= 9) break;
		}
		
		//got nothing so ignore the tNCS (should flag this)
		if (steps.size() == 0 ) {
			for ( int g = 1 ; g != 7 ; g += 1) {
				for ( int delta1 = -g; delta1 <= g; delta1 += 1 ) {
					for ( int delta2 = -g; delta2 <= g; delta2 += 1 ) {
						for ( int delta3 = -g; delta3 <= g; delta3 += 1 ) {
							if ( !(delta1==0 && delta2==0 && delta3==0) && (std::abs(delta1)==g || std::abs(delta2)==g || std::abs(delta3)==g) ) {
								clipper::HKL hkl2,hkl3;
								hkl2.h() = delta1;
								hkl2.k() = delta2;
								hkl2.l() = delta3;
								bool sys_abs(false);
								for ( int i = 1; i != spgr.num_symops(); ++i ) {
									hkl3 = hkl2.transform(spgr.symop(i));
									clipper::ftype shift = hkl2.sym_phase_shift(spgr.symop(i));
									if ( hkl3 == hkl2 ) {
										if ( cos(shift) < 0.999 ) {
											sys_abs = true; // flag sysabs
											break;
										}
									}
								}
								if (!sys_abs) {
									steps.push_back(hkl2);
								}
							}
						}
					}
				}
				if (steps.size() >= 9) break;
			}
			
		}
		
        double LT=0.0;
		double LT2=0.0;
		double NLT=0.0;
        
        for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() && !ih.hkl_class().centric() ) {
				if (_reso.contains(ih.invresolsq() ) ) {
					clipper::HKL hkl = ih.hkl();
					clipper::ftype weight = 2.0/ih.hkl_class().epsilonc();
					int h = hkl.h();
					int k = hkl.k();
					int l = hkl.l();
					for (std::vector<clipper::HKL>::const_iterator i = steps.begin() ; i != steps.end() ; ++i ) {
						clipper::HKL hkl2 = hkl + (*i);
						clipper::ftype I1 = isig[ih].I()/ih.hkl_class().epsilon();
						if (!isig[hkl2].missing() ) {
							clipper::ftype I2 = isig[hkl2].I()/spgr.hkl_class(hkl2).epsilon();
							clipper::ftype L(0.0);
							if ( I1 != 0.0 && I2 != 0.0 ) L = (I2-I1)/(I2+I1);
							if (fabs(L) < 1.0){
								LT += std::fabs(L)*weight;
								LT2 += L*L*weight;
								NLT += weight;
								for (int i=0;i != _cdf.size() ;++i) {
									if ( std::fabs(L) < (clipper::ftype(i+1))/double(_cdf.size()) ) _cdf[i]+=weight;
								}
							}
						}
					}
				}
			}
		}
			
        _NLT = NLT;
		_Lav = LT/NLT;
		_L2av = LT2/NLT;
		_alpha = std::max(std::max(this->estimateAlphafromL(),this->estimateAlphafromL2() ),this->estimateAlphafromNL() );
        return _Lav;
    }
	
     
    //-------H-test-------------------------------------------------------------
    
    /*! constructor for the H-test.
     \param isig Experimental intensities
     \param twinop Twinning operator to investigate
     \param reso Resolution limit for calculation (default 3.0A)
     \param nbins Number of bins for the CDF (default 20)
     \return type H_test */
    /*template<class T, template<class> class D> H_test::H_test( const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twin, clipper::Range<clipper::ftype> reso , int nbins) {
        _nbins = nbins;
        _reso = reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
        _twinop = twin;
        calc(isig,twin);
    }*/
    
    /*! recalculate H-test.
     \param reso Resolution limit for calculation (default 3.0A)
     \return alphas */
    /*template<class T, template<class> class D> const clipper::ftype H_test::operator() (const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twin, clipper::Range<clipper::ftype> reso ) { 
        _reso=reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
        _twinop=twin;
        calc(isig,twin);
        return _alpha;
    }*/
        
    /*! output name of operator
     */
    const std::string H_test::description() const {
        clipper::String s;
        clipper::Mat33<int> mat33=_twinop.rot();
        MatrixToString(mat33,s);
        return s;
    }
    
    /*! print summary from H-tests
     */
    void H_test::summary() const {
        printf("\nH test for twinning: (Yeates Acta Cryst. A44 142 (1980))\n");
        std::cout << "Twinning Operator " << description() << " with H " << std::setprecision(2) <<  _Hav;
        if (_alpha > 0.05) {
            printf(", fraction = %5.2f\n\n",_alpha );
        } else {
            std::cout << ", NO twinning" << std::endl << std::endl;
        }
    }
    
    /*! print loggraphs from H-tests
     */
    void H_test::loggraph() const {
        printf("\n");
        if (_alpha > 0.05) {
            printf("$TABLE: H test for twinning (operator %s) alpha = %5.2f:\n", description().c_str(), _alpha );
            printf("$GRAPHS");
            printf(": cumulative distribution function for |H|:0|1x0|1:1,2,3,4,5,6,7:\n");
            printf("$$ |H| 0.4 0.3 0.2 0.1 0.0 Observed $$\n$$\n");
            printf("0.000000 0.0 0.0 0.0 0.0 0.0 0.000000\n");
            for (int j=0;j!=_cdf.size();++j) {
                double x = (double(j+1))/(double((_cdf).size()));
                printf("%f  -   -   -   -   -  %f\n", x, _cdf[j]);
            }
            printf("1.000000 5.0 2.5 1.667 1.25 1.0 1.0\n");
            printf("$$\n\n");
        }
        
    }
    
    template<class T, template<class> class D>  void H_test::calc(const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twinop) {
        typedef clipper::HKL_data_base::HKL_reference_index HRI;
        
		clipper::Spacegroup spgr = isig.hkl_info().spacegroup();
		
        double HT(0.0), HT2(0.0), NT(0.0);
        double scalefac = 12.0;
        
        clipper::Vec3<int> jhkl, jhkl2;
        for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
            if ( !isig[ih].missing() && !ih.hkl_class().centric() ) {
				if (_reso.contains(ih.invresolsq() ) ) {
					clipper::HKL hkl = ih.hkl();
					jhkl[0] = hkl.h();
					jhkl[1] = hkl.k();
					jhkl[2] = hkl.l();
					clipper::HKL twin;
					jhkl2 = twinop*jhkl;
					twin.h() = jhkl2[0]/scalefac;
					twin.k() = jhkl2[1]/scalefac;
					twin.l() = jhkl2[2]/scalefac;
					if (!isig[twin].missing()) {
						double I1 = isig[ih].I()/ih.hkl_class().epsilon();
						double I2 = isig[twin].I()/spgr.hkl_class(twin).epsilon();
						double weight = 2.0/ih.hkl_class().epsilonc();
						double H = 0.0;
						if ( I1 != 0.0 && I2 != 0.0) H = (I2-I1)/(I2+I1);
						if (fabs(H) < 1){
							HT += fabs(H)*weight;
							HT2 += H*H*weight;
							NT += weight;
							for (int j=0; j != _cdf.size(); ++j) {
								if ( fabs(H) < (double(j+1))/_cdf.size() ) _cdf[j] += 1.0;
							}
						}
					}
                }
            }
		}
		for (int i = 0; i != _cdf.size(); ++i ) _cdf[i] /= _cdf.back();
		_Hav = HT/NT;
		_H2av = HT2/NT;
		_alpha = std::max(0.5*(1.0 - sqrt(3.0*_H2av)),0.0);
    }
    
    
    //-------Britton-test-------------------------------------------------------
    
    /*! constructor for the Britton-test.
     \param isig Experimental intensities
     \param twinop Twinning operator
     \param reso Resolution limit for calculation (default 3.0A)
     \param nbins Number of bins for the CDF (default 50)
     \return type Britton_test */
    /*template <class T> Britton_test::Britton_test( const clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, const clipper::Isymop& twinop, clipper::Range<clipper::ftype> reso , int nbins) {
        _reso=reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
        _twinop=twinop;
        calc(isig,twinop);        
    }*/
    
    /*! recalculate Britton-test.
     \param reso Resolution limit for calculation (default 3.0A)
     \return alphas */
    /*template <class T> const clipper::ftype Britton_test::operator() (const clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, const clipper::Isymop& twinop, clipper::Range<clipper::ftype> reso ) {
        _reso=reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
        _twinop=twinop;
        calc(isig,twinop);
        return _alpha;
    }*/
    
    /*! output name of operator
     */
    const std::string Britton_test::description() const {
        clipper::String s;
        clipper::Mat33<int> mat33=_twinop.rot();
        MatrixToString(mat33,s);
        return s;
    }

    /*! print summary from Britton-tests
     */
    void Britton_test::summary() const {
        printf("\nBritton/Murray-Rust test for twinning (Britton Acta Cryst. A28 296 (1972)): \n");
        std::cout << "Twinning Operator " << description();
        if (_alpha > 0.05) {
            printf(", fraction = %5.2f\n\n",_alpha );
        } else {
            std::cout << ", NO twinning" << std::endl << std::endl;
        }
    }
    
    /*! print loggraphs from Britton-tests
     */
    void Britton_test::loggraph() const {
        printf("\n");
        printf("$TABLE: Britton test for twinning (operator %s) alpha = %5.2f:\n", description().c_str(), _alpha );
        printf("$GRAPHS");
        printf(": Murray-Rust plot |Britton|:A:1,3,4:\n");
        printf(": Britton plot |Britton|:A:1,2:\n");
        printf("$$ alpha Britton Murray-Rust fit $$\n$$\n");
        
        double a(1.0/_beta), b(_alpha/_beta);
        for (int j=0;j!=_pdf.size();++j) {
            double x = (double(j))/(2.0*double((_pdf.size())));
            double y = a*x - b;
            printf("%f %f %f %f\n", x, double(_pdf[j]), double(_zpdf[j]), (y > 0 ) ? y : 0 );
        }
        printf("$$\n\n");
        //}
        
    }
    
    template<class T, template<class> class D>  void Britton_test::calc(const clipper::HKL_data<D<T> >& isig, const clipper::Isymop& twinop) {
        
        typedef clipper::HKL_data_base::HKL_reference_index HRI;
        
		clipper::Spacegroup spgr = isig.hkl_info().spacegroup();
		
        clipper::ftype scalefac(12.0);
        
        int nbins(_pdf.size() );
        
        std::vector<clipper::ftype> x(nbins);
        for (int j = 0; j != nbins ; ++j) x[j] = double(j)/double(2.0*nbins);
        std::vector<clipper::ftype> weights(nbins,1.0);
        // Britton test
        
		clipper::ftype NT(0.0);
        clipper::Vec3<int> jhkl, jhkl2;
        for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
            if ( !isig[ih].missing() && !ih.hkl_class().centric() ) {
				if (_reso.contains(ih.invresolsq() ) ) {
					clipper::HKL hkl = ih.hkl();
					jhkl[0] = hkl.h();
					jhkl[1] = hkl.k();
					jhkl[2] = hkl.l();
					clipper::HKL twin;
					jhkl2 = twinop*jhkl;
					twin.h() = jhkl2[0];
					twin.k() = jhkl2[1];
					twin.l() = jhkl2[2];
					twin.h() = jhkl2[0]/scalefac;
					twin.k() = jhkl2[1]/scalefac;
					twin.l() = jhkl2[2]/scalefac;
					//printf("%d %d %d %d %d %d\n",hkl.h(),hkl.k(),hkl.l(),twin.h(),twin.k(),twin.l());
					//printf("%d %d %d %d %d %d\n",jhkl[0],jhkl[1],jhkl[2],jhkl2[0],jhkl2[1],jhkl2[2]);
					if (!isig[twin].missing()) {
						double I1 = isig[ih].I()/ih.hkl_class().epsilon();
						double I2 = isig[twin].I()/spgr.hkl_class(twin).epsilon();
						NT += 2.0/ih.hkl_class().epsilonc();
						
						if ( I1 > 0.0 && I2 > 0.0 && I2 < I1 ) {
							double B = I2/(I1+I2);
							int bin = int(2*nbins*B);
							if (bin >= 0 && bin < nbins) _pdf[bin]++;
						}
						for (int j = 0; j != nbins ; ++j) {
							clipper::ftype alpha = 0.5*double(j)/double(nbins);
							clipper::ftype J = ((1.0-alpha)*I1-alpha*I2)/(1.0-2*alpha);
							if ( J < 0.0 ) _zpdf[j]++;
						}
					}
                }
            }
        }
        for (int j = 0; j != _pdf.size() ; ++j) _pdf[j] /= NT;
        for (int j = 0; j != _zpdf.size() ; ++j) _zpdf[j] /= NT;
        NT = 0.0;
        for (int j = 0; j != _zpdf.size()-1 ; ++j) {
            weights[j] = ( _zpdf[j] == 0.0 && _zpdf[j+1] == 0.0 ) ? 0.0 : 1.0 ;
            NT += 1.0;
        }

        clipper::ftype alpha(-1.0),alpha1(1.0),beta(0.0);
        do {
            clipper::ftype a, b, siga, sigb;
            straight_line_fit(_zpdf,x,weights,nbins,a,b,siga,sigb);
            if (b >= 0.5) {
                alpha = 0.5;
				beta = 0.0;
				break;
            } else if (b < 0.0 ) {
				alpha = 0.0;
				beta = 1.0;
				break;
			}
            alpha1 = alpha;
            alpha = b;
            beta = a;
            int limit = int((b)*2.0*float(nbins));
            for (int jj = 0 ; jj != limit ; ++jj) weights[jj] = 0.0;
            for (int jj = limit ; jj != nbins ; ++jj) weights[jj] = 1.0;
        } while (std::fabs(alpha1-alpha) > 0.0005);
        _alpha = alpha;
        _beta = beta;
        
        
    }
    
    //-------ML Britton-test-------------------------------------------------------

#include "intensity_target.h"
        
    /*! constructor for the ML Britton test.
     \param isig Experimental intensities
     \param twinop Twinning operator to be tested
     \param ncs tNCS vector as a Coord_frac (default NULL)
     \param reso Resolution limit for calculation (default 3.0A)
     \param nbins Number of bins for the CDF (default 20)
     \return type MLBritton_test */
    /*template <class T> MLBritton_test::MLBritton_test( const clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, const clipper::Isymop& twinop, clipper::Symop ncs, clipper::Range<clipper::ftype> reso, int nbins ) {
        _reso=reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
        _twinop = twinop;
        _ncs = ncs;
        //norm using lots of bins
        unsigned int nprm = 60;
        std::vector<double> params( nprm, 1.0 );
        clipper::BasisFn_spline basis( isig, nprm, 2.0 );
        TargetFn_meanInth<clipper::data32::I_sigI> target( isig, 1 );
        clipper::ResolutionFn norm( isig.hkl_info(), basis, target, params );
        
        calc(isig,twinop,norm,ncs);        
    } */
    
    /*! calculate Britton-test.
     \param isig Intensity data
     \param twinop Twinning operator as Isymop
     \param ncs tNCS vector as Coord_frac (default NULL)
     \param reso Resolution limit for calculation (default 3.0A)
     \return alpha */
    /*template <class T> const clipper::ftype MLBritton_test::operator() (const clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, const clipper::Isymop& twinop, clipper::Symop ncs, clipper::Range<clipper::ftype> reso ) {
        _reso=reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
        _twinop = twinop;
        _ncs = ncs;
        //norm using lots of bins
        unsigned int nprm = 60;
        std::vector<double> params( nprm, 1.0 );
        clipper::BasisFn_spline basis( isig, nprm, 2.0 );
        TargetFn_meanInth<clipper::datatypes::I_sigI<T> > target( isig, 1 );
        clipper::ResolutionFn norm( isig.hkl_info(), basis, target, params );
        
        calc(isig,twinop,norm,ncs);
        return _alpha;
    }*/
	
	/*! calculate Britton-test.
     \param isig Intensity data
     \param twinop Twinning operator as Isymop
     \param ncs tNCS vector as Coord_frac (default NULL)
     \param reso Resolution limit for calculation (default 3.0A)
     \return alpha */
    /*template <class T> const clipper::ftype MLBritton_test::operator() (const clipper::HKL_data< clipper::datatypes::I_sigI<T> >& isig, const clipper::Isymop& twinop, clipper::Range<clipper::ftype>& reso, clipper::Symop ncs ) {
        _reso=reso;
        if (_reso.max() < 0 ) _reso = clipper::Range<clipper::ftype>(isig.hkl_info().invresolsq_range() );
        _twinop = twinop;
        _ncs = ncs;
        //norm using lots of bins
        unsigned int nprm = 60;
        std::vector<double> params( nprm, 1.0 );
        clipper::BasisFn_spline basis( isig, nprm, 2.0 );
        TargetFn_meanInth<clipper::datatypes::I_sigI<T> > target( isig, 1 );
        clipper::ResolutionFn norm( isig.hkl_info(), basis, target, params );
        
        calc(isig,twinop,norm,ncs);
        return _alpha;
    }*/
	
        
    /*! output name of operator
     */
    const std::string MLBritton_test::description() const {
        clipper::String s;
        clipper::Mat33<int> mat33=_twinop.rot();
        MatrixToString(mat33,s);
        return s;
    }

    /*! print summary from Britton-tests
     */
    void MLBritton_test::summary() const {
        printf("\nML Britton test for twinning: \n");
        std::cout << "Twinning Operator " << description();
        if (_alpha > 0.05) {
            printf(", fraction = %5.2f (DeltaR = %5.2f)\n\n",_alpha,_dr );
        } else {
            std::cout << ", NO twinning" << std::endl << std::endl;
        }
    }
    
    /*! print loggraphs from Britton-tests
     */
    void MLBritton_test::loggraph() const {
        printf("\n");
        
        printf("$TABLE: ML Britton test for twinning (operator %s) alpha = %5.2f:\n", description().c_str(), _alpha );
        printf("$GRAPHS");
        printf(": ML Murray-Rust plot |Britton|:A:1,2:\n");
        printf("$$ alpha MLBritton$$\n$$\n");
        
        for (int j=0;j!=_pdf.size();++j) {
            double x = (double(j))/(2.0*double(_pdf.size() ));
            printf("%f %f\n", x, double(_pdf[j]) );
        }
        printf("$$\n\n");
        //}
        
    }
    
    template<class T, template<class> class D> clipper::ftype MLBritton_test::point(const clipper::HKL_data< D<T> >& isig, const clipper::ResolutionFn& norm, const clipper::Isymop& twinop, const clipper::ftype& alpha, const clipper::ftype& Dr) {
        typedef clipper::HKL_data_base::HKL_reference_index HRI;
        
        const clipper::HKL_info& hklinf = isig.hkl_info();
        const clipper::Spacegroup& spgr = hklinf.spacegroup();

        double width = 0.5;
        int nprm = (2.0*6.0/width)+1; //for the standard Simpson's rule the number of intervals must be even
        double a12 = (1.0 - 2.0*alpha);
        double a1 = 1.0 - alpha;
        double llk(0.0);
		const double SMALL(exp(-10.0));
        
        double scalefac=12.0;


        clipper::Vec3<int> jhkl, jhkl2;
        for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
            if ( !isig[ih].missing() && !ih.hkl_class().centric() ) {
				if (_reso.contains(ih.invresolsq() ) ) {
					clipper::HKL hkl = ih.hkl();
					jhkl[0] = hkl.h();
					jhkl[1] = hkl.k();
					jhkl[2] = hkl.l();
					clipper::HKL twin;
					jhkl2 = twinop*jhkl;
					//twin.h() = jhkl2[0];
					//twin.k() = jhkl2[1];
					//twin.l() = jhkl2[2];
					twin.h() = jhkl2[0]/scalefac;
					twin.k() = jhkl2[1]/scalefac;
					twin.l() = jhkl2[2]/scalefac;
					if (!isig[twin].missing()) {
						double DD = std::exp(-std::pow(clipper::Util::pi(),3.0)*std::pow(Dr,2.0)*hkl.invresolsq(hklinf.cell()));
						// Assume normal distribution of the coordinate errors Luzzatti (1952) Acta Cryst 5 802.
						//NT++;
						double eps1 = ih.hkl_class().epsilon();
						double eps2 = spgr.hkl_class(twin).epsilon(); 
						double w1 = 2.0/ih.hkl_class().epsilonc();
						double w2 = 2.0/spgr.hkl_class(twin).epsilonc();
						double I1 = isig[ih].I()/eps1;
						double I2 = isig[twin].I()/eps2;
						double S1 = isig[ih].sigI()/eps1;
						double S2 = isig[twin].sigI()/eps2;
						double sigma1 = norm.f(ih);
						double sigma2 = norm.f(ih);
						//double j1 = I1-S1*S1/(sigma1); // mean of J1 integral
						//double j2 = I2-S2*S2/(sigma2); // mean of J2 integral
						double j1 = I1-S1*S1/(sigma1*(1.0-DD*DD)); // mean of J1 integral
						double j2 = I2-S2*S2/(sigma2*(1.0-DD*DD)); // mean of J2 integral
						
						double n = clipper::Util::twopi()*sigma1*(1.0-DD*DD);
						//std::cout << "Norm: " << D << " " << sigma1 << " = " << n << std::endl;
						
						double num(0.0);
						double denom(0.0);
						
						double lrange2 = std::max(0.0,j2 - width*clipper::ftype(nprm/2)*S2);
						double urange2 = j2 + width*clipper::ftype(nprm/2)*S2;
						double interval2 = (urange2-lrange2)/clipper::ftype(nprm-1);
						double lrange1 = std::max(0.0,j1 - width*clipper::ftype(nprm/2)*S1);
						double urange1 = j1 + width*clipper::ftype(nprm/2)*S1;
						double interval1 = (urange1 - lrange1)/clipper::ftype(nprm-1);
						
						//denominator
						for (int ii = 0; ii != nprm ; ++ii ) {
							// integrals based on limits (1-a)I2/a and aI2/(1-a), adjusted for none zero mean and F&W gaussian approx
							double scale1 = ( ii == 0 || ii == (nprm-1) ) ? 1.0 : ((ii % 2 ) ? 4.0 : 2.0) ; //section width width*sigma, Simpsons Rule
							double x2 = lrange2+ii*interval2;
							
							double cov(0.0);
							double denom1(0.0);
							
							for (int jj = 0; jj != nprm ; ++jj ) {
								double scale = ( jj == 0 || jj == (nprm-1) ) ? 1.0 : ((jj % 2 ) ? 4.0 : 2.0) ; //section width interval, Simpsons Rule
								scale *= scale1;
								double x1 = lrange1+clipper::ftype(jj)*interval1;
								double i1 = a1*x1-alpha*x2;
								double i2 = a1*x2-alpha*x1;
								//cov = ( D > 0.0001 ) ?  D*std::sqrt(std::abs(i1*i2))/(a12*(1.0-D*D)*std::sqrt(sigma1*sigma2)) : 0.0;
								//cov = (i1*i2 > 0.0000 ) ? cov : 0.0;
								//cov = 0.0;
								denom1 += (scale/9.0)*interval1*interval2*exp(-0.5*std::pow( (x1-I1),2.0)/(S1*S1) -0.5*std::pow( (x2-I2),2.0)/(S2*S2)
																			  -x1/((1.0-DD*DD)*sigma1)
																			  -x2/((1.0-DD*DD)*sigma2));
								//                                             +2.0*cov);
							}
							denom += denom1;
							//limits
							double upper = ( alpha > 0.00001 ) ? (1.0-alpha)*x2/alpha : 99999.0*x2;
							double lower = alpha*x2/(1.0-alpha);
							if (upper < lrange1 || lower > urange1) continue;
							double lrange = std::max(lower,lrange1);
							double urange = std::min(upper,urange1);
							int nint = (urange - lrange)/interval1+1;
							if (nint <= 0 ) continue;
							nint = ( nint < 5 ) ? 5 : (nint/2)*2+1;
							double interval = (urange - lrange)/clipper::ftype(nint-1);
							
							double num1(0.0),num2(0.0);
							for (int jj = 0; jj != nint ; ++jj ) {
								double scale = ( jj == 0 || jj == (nint-1) ) ? 1.0 : ((jj % 2 ) ? 4.0 : 2.0) ; //section width interval, Simpsons Rule
								scale *= scale1;
								double x1 = lrange+clipper::ftype(jj)*interval;
								double i1 = (a1*x1-alpha*x2)/a12;
								double i2 = (a1*x2-alpha*x1)/a12;
								double i1i2 = i1*i2;
								cov = ( DD > 0.00001 && i1i2 > 0.00001) ? DD*std::sqrt(i1i2)/((1.0-DD*DD)*std::sqrt(sigma1*sigma2)) : 0.0;
								double v = -0.5*std::pow( (x1-I1),2.0)/(S1*S1) -0.5*std::pow( (x2-I2),2.0)/(S2*S2)
								-x1/((1.0-DD*DD)*sigma1)
								-x2/((1.0-DD*DD)*sigma2);
								num1 += (scale/9.0)*interval*interval2*exp(v+2.0*cov);
								num2 += (scale/9.0)*interval*interval2*exp(v);
							}
							num += num1;
							denom += num1 - num2;
							
						}
						double m = std::abs(num)/std::abs(denom);
						llk += ( m > SMALL ) ? w1*w2*(-log(m)+log(a12)) : w1*w2*(10.0+log(a12));
					}
				}
			}
		}
        return clipper::ftype(llk);
    }
    
    
    template<class T, template<class> class D>  void MLBritton_test::calc(const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twinop, const clipper::ResolutionFn& norm, const::clipper::Symop& ncs) {
        typedef clipper::HKL_data_base::HKL_reference_index HRI;
        clipper::Spacegroup spgr = isig.hkl_info().spacegroup();
        
        clipper::ftype width(0.5);
        int nprm(19); //for the standard Simpson's rule the number of intervals must be even
        //using 4.5 sigma
        clipper::ftype scalefac(12.0);
        
        int nbins=_pdf.size();
        // Britton test
        clipper::Vec3<int> jhkl, jhkl2;
        for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
            if ( !isig[ih].missing() && !ih.hkl_class().centric() ) {
				if (_reso.contains(ih.invresolsq() ) ) {
					clipper::HKL hkl = ih.hkl();
					jhkl[0] = hkl.h();
					jhkl[1] = hkl.k();
					jhkl[2] = hkl.l();
					clipper::HKL twin;
					jhkl2 = twinop*jhkl;
					//twin.h() = jhkl2[0];
					//twin.k() = jhkl2[1];
					//twin.l() = jhkl2[2];
					twin.h() = jhkl2[0]/scalefac;
					twin.k() = jhkl2[1]/scalefac;
					twin.l() = jhkl2[2]/scalefac;
					
					if (!isig[twin].missing()) {
						if (!isig[twin].missing()) {
							//NT++;
							double eps1 = ih.hkl_class().epsilon();
							double eps2 = spgr.hkl_class(twin).epsilon();
							double w1 = 2.0/ih.hkl_class().epsilonc();
							double w2 = 2.0/spgr.hkl_class(twin).epsilonc();
							double I1 = isig[ih].I()/eps1;
							double I2 = isig[twin].I()/eps2;
							double S1 = isig[ih].sigI()/eps1;
							double S2 = isig[twin].sigI()/eps2;
							clipper::ftype j1 = I1-S1*S1/norm.f(ih); // mean of J1 integral
							clipper::ftype j2 = I2-S2*S2/norm.f(ih); // mean of J2 integral
                            clipper::ftype r3 = ( 0.0 > j1 ) ? 1.0 + erf(-j1/(std::sqrt(2.0)*S1)) : erfc(j1/(std::sqrt(2.0)*S1));
                            clipper::ftype r4 = ( 0.0 > j2 ) ? 1.0 + erf(-j2/(std::sqrt(2.0)*S2)) : erfc(j2/(std::sqrt(2.0)*S2));
                            clipper::ftype d3 = 2.0-r3;
                            clipper::ftype d4 = std::sqrt(clipper::Util::twopi())*S2*0.5*(2.0-r4);
                            //std::cout << j1 << " " << j2 << " " << d3 << " " << d4 << std::endl;
							//alphas
							for (int j = 0; j != nbins ; ++j) {
								clipper::ftype alpha = 0.5*double(j)/double(nbins);
								clipper::ftype a12 = (1.0 - 2.0*alpha);
								clipper::ftype a1 = (1.0 - alpha);
								clipper::ftype tmp = 0.0; //store integral
								// sum around I2
								clipper::ftype lrange2 = std::max(0.0,j2 - width*clipper::ftype(nprm/2)*S2);
								clipper::ftype urange2 = j2 + width*clipper::ftype(nprm/2)*S2;
								clipper::ftype interval2 = (urange2-lrange2)/clipper::ftype(nprm-1);
								for (int ii = 0; ii != nprm ; ++ii ) {
									// integrals based on limits (1-a)I2/a and aI2/(1-a), adjusted for none zero mean and F&W gaussian approx
									clipper::ftype scale = ( ii == 0 || ii == (nprm-1) ) ? 1.0 : ((ii % 2 ) ? 4.0 : 2.0) ; //section width width*sigma, Simpsons Rule
									clipper::ftype x2 = lrange2+ii*interval2;                                    
									clipper::ftype upper = (1.0-alpha)*x2/alpha;
									clipper::ftype lower = alpha*x2/(1.0-alpha);
                                    //std::cout << lower << " - " << x2 << " - " << upper << " " << scale << std::endl;
									clipper::ftype r2 = ( upper >= j1 ) ? 1.0 + erf((upper-j1)/(std::sqrt(2.0)*S1)) : erfc((j1-upper)/(std::sqrt(2.0)*S1));
									clipper::ftype r1 = ( lower >= j1 ) ? 1.0 + erf((lower-j1)/(std::sqrt(2.0)*S1)) : erfc((j1-lower)/(std::sqrt(2.0)*S1));
									//clipper::ftype r3 = ( 0.0 > j1 ) ? 1.0 + erf(-j1/(std::sqrt(2.0)*S1)) : erfc(j1/(std::sqrt(2.0)*S1));
									//clipper::ftype r4 = ( 0.0 > j2 ) ? 1.0 + erf(-j2/(std::sqrt(2.0)*S2)) : erfc(j2/(std::sqrt(2.0)*S2));
                                   // std::cout << r1 << " " << r2 << std::endl;
                                    clipper::ftype d1 = r2 -r1;
                                    clipper::ftype d2 = ( d1 > 2.0 ) ? 2.0 : ( d1 < 0.0) ? 0.0 : d1;
									tmp += (scale/3.0)*((d2)/d3)*interval2*exp(-0.5*std::pow( (x2-j2),2.0)/(S2*S2) )/d4; // need quadrature for integtral (numerator and denominator contain S2)
								}
                                //std::cout << "MLBritton::calc " << j << " " << tmp << std::endl;
                                _pdf[j] += ( tmp != clipper::Util::nan() && tmp > 0.0 ) ? w1*w2*(-log(tmp)+log(a12)) : 0.0;
                                //std::cout << "MLBritton::calc " << j << " " << _pdf[j] << std::endl;
							}
						}
					}
                }
            }
        }
        //throw;
        clipper::ftype alpha(0.0),alpha1(0.0),llk(0.0),dr(100.0);
        llk = _pdf[0];
        for (int ii=0 ; ii != _pdf.size() ; ++ii ) {
            if ( _pdf[ii] < llk ) {
                llk = _pdf[ii];
                alpha = double(ii)/(2.0*_pdf.size());
            }
        }
        _alpha = alpha;
        _llk = llk;
        _dr = dr;

        // coincidence of ncs rotation operator and twinop
        clipper::ftype product(0.0);
        if (clipper::Coord_frac(ncs.trn()).lengthsq(isig.hkl_info().cell()) > 0.1 ) {
            clipper::Vec3<clipper::ftype> vect = ncs.rtop_orth(isig.hkl_info().cell() ).trn();
            clipper::Mat33<int> tmp = twinop.rot(); 
            clipper::Rotation rot(clipper::Mat33<clipper::ftype>(
                                                                 clipper::ftype(tmp(0,0) )/12.0,
                                                                 clipper::ftype(tmp(0,1) )/12.0,
                                                                 clipper::ftype(tmp(0,2) )/12.0,
                                                                 clipper::ftype(tmp(1,0) )/12.0,
                                                                 clipper::ftype(tmp(1,1) )/12.0,
                                                                 clipper::ftype(tmp(1,2) )/12.0,
                                                                 clipper::ftype(tmp(2,0) )/12.0,
                                                                 clipper::ftype(tmp(2,1) )/12.0,
                                                                 clipper::ftype(tmp(2,2) )/12.0));
            clipper::Vec3<clipper::ftype> vecr(rot.x(),rot.y(),rot.z() );
            product = std::fabs(vecr.unit()*vect.unit());
        }

        if ( product < 0.95 ) return;
        
        //grid search time
        clipper::ftype rcalc[] = {3.0,2.0,1.0,0.5,0.1};
        _llk = point(isig,norm,twinop,_alpha,_dr);

        int ni = sizeof(rcalc)/sizeof(clipper::ftype);
        clipper::ftype step = 0.5/nbins;
        
        for (int i = 0 ; i != ni ; ++i ) {
            clipper::ftype llk1(0.0), llk(0.0), tmp(0.0), a(alpha);
            int j(0);
            do {
                a = alpha +clipper::ftype(j)*step;
				if (a > 0.475) a = 0.475;
                llk1 = point(isig,norm,twinop,a,rcalc[i]);
                if ( llk > llk1 ) {
                    llk = llk1;
                    alpha = a;
                } 
                --j;
            } while ( llk == llk1 &&  a > 0.0);
            if ( llk < _llk ) {
                _llk = llk;
                _alpha = alpha;
                _dr = rcalc[i];
            } else if ( std::abs(llk-_llk) > _llk*0.00001) {
                break;
            }
        }
        return;
    }

    //-------R-test-------------------------------------------------------------
    
	
    /*! output name of operator
     */
    const std::string R_test::description() const {
        clipper::String s;
        clipper::Mat33<int> mat33=_twinop.rot();
        MatrixToString(mat33,s);
        return s;
    }
    
    /*! print summary from H-tests
     */
    void R_test::summary() const {
        printf("\nR value for twinning: \n");
        std::cout << "Twinning Operator " << description() << " with R " << std::setprecision(2) <<  _Rav;
    }
    
    /*! print loggraphs from H-tests
     */
    void R_test::loggraph() const {
		const int n(size() );
        printf("\n");
		printf("$TABLE: R test for twinning (operator %s) R = %5.2f:\n", description().c_str(), _Rav );
		printf("$GRAPHS");
		printf(":|R| vs resoltion :A:1,2:\n");
		for (int i=0; i!= n; ++i) {
			double res = ResolStats_base::operator[](i);
			printf("%10.6f %10.6f\n", res, _r[i]);
		}
		printf("$$\n\n");        
    }
    
   /* template<class T, template<class> class D>  void R_test::calc(const clipper::HKL_data< D<T> >& isig, const clipper::Isymop& twinop, const clipper::Range& reso) {
        typedef clipper::HKL_data_base::HKL_reference_index HRI;
        init(isig,false);
			 
        double HT(0.0), NT(0.0);
        double scalefac = 12.0;
        _r.resize(this->size(),0.0);
		nt.resize(this->size(),0.0);
        clipper::Vec3<int> jhkl, jhkl2;
        for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
            if ( !isig[ih].missing() && !ih.hkl_class().centric() ) {
				if (reso.contains(ih.invresolsq() ) ) {
					clipper::HKL hkl = ih.hkl();
					jhkl[0] = hkl.h();
					jhkl[1] = hkl.k();
					jhkl[2] = hkl.l();
					clipper::HKL twin;
					jhkl2 = twinop*jhkl;
					twin.h() = jhkl2[0]/scalefac;
					twin.k() = jhkl2[1]/scalefac;
					twin.l() = jhkl2[2]/scalefac;
					if (!isig[twin].missing()) {
						double I1 = isig[ih].I();
						double I2 = isig[twin].I();
						clipper::ftype eps = (this->is_intensity() ) ? 1.0/ih.hkl_class().epsilonc() : 1.0/std::sqrt(ih.hkl_class().epsilonc());
						int bin = clipper::Util::bound( 0,clipper::Util::intf( clipper::ftype(this->size() ) * ResolStats_base::operator(i) ), this->size() );
						double H = 0.0;
						if ( I1 != 0.0 && I2 != 0.0) H = (I2-I1)/(I2+I1);
						if (fabs(H) < 1){
							HT +=fabs(H)*eps;
							NT += eps;
							_r[bin] += fabs(H)*eps;
							nt[bin] += eps;
						}
					}
                }
            }
		}
		for (int i = 0; i != this->size(); ++i ) _r[i] /= nt[i];
		_Rav = HT/NT;
    }*/
    	
	//******TwinAnalysis**************************************************************************
	
	//-------constructor--------------------------------------------------------
	
	/*template<class T, template<class> class D> TwinAnalysis::TwinAnalysis(clipper::HKL_data<D<T> >& hkldata)
	{
		std::vector<clipper::Symop> ncsops;
		clipper::Range<clipper::ftype> range;
		init(hkldata,ncsops,range);
	}*/
	
	/*template<class T, template<class> class D> TwinAnalysis::TwinAnalysis(clipper::HKL_data<D<T> >& hkldata, std::vector<clipper::Symop>& ncsops)
	{
		clipper::Range<clipper::ftype> range;
		init(hkldata,ncsops,range);
	}*/
	
	
	//template<class T, template<class> class D, class TT> TwinAnalysis::TwinAnalysis(clipper::HKL_data<D<T> >& hkldata, std::vector<clipper::Symop>& ncsops, clipper::Range<TT>& range)
	/*template<class D, class TT> TwinAnalysis::TwinAnalysis(clipper::HKL_data<D >& hkldata, std::vector<clipper::Symop>& ncsops, clipper::Range<TT>& range)
	{
		init(hkldata,ncsops,range);
	}*/
		
	/*template<class T, template<class> class D, class TT> void TwinAnalysis::init(clipper::HKL_data<D<T> >& hkldata, std::vector<clipper::Symop>& ncsops, clipper::Range<TT>& range)
	{
		_base = (&hkldata);
		_range = (range.min() > range.max() ) ? _base->invresolsq_range() : range;
		_ncsops = ncsops;
		
		TwinSymops ts(_base->cell(),_base->spacegroup() );
		
		_ltest(hkldata,_ncsops,_range);
		
		_moments(hkldata,_range);

		_htests.resize(ts.size() );
		for (int i = 0; i != ts.size() ; ++i ) _htests[i](hkldata,ts[i],_range);
		
		_btests.resize(ts.size() );
		for (int i = 0; i != ts.size() ; ++i ) _btests[i](hkldata,ts[i],_range);
		
		_mlbtests.resize(ts.size() );
		for (int i = 0; i != ts.size() ; ++i ) {
			clipper::ftype product(0.0);
			int jp;
			if ( _ncsops.size() ) {
				for (int j=0; j != _ncsops.size() ; ++j) {
					clipper::Vec3<clipper::ftype> vect = _ncsops[i].rtop_orth(_base->cell() ).trn();
					clipper::Mat33<int> tmp = ts[i].rot(); 
					clipper::Rotation rot(clipper::Mat33<clipper::ftype>(
																		 clipper::ftype(tmp(0,0) )/12.0,
																		 clipper::ftype(tmp(0,1) )/12.0,
																		 clipper::ftype(tmp(0,2) )/12.0,
																		 clipper::ftype(tmp(1,0) )/12.0,
																		 clipper::ftype(tmp(1,1) )/12.0,
																		 clipper::ftype(tmp(1,2) )/12.0,
																		 clipper::ftype(tmp(2,0) )/12.0,
																		 clipper::ftype(tmp(2,1) )/12.0,
																		 clipper::ftype(tmp(2,2) )/12.0));
					clipper::Vec3<clipper::ftype> vecr(rot.x(),rot.y(),rot.z() );
					clipper::ftype p = vecr.unit()*vect.unit();
					if (p > product) {
						product = p;
						jp = j;
					}
				}
			}
			if (product > 0.95 ) _mlbtests[i](hkldata,ts[i],_ncsops[jp],_range);
			else _mlbtests[i](hkldata,ts[i],_range);
		}
	}*/
	
	
	//-------output to std::out-------------------------------------------------
	
	void TwinAnalysis::output() {
		printf("\nTWINNING ANALYSIS:\n\n");
		
        
        std::cout << "Global twinning statistics.\n\nThese tests rely on the fact that it is highly improbably that very weak or very strong reflections will coincide, therefore, the tails for the distribution of twinned datasets will be less pronounced\n";
        
		printf("\nData truncated to %6.2f - %6.2f A resolution\n",1.0/std::sqrt(_range.min()), 1.0/std::sqrt(_range.max()));
		
		_ca.output();
		
		std::cout << std::endl;
        std::cout << "The culmulative intensity, N(Z), plot is diagnostic for both twinning and tNCS.  For twinned data there are fewer weak reflections, therefore, N(Z) is sigmoidal for twinned data.  However, if both twinning and tNCS are present, the effects may cancel each out. Therefore the results of the L-test and patterson test should be consulted" << std::endl << std::endl;
		
		_ltest.summary();
		_ltest.loggraph();
		
        std::cout << std::endl << "The Cumulative |L| plot for acentric data, where L = (I1-I2)/(I1+I2). This depends on the local difference in intensities.  The difference operators used link to the neighbouring reflections taking into account possible tNCS operators." << std::endl;
        std::cout << "Note that this estimate is not as reliable as obtained via the H-test or ML Britton test if twin laws are available.  However, it is less prone to the effects of anisotropy than the H-test" << std::endl;
        std::cout << std::endl << "Reference: Padilla, Yeates. A statistic for local intensity differences: robustness to anisotropy and pseudo-centering and utility for detecting twinning. Acta Cryst. D59, 1124-30, 2003." << std::endl << std::endl;
        
		printf("\nMean acentric moments I from input data:\n\n");
		printf("  <I^2>/<I>^2 = %6.3f (Expected = %6.3f, Perfect Twin = %6.3f)\n", _moments.acentric_second(), _moments.theo_untwinned_acentric_second(), _moments.theo_perfect_acentric_second() );
		printf("  <I^3>/<I>^3 = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", _moments.acentric_third(), _moments.theo_untwinned_acentric_third(), _moments.theo_perfect_acentric_third() );
		printf("  <I^4>/<I>^4 = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", _moments.acentric_fourth(), _moments.theo_untwinned_acentric_fourth(), _moments.theo_perfect_acentric_fourth());
		std::cout << std::endl;
		
		_moments.loggraph();
		
		/*if (anisobysymm && anisodemo) {
			Moments<data32::I_sigI> mc(ianiso,range_Twin);
			printf("\n\nMean acentric moments I from anisotropically corrected data:\n\n");
			printf("  <I^2>/<I>^2 = %6.3f (Expected = %6.3f, Perfect Twin = %6.3f)\n", mc.acentric_second(), mc.theo_untwinned_acentric_second(), mc.theo_perfect_acentric_second() );
			printf("  <I^3>/<I>^3 = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", mc.acentric_third(), mc.theo_untwinned_acentric_third(), mc.theo_perfect_acentric_third() );
			printf("  <I^4>/<I>^4 = %6.3f (Expected value = %6.3f, Perfect Twin = %6.3f)\n", mc.acentric_fourth(), mc.theo_untwinned_acentric_fourth(), mc.theo_perfect_acentric_fourth());	
			m_fraction = mc.fraction();
		}
		std::cout << std::endl << std::endl; */
		
		_ts.output();
		
        std::cout << std::cout << "The appearance of twinning operators only indicates that the crystal symmetry and lattice symmetry permit twinning.  It does not mean that there is twinning present.  Only the presence of statistics consistent with twinning gives a strong indicator." << std::endl << std::endl;
		//htest plots
		if (_htests.size()) {
            std::cout << "Twinning operator based tests:" << std::endl << std::endl;
            
            std::cout << "H-test: Cumulative plot of H=|I-T(I)|/(I-T(I)) for twin related reflections.  This should be linear with slope 1/(1-2a)." << std::endl << std::endl;
			printf("\n$TABLE: H test for twinning\n");
			printf("$GRAPHS");
			for (int i = 0; i != _htests.size() ; ++i ) {
				printf(": cumulative distribution function for |H| (operator %s) alpha = %5.2f:0|1x0|1:1,2,3,4,5,6,%d:\n",_htests[i].description().c_str(),_htests[i].fraction(),i+7);
			}
			std::cout << "$$ |H| 0.4 0.3 0.2 0.1 0.0";
			for (int i = 0; i != _htests.size() ; ++i ) {
				std::string tmp = _htests[i].description();
				tmp.erase(std::remove(tmp.begin(), tmp.end(), ' '), tmp.end());
				std::cout << " " << tmp; 
			}
			std::cout << "$$\n$$\n";
			std::cout << "0.00 0.0 0.0 0.0 0.0 0.0";
			for (int i = 0; i != _htests.size() ; ++i ) {
				std::cout << " 0.00"; 
			}
			std::cout << std::endl;
			int dim = _htests[0].array().size();
			for (int j = 0 ; j != dim-1 ; ++j) {
				std::cout << double(j+1)/double(dim) << "  -   -   -   -   -";
				for (int i = 0; i != _htests.size() ; ++i ) {
					std::cout << std::fixed << std::setw(4) << std::setprecision(2) << " " << (_htests[i].array())[j];
				}
				std::cout << std::endl;
			}
			std::cout << "1.00 5.0 2.5 1.67 1.25 1.0";
			for (int i = 0; i != _htests.size() ; ++i ) {
				std::cout << " 1.00"; 
			}
			std::cout << "\n$$\n\n";
		}
			
		/*for (int i = 0; i != _htests.size() ; ++i ) {
			_htests[i].loggraph();
		}*/
		
		if (_btests.size()) {
            std::cout << "Britton plot: Plot of number of negative detwinned intensities." << std::endl << std::endl;
			printf("\n$TABLE: Britton plot for twinning\n");
			printf("$GRAPHS");
			for (int i = 0; i != _btests.size() ; ++i ) {
				printf(": aI1+(1-a)I2 > 0 (operator %s) alpha = %5.2f:A:1,%d:\n",_btests[i].description().c_str(),_btests[i].fraction(),i+2);
			}
			std::cout << "$$ alpha";
			for (int i = 0; i != _btests.size() ; ++i ) {
				std::string tmp = _btests[i].description();
				tmp.erase(std::remove(tmp.begin(), tmp.end(), ' '), tmp.end());
				std::cout << " " << tmp;
			}
			std::cout << "$$\n$$\n";
			int dim = _btests[0].array().size();
			for (int j = 0 ; j != dim ; ++j) {
				std::cout << double(j)/(2.0*double(dim));
				for (int i = 0; i != _btests.size() ; ++i ) {
					std::cout << "  " << std::right << std::fixed << std::setw(8) << std::setprecision(2) << (_btests[i].array())[j];
				}
				std::cout << std::endl;
			}
			std::cout << "$$\n\n";
		}
		
		/*for (int i = 0; i != _btests.size() ; ++i ) {
			_btests[i].loggraph();
		}*/
		
		if (_mlbtests.size()) {
            std::cout << "ML-Britton: Plot of number of negative detwinned intensities.  The ML element corrects for the sigma in the observed intensity and for the effects of a single tNCS operator, if it is present." << std::endl << std::endl;
			printf("\n$TABLE: ML-Britton test for twinning\n");
			printf("$GRAPHS");
			for (int i = 0; i != _mlbtests.size() ; ++i ) {
				printf(": aI1+(1-a)I2 > 0 (operator %s) alpha = %5.2f:A:1,%d:\n",_mlbtests[i].description().c_str(),_mlbtests[i].fraction(),i+2);
			}
			std::cout << "$$ alpha";
			for (int i = 0; i != _mlbtests.size() ; ++i ) {
				std::string tmp = _mlbtests[i].description();
				tmp.erase(std::remove(tmp.begin(), tmp.end(), ' '), tmp.end());
				std::cout << "  " << tmp;
			}
			std::cout << "$$\n$$\n";
			int dim = _mlbtests[0].array().size();
			for (int j = 0 ; j != dim ; ++j) {
				std::cout << std::right << std::fixed << std::setw(4) << std::setprecision(2) << double(j)/(2.0*double(dim));
				for (int i = 0; i != _mlbtests.size() ; ++i ) {
					std::cout << "  " << std::right << std::fixed << std::setw(10) << std::setprecision(2) << (_mlbtests[i].array())[j];
				}
				std::cout << std::endl;
			}
			std::cout << "$$\n\n";
		}
		
		/*for (int i = 0; i != _mlbtests.size() ; ++i ) {
		 _mlbtests[i].loggraph();
		 }*/
		
		
		std::cout << "Twin fraction estimates based on global statistics:" << std::endl;
		std::cout << "  Twin fraction estimate from L-test:  " << std::setw(4) << std::setprecision(2) << _ltest.fraction() << std::endl;
		std::cout << "  Twin fraction estimate from moments: " << std::setw(4) << std::setprecision(2) << _moments.fraction() << std::endl << std::endl;
        
		std::cout << "Twin fraction estimates by twinning operator" << std::endl << std::endl;
		if ( _htests.size() > 0 ) {
            
            std::cout << "The following operator based twinning estimates analyse data with each of the possible twin operators.  If twinning is present the most likely operator will have a low RTwin score (&lt;I-T(I)&gt;/&lt;I+T(I)&gt;) and estimates of the twin fraction above 0." << std::endl << std::endl;
            
			std::cout << "-------------------------------------------------------------------------------------------------" << std::endl;
			std::cout << "| " << std::setw(40) << "operator" <<         " | L-test | |Rtwin| | H-test | Britton | ML Britton   |" << std::endl;
			std::cout << "-------------------------------------------------------------------------------------------------" << std::endl;
			for (int i = 0; i != _htests.size() ; ++i ) {
				std::cout << "| " << std::setw(40) << _htests[i].description() << " |  " << ((0.1 <= _ltest.statistic() && _ltest.statistic() < 0.440) ? " Yes " : " No  " ) 
				<< " |  " << std::setw(4) << std::setprecision(2) << (_htests[i]).statistics() << "   |  " << (_htests[i]).fraction() << "  |  " << (_btests[i]).fraction() << "  |  " << (_mlbtests[i]).fraction() << " (";
				if ( (_mlbtests[i]).deltaR() == 100.0 ) std::cout << " N/A ";
				else std::cout << std::setw(5) << (_mlbtests[i]).deltaR();
				std::cout << ") |" << std::endl;
			}
			std::cout << "-------------------------------------------------------------------------------------------------" << std::endl;
		} else {
			std::cout << "No operators found" << std::endl << std::endl;
		}
		
		clipper::ftype lval = _ltest.statistic();
		clipper::ftype alpha(0.0);
		for (int i = 0; i != _htests.size() ; ++ i) {
			alpha = std::max(alpha,_htests[i].fraction() );
		}
		
		std::cout << std::endl;
			printf("TWINNING SUMMARY\n\n");
			if (_htests.size() > 0 ) printf("Twinning fraction from H-test: %6.2f\n",alpha );
			printf("Twinning fraction from L-Test: %6.2f\n\n",_ltest.fraction() );
			
		if ( alpha <= 0.0001 ) {
            if ( 0.440 <= lval && lval <= 0.500 ) {
                printf("NO Twinning detected\n\n");
            } else if ( 0.375 <= lval && lval < 0.440 ) {
                printf("L-test statistics indicate partial twinning\n");
                printf("   It is quite likely that your data were merged into a HIGHER symmetry space group than the\n   true space group.\n");
                printf("   Please revise the space group assignment if there are problems with model building or refinement.\n");
                printf("   (Very week data in higher resolution shell may be a reason of this L-value. Run twinning tests\n  with resolution cut off 3A.)\n\n");
			}
        } else if ( 0.0001 < alpha && alpha <= 0.1 ) {
            printf("No twinning or very low twinning fraction.\n\n");
            printf("   Twinning, if any, can be safely be ignored. However, twin refinement may be attempted, but not before the\n model is completely build.\n\n");
        } else if ( 0.1 < alpha && alpha < 0.4 ) {
            printf("It is highly probable that your crystal is TWINNED.\n\n");
            printf("   Please use twin refinement after your model is almost completed and R-free is below 40%%.\n\n");
        } else {
            if ( 0.440 <= lval && lval <= 0.500 ) {
                printf("Your data might have been scaled in a LOWER symmetry space group than the true space group.\n\n");
                printf("   Please run pointless, aimless (scala) and ctruncate to revise space group assignment.\n\n");
            } else if ( 0.375 <= lval && lval < 0.440 ) {
                printf("It is highly probable that your crystal is TWINNED.\n\n");
                printf("   Please use twin refinement after your model is almost completed and R-free is below 40%%\n\n");
            }
        }
		
		
	}

	//-------output to xml------------------------------------------------------
	std::stringstream& TwinAnalysis::xml_output(std::stringstream& ss) {
		clipper::ftype lval = _ltest.statistic();
		clipper::ftype alpha(0.0);
        for (int i = 0; i != _htests.size() ; ++ i) {
            alpha = std::max(alpha,_htests[i].fraction() );
        }
		ss << "<Twinning>" << std::endl;
		_ts.xml_output(ss);  //ACCESS to elements so can write fully here
        ss << "<Comment id=\"IfTwinOpsPresent\">" << std::endl;
        ss << "The appearance of twinning operators only indicates that the crystal symmetry and lattice symmetry permit twinning.  It does not mean that there is twinning present.  Only the presence of statistics consistent with twinning gives a strong indicator." << std::endl;
        ss << "</Comment>" << std::endl;
		ss << "  <L-test>" << std::endl;
		ss << "     <Twinned>" << ( (_ltest.statistic() < 0.44) ? "Yes" : "No") << "</Twinned>" << std::endl;
		ss << "     <Result>" << _ltest.statistic() << "</Result>" << std::endl;
		ss << "     <TwinFraction>" << _ltest.fraction() << "</TwinFraction>" << std::endl;
		ss << "     <Comment id=\"L-test\">Relation between L statistics and twinning fraction:\n  Twinning fraction = 0.000  L statistics = 0.500:\n  Twinning fraction = 0.100  L statistics = 0.440:\n  Twinning fraction = 0.500  L statistics = 0.375:</Comment>" << std::endl;
		ss << "  </L-test>" << std::endl;
		_moments.xml_output(ss);
		if ( _htests.size() > 0 )             {
			for (int i = 0; i != _htests.size() ; ++i ) {
				ss << "  <Rtwin operator=\"" << _htests[i].description() << "\">" ;
				ss << std::setw(4) << std::setprecision(2) << _htests[i].statistics();
				ss << " </Rtwin>" << std::endl;
			}
			for (int i = 0; i != _htests.size() ; ++i ) {
				ss << "  <H-test operator=\"" << _htests[i].description() << "\">" << std::endl;
				ss << "    <TwinFraction>" << std::setw(4) << std::setprecision(2) << _htests[i].fraction() << "</TwinFraction>" << std::endl;
				ss << " </H-test>" << std::endl;
			}
			for (int i = 0; i != _mlbtests.size() ; ++i ) {
				ss << "  <ML-Britton operator=\"" << _mlbtests[i].description() << "\">" << std::endl;
				ss << "    <TwinFraction>" << std::setw(4) << std::setprecision(2) << _mlbtests[i].fraction() << "</TwinFraction>" << std::endl;
				ss << " </ML-Britton>" << std::endl;
			}
		}
		ss << "  <Comment id=\"TwinningSummary\">" << std::endl;
		if ( alpha <= 0.0001 ) {
            if ( 0.440 <= lval && lval <= 0.500 ) {
				ss << "NO Twinning detected";
            } else if ( 0.375 <= lval && lval < 0.440 ) {
				ss << "L-test statistics indicate partial twinning\n";
                ss << "It is quite likely that your data were merged into a HIGHER symmetry space group than the\n   true space group.\n";
                ss << "Please revise the space group assignment if there are problems with model building or refinement.\n";
                ss << "(Very week data in higher resolution shell may be a reason of this L-value. Run twinning tests\n  with resolution cut off 3A.)";
			}
        } else if ( 0.0001 < alpha && alpha <= 0.1 ) {
            ss << "No twinning or very low twinning fraction.\n";
            ss << "Twinning, if any, can be safely be ignored. However, twin refinement may be attempted, but not before the\n model is completely build.";
        } else if ( 0.1 < alpha && alpha < 0.4 ) {
            ss << "It is highly probable that your crystal is TWINNED.\n";
            ss << "Please use twin refinement after your model is almost completed and R-free is below 40%%.";
        } else {
            if ( 0.440 <= lval && lval <= 0.500 ) {
                ss << "Your data might have been scaled in a LOWER symmetry space group than the true space group.\n";
                ss << "Please run pointless, aimless (scala) and ctruncate to revise space group assignment.";
            } else if ( 0.375 <= lval && lval < 0.440 ) {
                ss << "It is highly probable that your crystal is TWINNED.\n";
                ss <<"Please use twin refinement after your model is almost completed and R-free is below 40%%";
            }
        }
		ss << std::endl;
		ss << "  </Comment>" << std::endl;
        ss << "  <Comment id=\"TwinningGlobal\">" << std::endl;
        ss << " These tests rely on the fact that it is highly improbably that very weak or very strong reflections will coincide, therefore, the tails for the distribution of twinned datasets will be less pronounced" << std::endl;
        ss << " </Comment>" << std::endl;
		_ca.xml_output(ss);
        
        ss << "<Comment id=\"NZ\">" << std::endl;
        ss << "The culmulative intensity, N(Z), plot is diagnostic for both twinning and tNCS.  For twinned data there are fewer weak reflections, therefore, N(Z) is sigmoidal for twinned data.  However, if both twinning and tNCS are present, the effects may cancel each out. Therefore the results of the L-test and patterson test should be consulted" << std::endl;
        ss << "</Comment>" << std::endl;
		_ltest.xml_graph(ss);
		
        ss << "<Comment id=\"LTestDesc\">" << std::endl;
        ss << "The Cumulative |L| plot for acentric data, where L = (I1-I2)/(I1+I2). This depends on the local difference in intensities.  The difference operators used link to the neighbouring reflections taking into account possible tNCS operators.  Note that this estimate is not as reliable as obtained via the H-test or ML Britton test if twin laws are available.  However, it is less prone to the effects of anisotropy than the H-test." << std::endl;
        ss << "</Comment>" << std::endl;
        
        ss << "<Reference id=\"LTest\">" << std::endl;
        ss << "Reference: Padilla, Yeates. A statistic for local intensity differences: robustness to anisotropy and pseudo-centering and utility for detecting twinning. Acta Cryst. D59, 1124-30, 2003" << std::endl;
        ss << "</Reference>" << std::endl;
        
        if ( _htests.size() != 0 ) {
            ss << "<Comment id=\"TwinningOperatorBased\">" << std::endl;
            ss << "The following operator based twinning estimates analyse data with each of the possible twin operators.  If twinning is present the most likely operator will have a low RTwin score (&lt;I-T(I)&gt;/&lt;I+T(I)&gt;) and estimates of the twin fraction above 0." << std::endl;
            ss << "</Comment>" << std::endl;
            ss << "<CCP4Table groupID=\"graph\" id=\"H-test\" title=\"cumulative distribution function for |H|\">" << std::endl;
            for (int i = 0; i != _htests.size() ; ++i ) {
                ss << "<plot>" << std::endl;
                ss << "<title>cumulative distribution function for |H| (operator " << _htests[i].description() << ") alpha = " << std::fixed << std::setw(5) << std::setprecision(3) << _htests[i].fraction() << " </title>" << std::endl;
                ss << "<xrange min=\"0\" max=\"1\"/>" << std::endl;
                ss << "<yrange min=\"0\" max=\"1\"/>" << std::endl;
                ss << "<plotline xcol=\"1\" ycol=\""<< i+2 << "\" >" <<  std::endl;
                ss << "<symbolsize>  2</symbolsize>" << std::endl;
                ss << "<linestyle>-</linestyle>" << std::endl;
                ss << "<colour>black</colour>" << std::endl;
                ss << "</plotline>" << std::endl;
                //5.0 2.5 1.67 1.25 1.0
                ss << "<line x1=\"   0.0000\" x2=\"   1.0000\" y1=\"   0.0000\" y2=\"   5.0000\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"   1.0000\" y1=\"   0.0000\" y2=\"   2.5000\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"   1.0000\" y1=\"   0.0000\" y2=\"   1.6700\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"   1.0000\" y1=\"   0.0000\" y2=\"   2.2500\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "<line x1=\"   0.0000\" x2=\"   1.0000\" y1=\"   0.0000\" y2=\"   1.0000\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "</plot>" << std::endl;
            }
            ss <<"<headers separator=\" \">\n |H| ";
            for (int i = 0; i != _htests.size() ; ++i ) {
                std::string tmp = _htests[i].description();
                tmp.erase(std::remove(tmp.begin(), tmp.end(), ' '), tmp.end());
                ss << " " << tmp;
            }
            ss << "\n </headers>" << std::endl;
            ss << "<data>" << std::endl;
            int dim = _htests[0].array().size();
            for (int j = 0 ; j != dim-1 ; ++j) {
                ss << std::fixed << std::setw(4) << std::setprecision(2) << double(j+1)/double(dim);
                for (int i = 0; i != _htests.size() ; ++i ) {
                    ss << "    " << std::fixed << std::setw(4) << std::setprecision(2)  << (_htests[i].array())[j];
                }
                ss << std::endl;
            }
            ss << "1.00";
            for (int i = 0; i != _htests.size() ; ++i ) {
                ss << "    1.00"; 
            }
            ss << std::endl;
            ss << "</data>" << std::endl;
            ss << "</CCP4Table>" << std::endl;
            ss << "<Comment id=\"HTest\">" << std::endl;
            ss << "Cumulative plot of H=|I-T(I)|/(I-T(I)) for twin related reflections.  This should be linear with slope 1/(1-2a)." << std::endl;
            ss << "</Comment>" << std::endl;
        }
        
        if (_mlbtests.size() != 0 ) {
            int dim = _mlbtests[0].array().size();
            ss << "<CCP4Table groupID=\"graph\" id=\"ML-Britton\" title=\"count of none negative reflections\">" << std::endl;
            for (int i = 0; i != _mlbtests.size() ; ++i ) {
                ss << "<plot>" << std::endl;
                ss << "<title>aI1+(1-a)I2 > 0  (operator " << _mlbtests[i].description() << ") alpha = " << std::fixed << std::setw(5) << std::setprecision(3) << _mlbtests[i].fraction() << " </title>" << std::endl;
                ss << "<xrange min=\"0\" max=\"0.5\"/>" << std::endl;
                ss << "<yrange min=\"None\" max=\"None\"/>" << std::endl;
                ss << "<plotline xcol=\"1\" ycol=\""<< i+2 << "\" >" <<  std::endl;
                ss << "<symbolsize>  2</symbolsize>" << std::endl;
                ss << "<linestyle>-</linestyle>" << std::endl;
                ss << "<colour>black</colour>" << std::endl;
                ss << "</plotline>" << std::endl;
                ss << "<line x1=\"   " << std::fixed << std::setw(6) << std::setprecision(4) << _mlbtests[i].fraction() <<" \" x2=\"   0.5000\" y1=\"   0.0000\" y2=\"   " <<  (_mlbtests[i].array())[dim-1] << "\" linestyle=\"-\" linecolour=\"black\"/>" << std::endl;
                ss << "</plot>" << std::endl;
            }
            ss <<"<headers separator=\" \">\n alpha ";
            for (int i = 0; i != _mlbtests.size() ; ++i ) {
                std::string tmp = _mlbtests[i].description();
                tmp.erase(std::remove(tmp.begin(), tmp.end(), ' '), tmp.end());
                ss << " " << tmp;
            }
            ss << "\n </headers>" << std::endl;
            ss << "<data>" << std::endl;
            for (int j = 0 ; j != dim-1 ; ++j) {
                ss << std::fixed << std::setw(4) << std::setprecision(2) << double(j)/(2.0*double(dim));;
                for (int i = 0; i != _mlbtests.size() ; ++i ) {
                    ss << " " << std::fixed << std::setw(10) << std::setprecision(2)  << (_mlbtests[i].array())[j];
                }
                ss << std::endl;
            }
            ss << "</data>" << std::endl;
            ss << "</CCP4Table>" << std::endl;
            ss << "<Comment id=\"MLBritton\">" << std::endl;
            ss << "Plot of number of negative detwinned intensities.  The ML element corrects for the sigma in the observed intensity and for the effects of a single tNCS operator, if it is present." << std::endl;
            ss << "</Comment>" << std::endl;
        }
		
		ss << "</Twinning>\n";
		return ss;
	}
	
	//-------CumulativeAnalysis
	double CumulativeAnalysis::_acentricideal[] = {0.0,
		0.0392106, 0.0768837, 0.1130796, 0.1478562, 0.1812692, 0.2133721, 0.2442163, 0.2738510, 0.3023237, 0.3296800,
		0.3559636, 0.3812166, 0.4054795, 0.4287909, 0.4511884, 0.4727076, 0.4933830, 0.5132477, 0.5323336, 0.5506710,
		0.5682895, 0.5852171, 0.6014810, 0.6171071, 0.6321206, 0.6465453, 0.6604045, 0.6737202, 0.6865138, 0.6988058,
		0.7106158, 0.7219627, 0.7328647, 0.7433392, 0.7534030, 0.7630722, 0.7723623, 0.7812881, 0.7898639, 0.7981035,
		0.8060200, 0.8136260, 0.8209339, 0.8279551, 0.8347011, 0.8411826, 0.8474099, 0.8533930, 0.8591416, 0.8646647};
	
	double CumulativeAnalysis::_centricideal[] = {0.0,
		0.1585194, 0.2227026, 0.2709655, 0.3108435, 0.3452792, 0.3757939, 0.4032988, 0.4283924, 0.4514938, 0.4729107,
		0.4928775, 0.5115777, 0.5291583, 0.5457398, 0.5614220, 0.5762892, 0.5904133, 0.6038561, 0.6166715, 0.6289066,
		0.6406032, 0.6517983, 0.6625250, 0.6728131, 0.6826895, 0.6921785, 0.7013024, 0.7100815, 0.7185345, 0.7266783,
		0.7345289, 0.7421010, 0.7494079, 0.7564625, 0.7632764, 0.7698607, 0.7762255, 0.7823805, 0.7883348, 0.7940968,
		0.7996745, 0.8050755, 0.8103070, 0.8153755, 0.8202875, 0.8250491, 0.8296659, 0.8341433, 0.8384867, 0.8427008};
	
	double CumulativeAnalysis::_acentrictwin[] = {0.0,
		0.0030343, 0.0115132, 0.0245815, 0.0414833, 0.0615519, 0.0842006, 0.1089139, 0.1352404, 0.1627861, 0.1912079, 
		0.2202081, 0.2495299, 0.2789524, 0.3082868, 0.3373727, 0.3660750, 0.3942806, 0.4218963, 0.4488460, 0.4750691, 
		0.5005177, 0.5251562, 0.5489585, 0.5719077, 0.5939942, 0.6152149, 0.6355726, 0.6550744, 0.6737317, 0.6915590, 
		0.7085736, 0.7247951, 0.7402450, 0.7549459, 0.7689218, 0.7821971, 0.7947971, 0.8067470, 0.8180725, 0.8287987, 
		0.8389511, 0.8485543, 0.8576328, 0.8662106, 0.8743109, 0.8819565, 0.8891694, 0.8959710, 0.9023818, 0.9084218}; 


	void CumulativeAnalysis::output() const {
		printf("$TABLE: Cumulative intensity distribution:\n");
		printf("$GRAPHS");
		printf(": Cumulative intensity distribution (Acentric and centric):N:1,2,3,4,5,6:\n$$");
		printf(" Z Acent_theor Acent_twin Acent_obser Cent_theor Cent_obser $$\n$$\n");
		double x = 0.0;
		double deltax=0.04;
		bool centric(_intensity_ord_c.ordinal(2.0) >0.1 );
		for (int i=0; i!=51; ++i) {
			if (centric) printf("%10.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n", x, _acentricideal[i], _acentrictwin[i], _intensity_ord_a.ordinal(x), _centricideal[i], _intensity_ord_c.ordinal(x));
			else printf("%10.5f %8.5f %8.5f %8.5f %8.5f -\n", x, _acentricideal[i], _acentrictwin[i], _intensity_ord_a.ordinal(x), _centricideal[i]);
			x += deltax;
		}
		printf("$$\n\n");
		
	}
	
	std::stringstream& CumulativeAnalysis::xml_output(std::stringstream& ss) const {
		ss << "<CCP4Table groupID=\"graph\" id=\"Cumulative intensity distribution\" title=\"Cumulative intensity distribution\">" << std::endl;
		ss << "<plot>" << std::endl;
		ss << "<title>Cumulative intensity distribution (Acentric and centric)</title>" << std::endl;
		ss << "<yrange min=\"0\" max=\"None\"/>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"6\" >" << std::endl;
		ss << "<symbolsize>  2</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>red</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"5\" >" << std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>red</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"4\" >" <<  std::endl;
		ss << "<symbolsize>  2</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>black</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"1\" ycol=\"3\" >" << std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>blue</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "<plotline xcol=\"  1\" ycol=\"2\" >" << std::endl;
		ss << "<symbolsize>  0</symbolsize>" << std::endl;
		ss << "<linestyle>-</linestyle>" << std::endl;
		ss << "<colour>black</colour>" << std::endl;
		ss << "</plotline>" << std::endl;
		ss << "</plot>" << std::endl;
		ss << "<headers separator=\" \">\n Z Acent_theor Acent_twin Acent_obser Cent_theor Cent_obser \n </headers>" << std::endl;
		ss << "<data>" << std::endl;
		double x = 0.0;
		double deltax=0.04;
		bool centric(_intensity_ord_c.ordinal(2.0) >0.1 );
		for (int i=0; i!=51; ++i) {
			ss << std::fixed << std::setw(6) << std::setprecision(4) << " " << x << " ";
			ss << std::fixed << std::setw(6) << std::setprecision(4) << _acentricideal[i] << " " << std::setw(6) << _acentrictwin[i] << " " << std::setw(6) << _intensity_ord_a.ordinal(x) << " " << std::setw(6) <<  _centricideal[i];
			if (centric) ss << std::fixed << std::setw(6) << std::setprecision(4) << _intensity_ord_c.ordinal(x);
			else ss << " - ";
			ss << std::endl;
			x += deltax;
		}
		ss << "</data>" << std::endl;
		ss << "</CCP4Table>" << std::endl;
		
		return ss;
	}
	
    //-------instantiate templates----------------------------------------------
    
    //template clipper::ftype L_test::operator()<clipper::ftype32>(clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >&, clipper::Range<clipper::ftype>);
    //template clipper::ftype L_test::operator()<clipper::ftype64>(clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >&, clipper::Range<clipper::ftype>);
	//template clipper::ftype L_test::operator()<clipper::ftype32>(clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >&, std::vector<clipper::Symop>&, clipper::Range<clipper::ftype>);
    //template clipper::ftype L_test::operator()<clipper::ftype64>(clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >&, std::vector<clipper::Symop>&, clipper::Range<clipper::ftype>);
    //template class L_test<clipper::ftype32>;
    //template class L_test<clipper::ftype64>;
    template clipper::ftype L_test::calc(clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >&, std::vector<clipper::Symop>&);
    
    //template clipper::ftype H_test::operator()<clipper::ftype32>(const clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >&, const clipper::Isymop&, clipper::Range<clipper::ftype>);
    //template clipper::ftype H_test::operator()<clipper::ftype64>(const clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >&, const clipper::Isymop&, clipper::Range<clipper::ftype>);
    //template class H_test<clipper::ftype32>;
    //template class H_test<clipper::ftype64>;
    template void H_test::calc(const clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >&, const clipper::Isymop&);
    
    //template class Britton_test<clipper::ftype32>;
    //template class Britton_test<clipper::ftype64>;
    //template clipper::ftype Britton_test::operator()<clipper::ftype32>(const clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >&, const clipper::Isymop&, clipper::Range<clipper::ftype>);
    //template clipper::ftype Britton_test::operator()<clipper::ftype64>(const clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >&, const clipper::Isymop&, clipper::Range<clipper::ftype>);
    template void Britton_test::calc(const clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >&, const clipper::Isymop&);

    //template class MLBritton_test<clipper::ftype32>;
    //template class MLBritton_test<clipper::ftype64>;
    //template clipper::ftype MLBritton_test::operator()<clipper::ftype32>(const clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >&, const clipper::Isymop&, clipper::Symop, clipper::Range<clipper::ftype>);
    //template clipper::ftype MLBritton_test::operator()<clipper::ftype64>(const clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >&, const clipper::Isymop&, clipper::Symop, clipper::Range<clipper::ftype>);
	//template clipper::ftype MLBritton_test::operator()<clipper::ftype32>(const clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >&, const clipper::Isymop&, clipper::Range<clipper::ftype>&, clipper::Symop);
	//template clipper::ftype MLBritton_test::operator()<clipper::ftype64>(const clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype64> >&, const clipper::Isymop&, clipper::Range<clipper::ftype>&, clipper::Symop);

    template void MLBritton_test::calc(const clipper::HKL_data<clipper::datatypes::I_sigI<clipper::ftype32> >&, const clipper::Isymop&, const clipper::ResolutionFn&, const::clipper::Symop&);
    
	//template void TwinAnalysis::init<float, clipper::datatypes::I_sigI, double>(clipper::HKL_data<clipper::datatypes::I_sigI<float> >&, std::vector<clipper::Symop, std::allocator<clipper::Symop> >&, clipper::Range<double>&);
}

