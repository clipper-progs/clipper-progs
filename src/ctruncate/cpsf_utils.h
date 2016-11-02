#ifndef __CPSF_UTILS_H
#define __CPSF_UTILS_H

#include "clipper/clipper.h"
#include "clipper/clipper-contrib.h"
#include "clipper/clipper-ccp4.h"
#include "clipper/clipper-minimol.h"
#include <map>
#include <algorithm>
#if defined (__sgi) && defined (_MIPS_ISA)
extern "C" {
#include <math.h>
}
#else
#include <cmath>
#endif

/* extension to clipper namespace for alt_origins 
*/
namespace clipper {
	
	template<class T> inline Mat33<T> operator -(const Mat33<T>& m1, const Mat33<T>& m2)
{ return Mat33<T>( m1(0,0) - m2(0,0), 
				   m1(0,1) - m2(0,1), 
				   m1(0,2) - m2(0,2), 
				   m1(1,0) - m2(1,0), 
				   m1(1,1) - m2(1,1), 
				   m1(1,2) - m2(1,2), 
				   m1(2,0) - m2(2,0), 
				   m1(2,1) - m2(2,1), 
				   m1(2,2) - m2(2,2) ); 
}

template<typename T> inline RTop<T> operator -( const RTop<T>& r1, const RTop<T>& r2 ) { 
	return RTop<T>( r1.rot()-r2.rot(), r1.trn()-r2.trn() );  
}
template<typename T> inline RTop<T> operator +( const RTop<T>& r1, const RTop<T>& r2 ) { 
	return RTop<T>( r1.rot()+r2.rot(), r1.trn()+r2.trn() ); 
}
}


/* Peak searching routine using KC neighbours
*/
class PeakSearch {
	class Neighbours {
public:
		Neighbours( const clipper::Xmap_base &map, const float min_dist = 0.5, const float max_dist = 2.5 );
		clipper::Coord_grid operator[] (int i) const { return nlist[i]; }
		int size() const { return nlist.size(); }
private:
		std::vector<clipper::Coord_grid> nlist;
	};
	
	class PeakInterp {
		struct IllegalStep { };
		
	private:
		template<class T> inline clipper::Coord_map interp( const clipper::Xmap<T> &, const clipper::Coord_map& );
	public:
		template<class T> clipper::Coord_map operator() ( const clipper::Xmap<T> &map, const clipper::Coord_map &pos ) {
			clipper::Coord_map res;
			try {
				res = interp( map , pos );
			}
			catch (IllegalStep) {
				return pos;
			}
			return res;
		}
		clipper::Coord_map operator() ( const clipper::Xmap<float> &map, const int grid ) {
			clipper::Coord_map res;
			try {
				res = interp( map , clipper::Coord_map( map.coord_of( grid ) ) );
			}
			catch (IllegalStep) {
				return clipper::Coord_map( map.coord_of( grid ) );
			}
			return res;
		}
	};
	
public:
	PeakSearch(clipper::ftype zfrac=0.1) : _zfrac(zfrac), _rho0(0.0), _sigma(0.0), _sig(0.0) {}
	~PeakSearch() {}
	template<class T> const std::vector<int>& operator() ( const clipper::Xmap<T>& map, float sig=3.0 );
	clipper::ftype zero() { return _rho0; }
	int operator[](int i) { return _peaks[i]; }
	int npeaks() { return _peaks.size(); }
private:
	clipper::Xmap_base* _xmap;
	clipper::ftype _rho0;
	clipper::ftype _sigma;
	clipper::ftype _zfrac;
	clipper::ftype _sig;
	std::vector<int> _peaks;
};


template<class T> clipper::Coord_map PeakSearch::PeakInterp::interp ( const clipper::Xmap<T> &map, const clipper::Coord_map &pos ) {
	float val;
	clipper::Grad_map<T> grad;
	clipper::Curv_map<T> curv;
	
	clipper::Interp_cubic::interp_curv( map, pos, val, grad, curv );
	
	if ( std::fabs( curv.det() ) < 1.0e-6 ) throw IllegalStep(); // defend against ill conditioned matrix
	
	clipper::Curv_map<float> inv_curv(curv.inverse() );
	
	clipper::Cell cell = map.cell();
	clipper::Grid_sampling grid = map.grid_sampling();
	
	clipper::Vec3<float> tmp = inv_curv * grad;
	clipper::Coord_map step( (grid.nu() == 1) ? 0 : tmp[0], (grid.nv() == 1) ? 0 : tmp[1], (grid.nw() == 1) ? 0 : tmp[2] );
	
	if ( ( tmp[0] > cell.a()/grid.nu() ) ||
		 ( tmp[1] > cell.b()/grid.nv() ) ||
		 ( tmp[2] > cell.c()/grid.nw() ) ) throw IllegalStep();  // reject bad steps
	
	return pos - step;
}

std::vector<clipper::Coord_frac> alt_origin( clipper::Spacegroup& spgr );

/* determine if the spacegroup is polar
*/
class IsPolar {
public:
	std::vector<bool> operator() ( clipper::Spacegroup& );
};

/* determine harker sections
*/
class harker {
public:
	harker() {}
	harker( clipper::Spacegroup &spgr) { init( spgr ); };
	~harker() {}
	
	void setSpacegroup( clipper::Spacegroup &spgr) { operators.clear(); init( spgr ); }
	int num_harker() { return operators.size(); }
	clipper::RTop_frac harkerop( int n ) { return operators[n]; } 
	bool is_harker( clipper::Coord_frac &coord, float tol=0.01);
	
private:
		void init( clipper::Spacegroup& spgr );
	
	std::vector< clipper::RTop_frac > operators;
	
};

class CoordsEquiv {
public:
	CoordsEquiv( clipper::Spacegroup& spg, clipper::Cell& cell, float accuracy = 1.0f) : spacegroup_(spg), cell_(cell), accuracy_(accuracy) {};
	bool operator() ( const clipper::Atom& a0, const clipper::Atom& a1) const {
		clipper::Coord_orth c0(a0.coord_orth() ), c1(a1.coord_orth() );
		for (int sym = 0 ; sym != spacegroup_.num_symops() ; ++sym ) {
			clipper::RTop_orth trans( spacegroup_.symop( sym ).rtop_orth(cell_) );
			clipper::Coord_orth c2 = c0 - ( trans * c1 );
			clipper::Coord_frac f1 = c2.coord_frac( cell_ );
			double u = std::fmod( f1.u() + 2.0, 1.0 );
			double v = std::fmod( f1.v() + 2.0, 1.0 );
			double w = std::fmod( f1.w() + 2.0, 1.0 );
			clipper::Coord_frac f2( std::min( std::fabs( u ), std::fabs( u-1.0 ) ), 
									std::min( std::fabs( v ), std::fabs( v-1.0 ) ),
									std::min(  std::fabs( w ), std::fabs( w-1.0 ) ) );
			if ( f2.lengthsq( cell_ ) < accuracy_ ) return true;
		}
		return false;
	}	
	bool operator() ( const clipper::Coord_orth& c0, const clipper::Coord_orth& c1) const {
		for (int sym = 0 ; sym != spacegroup_.num_symops() ; ++sym ) {
			clipper::RTop_orth trans( spacegroup_.symop( sym ).rtop_orth(cell_) );
			clipper::Coord_orth c2 = c0 - ( trans * c1 );
			clipper::Coord_frac f1 = c2.coord_frac( cell_ );
			double u = std::fmod( f1.u() + 2.0, 1.0 );
			double v = std::fmod( f1.v() + 2.0, 1.0 );
			double w = std::fmod( f1.w() + 2.0, 1.0 );
			clipper::Coord_frac f2( std::min( std::fabs( u ), std::fabs( u-1.0 ) ), 
									std::min( std::fabs( v ), std::fabs( v-1.0 ) ),
									std::min(  std::fabs( w ), std::fabs( w-1.0 ) ) );
			if ( f2.lengthsq( cell_ ) < accuracy_ ) return true;
		}
		return false;
	}
	
private:
		clipper::Spacegroup& spacegroup_;
	clipper::Cell& cell_;
	float accuracy_;
};

template<class T> const std::vector<int>& PeakSearch::operator() ( const clipper::Xmap<T>& map, float sig ) {
	clipper::Coord_grid c, x;
	
	_xmap = const_cast<clipper::Xmap<T> *>( &map );
	_sig = sig;
	
	{
		float rho(0.0);
		
		std::vector<T> vals;
		vals.reserve(map.grid_asu().size());
		for ( clipper::Xmap_base::Map_reference_index index = map.first() ; !index.last() ; index.next() ) {
			vals.push_back(map[index]); }
		std::sort( vals.begin(), vals.end() );
		_rho0  = ( vals[ int(_zfrac*float(vals.size())) ] );
		
		
		int n(0);
		for ( int index = 0 ; index != vals.size() ; ++index ) {
			_sigma += std::pow(vals[index],2);
			rho += vals[index];
			++n;
		}
		_sigma -= rho*rho/n;
		_sigma /= n;
		_sigma = std::sqrt(_sigma);
	}
	
	std::vector<int> peaklist;
	// loop over skeleton neighbours structure
	float sig0 = sig*_sigma+_rho0;
	Neighbours neigh( map );
	clipper::Coord_grid g0, g1;
	clipper::Xmap_base::Map_reference_coord i0;
	for ( clipper::Xmap_base::Map_reference_index index = map.first() ; !index.last() ; index.next() ) {
		float value = map[index];
		if (value > sig0) {
			bool peak = true;
			g0 = index.coord();
			for (int i = 0; i < neigh.size(); i++ ) {
				g1 = g0 + neigh[i];
				i0 = clipper::Xmap_base::Map_reference_coord( map, g1 );
				if ( map[i0] > value ) {
					peak = false;
					break;
				}
			}
			if ( peak ) peaklist.push_back( index.index() );
		}
	}
	
	clipper::Map_index_sort::sort_decreasing( map, peaklist );
	_peaks=peaklist;
	
	/*PeakInterp pki;
	 for (std::vector<int>::const_iterator i=peaklist.begin() ; i != peaklist.end() ; ++i) {
	 _peaks.push_back(map.coord_of( peaklist[*i] ) ); */
	/*try {
	 _peaks.push_back(pki(map, peaklist[*i]).coord_frac(map.grid_sampling() ) );
	 } catch (...) {
	 _peaks.push_back(map.coord_of( peaklist[*i] ).coord_frac(map.grid_sampling() ) );
	 }*/
	//}
	return _peaks;
}

#endif
