#include "cpsf_utils.h"

PeakSearch::Neighbours::Neighbours( const clipper::Xmap_base &map, const float min_distsq, const float max_distsq )
{
	clipper::Cell_descr rcd( map.cell().descr() );
	clipper::Cell_descr vcd( 1.0,1.0,1.0, rcd.alpha(), rcd.beta(), rcd.gamma() );
	clipper::Cell vcell( vcd );
	
	clipper::Coord_grid g0(-1,-1,-1);
	clipper::Coord_grid g1( 1, 1, 1);
	clipper::Grid_sampling vgrid( 1, 1, 1 );
	
	clipper::Coord_grid iu, iv, iw;
	float thisd2;
	for ( iu = g0; iu.u() <= g1.u(); iu.u()++ ) {
		for ( iv = iu; iv.v() <= g1.v(); iv.v()++ ) {
			for ( iw = iv; iw.w() <= g1.w(); iw.w()++ ) {
				thisd2 = iw.coord_frac( vgrid ).lengthsq( vcell );
				if (thisd2 > min_distsq && thisd2 < max_distsq) nlist.push_back( iw );
			}
		}
	}
}

// operate peak search, returning a vector of indices.
std::vector<int> PeakSearch::operator() ( const clipper::Xmap<float>& map ) {
	clipper::Coord_grid c, x;
	std::vector<int> peaklist;
	xmap = const_cast<clipper::Xmap<float> *>( &map );
	
	// loop over skeleton neighbours structure
	Neighbours neigh( map );
	clipper::Coord_grid g0, g1;
	clipper::Xmap<float>::Map_reference_coord i0;
	for ( clipper::Xmap<float>::Map_reference_index index = map.first() ; !index.last() ; index.next() ) {
		float value = map[index];
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
	
	// std::cout << "finished search with " << peaklist.size() << " peaks" << std::endl;
	clipper::Map_index_sort::sort_decreasing( map, peaklist );
	return peaklist;
}


std::vector<bool> IsPolar::operator() ( clipper::Spacegroup& spgr) {
	std::vector<bool> ispolar( 3, 1 );
	clipper::Coord_frac c0( 0.13, 0.17, 0.19 ), c1;
	for ( int sym = 1; sym != spgr.num_primitive_symops() ; ++sym ) {
		c1 = spgr.primitive_symop( sym )*c0;
		if ( fabs( c0.u() - c1.u() ) > 0.001f ) ispolar[0] = false;
		if ( fabs( c0.v() - c1.v() ) > 0.001f ) ispolar[1] = false;
		if ( fabs( c0.w() - c1.w() ) > 0.001f ) ispolar[2] = false;
	}
	return ispolar;
}


/* determine alternative origins based upon limited possiblities.
*/
std::vector<clipper::Coord_frac> alt_origin( clipper::Spacegroup& spgr ) {
	static float alts[] = {0.0f, 0.25f, 0.5f, 0.75f, 0.3333333f, 0.6666666f};
	
	IsPolar ispolar;
	std::vector<bool> polar = ispolar( spgr );
	
	std::vector<clipper::Coord_frac> origins;
	
	clipper::Coord_frac c0, c1;
	
	for (int u = 0 ; u != 6 ; ++u ) {
		for (int v = 0 ; v != 6 ; ++v ) {
			for (int w = 0 ; w != 6 ; ++w ) {
				c0 = clipper::Coord_frac( alts[u], alts[v], alts[w] );
				if ( ( c0.u() != 0.0f && polar[0] ) || ( c0.v() != 0.0f && polar[1] )  || ( c0.w() != 0.0f && polar[2] ) ) continue;
				int sym = 1;
				for ( ; sym != spgr.num_symops() ; ++sym ) {
					clipper::RTop<double> tmp = spgr.symop(sym) - clipper::RTop<double>::identity();
					clipper::RTop_frac t( tmp.rot() );
					c1 = t * c0;
					if ( std::fmod( c1.u(), 1.0 ) || std::fmod( c1.v(), 1.0 ) || std::fmod( c1.w(), 1.0 ) ) break;
				}
				if ( sym == spgr.num_symops() ) origins.push_back(c0);
			}
		}
	}
	return origins;
}


void harker::init( clipper::Spacegroup& spgr ) {
	for ( int isym = 0; isym != spgr.num_symops() ; ++isym ) {
		clipper::RTop_frac sym1( spgr.symop( isym ) );
		for ( int ksym = isym+1; ksym != spgr.num_symops() ; ++ksym ) {
			clipper::RTop_frac sym2( spgr.symop( ksym ) );
			if ( sym2 .equals( clipper::RTop_frac::identity(), 0.01 ) ) continue;
			clipper::RTop_frac testsym(sym1 - sym2);
			clipper::Mat33<float> mat33(testsym.rot());
			std::cout << isym << "\n" << sym1.format() << "\n" << sym2.format() << std::endl;
			for (int i = 0; i != 3 ; ++i) {
				if ( mat33(i,0)*mat33(i,0) +
					 mat33(i,1)*mat33(i,1) +
					 mat33(i,2)*mat33(i,2) < 0.1f ) {
					testsym.trn()[0] = std::fmod(testsym.trn()[0]+10.0,1.0);
					testsym.trn()[1] = std::fmod(testsym.trn()[1]+10.0,1.0);
					testsym.trn()[2] = std::fmod(testsym.trn()[2]+10.0,1.0);
					operators.push_back(testsym);
					std::cout << testsym.format() << std::endl;
					break;
				}
			}
		}
	}
}

bool harker::is_harker( clipper::Coord_frac &coord, float tol) {
	bool val(false);
	std::cout << coord.format() << " " ;
	for ( int k = 0 ; k != operators.size() ; ++k ) {
		clipper::Mat33<float> mat33(operators[k].rot() );
		clipper::Vec3<float> trn(operators[k].trn() );
		for (int i = 0; i != 3 ; ++i) {
			if ( mat33(i,0)*mat33(i,0) +
				 mat33(i,1)*mat33(i,1) +
				 mat33(i,2)*mat33(i,2) < 0.1f ) {
				if ( std::fabs(coord[i] - trn[i] ) < tol ) {
					val = true; 
				} else {
					val = false;
				}
			}
		}
		if ( val ) std::cout << val << std::endl;
		if ( val ) return val;
	}
	std::cout << val << std::endl;
	return val;
}

