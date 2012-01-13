
//
//    TWINLAWS
//    Copyright (C) 2004-2008 Andrey Lebedev
//
//    This code is distributed under the terms and conditions of the
//    CCP4 Program Suite Licence Agreement as a CCP4 Application.
//    A copy of the CCP4 licence can be obtained by writing to the
//    CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//
//
//    Code for calculating twinning operators from first principles
//    extracted from Sfcheck
//    and converted from FORTRAN to C
//

#include "twinlaws.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <cstdio>


// -----------------------------------------------------------------------------

void yyy_cell2met( double cell[6], double met[3][3] )
{

	double d2r = std::atan( 1.0 )/ 45.0 ;
	int m3[6] = { 0, 1, 2, 0, 1, 2 } ;

	bool rad = true ;
	bool deg = true ;
	for( int i = 3 ; i != 6 ; ++i ){
		rad &= ( cell[i] < 3.2 ) & ( cell[i] > 0 ) ;
		deg &= ( cell[i] > 3.2 ) & ( cell[i] < 180 ) ;
	}
	double coef = 1.0 ;
	if( deg ){
		coef = d2r ;
	}
	else if( ! rad ){
										throw "yyy_cell2met_:a" ;
	}

	for( int i = 0 ; i != 3 ; ++i ){
		int i1 = m3[i+1] ;
		int i2 = m3[i+2] ;
		met[i][i] = 1.0 ;
		met[i1][i2] = std::cos( cell[i+3]* coef ) ;
		met[i2][i1] = met[i1][i2] ;
	}

	for( int i = 0 ; i != 3 ; ++i ){
		for( int j = 0 ; j != 3 ; ++j ){
			met[j][i] = cell[j]* met[j][i]* cell[i] ;
		}
	}
}
// -----------------------------------------------------------------------------

void yyy_met2cell( double cell[6], double met[3][3] )
{

	double r2d = 45.0/ std::atan( 1.0 ) ;
	int m3[6] = { 0, 1, 2, 0, 1, 2 } ;

	for( int i = 0 ; i != 3 ; ++i ){
		if( met[i][i] <= 0 )						throw "yyy_met2cell_:a" ;
		cell[i] = sqrt( met[i][i] ) ;
	}

	for( int i = 0 ; i != 3 ; ++i ){
		int i1 = m3[i+1] ;
		int i2 = m3[i+2] ;
		cell[i+3] = met[i1][i2]/ ( cell[i1]* cell[i2] ) ;

		int tol = 10000 ;
		double test = cell[i+3] - met[i2][i1]/ ( cell[i1]* cell[i2] ) + tol ;
		if( test != tol* 1.0 )						throw "yyy_met2cell_:b" ;
	}

	for( int i = 3 ; i != 6 ; ++i ){
		cell[i] = r2d* std::acos( cell[i] ) ;
	}
}
// -----------------------------------------------------------------------------

void yyy_divide_int( int &n, int v[], int &d, bool &ok )
{

	ok = false ;
	if( d == 0 )								return ;
	for( int iv = 0 ; iv != n ; ++iv ){
		int j = v[iv]/ d ;
		if( v[iv] != d* j )						return ;
		v[iv] = j ;
	}
	ok = true ;
}
// -----------------------------------------------------------------------------

void yyy_invert_int( int b[3][3], int r[3][3], int &d )
{

	int c3[5] = { 0, 1, 2, 0, 1 } ;

	for( int i1 = 0 ; i1 != 3 ; ++i1 ){
		int i2 = c3[i1+1] ;
		int i3 = c3[i1+2] ;
		for( int j1 = 0 ; j1 != 3 ; ++j1 ){
			int j2 = c3[j1+1] ;
			int j3 = c3[j1+2] ;
			r[i1][j1] = b[j2][i2]* b[j3][i3] - b[j3][i2]* b[j2][i3] ;
		}
	}
	d = 0 ;
	for( int i1 = 0 ; i1 != 3 ; ++i1 ){
		d = d + b[0][i1]* r[i1][0] ;
	}
}
// -----------------------------------------------------------------------------

void yyy_find_base( int &ng, int uu_g[][3][3], int u_g[][3], int u_b[3][3] )
{

	int u[3], v[3], r[3][3] ;

	for( int i = 0 ; i != 3 ; ++i ){
		for( int j = 0 ; j != 3 ; ++j ){
			u_b[i][j] = 0 ;
			r[i][j] = 0 ;
		}
		u_b[i][i] = 12 ;
		r[i][i] = 1 ;
	}
	int d = 12 ;
	int ig = 0 ;
	bool found = true ;
	while( found & ( ig < ng ) ){
		int tr = 0 ;
		for( int i = 0 ; i != 3 ; ++i ){
			tr = tr + uu_g[ig][i][i] ;
		}
		if( tr == 3 ){
			for( int i = 0 ; i != 3 ; ++i ){
				v[i] = 0 ;
				for( int j = 0 ; j != 3 ; ++j ){
					v[i] = v[i] + r[j][i]* u_g[ig][j] ;
				}
			}
			int ib = 0 ;
			int ilen = 3 ;
			yyy_divide_int( ilen, v, d, found ) ;
			while( ! found & ( ib < 3 ) ){
				for( int i = 0 ; i != 3 ; ++i ){
					u[i] = u_b[ib][i] ;
					u_b[ib][i] = u_g[ig][i] ;
				}
				yyy_invert_int( u_b, r, d ) ;
				for( int i = 0 ; i != 3 ; ++i ){
					v[i] = 0 ;
					for( int j = 0 ; j != 3 ; ++j ){
						v[i] = v[i] + r[j][i]* u[j] ;
					}
				}
				int ilen = 3 ;
				yyy_divide_int( ilen, v, d, found ) ;
				if( ! found ){
					for( int i = 0 ; i != 3 ; ++i ){
						u_b[ib][i] = u[i] ;
					}
				}
				ib = ib + 1 ;
			}
		}
		ig = ig + 1 ;
	}
	if( ! found )								throw "yyy_find_base_:a" ;
}
// -----------------------------------------------------------------------------

void yyy_transform_m( int db[3][3], double m[3][3] )
{

	double mdb[3][3] ;

	for( int i = 0 ; i != 3 ; ++i ){
		for( int j = 0 ; j != 3 ; ++j ){
			mdb[i][j] = 0 ;
			for( int k = 0 ; k != 3 ; ++k ){
				mdb[i][j] = mdb[i][j] + m[k][j]* db[i][k] ;
			}
		}
	}
	for( int i = 0 ; i != 3 ; ++i ){
		for( int j = 0 ; j != 3 ; ++j ){
			m[i][j] = 0 ;
			for( int k = 0 ; k != 3 ; ++k ){
				m[i][j] = m[i][j] + db[j][k]* mdb[i][k] ;
			}
		}
	}
}
// -----------------------------------------------------------------------------

void yyy_transform_b( int db[3][3], int b[3][3] )
{

	int bdb[3][3] ;

	for( int i = 0 ; i != 3 ; ++i ){
		for( int j = 0 ; j != 3 ; ++j ){
			bdb[i][j] = 0 ;
			for( int k = 0 ; k != 3 ; ++k ){
				bdb[i][j] = bdb[i][j] + b[k][j]* db[i][k] ;
			}
		}
	}
	for( int i = 0 ; i != 3 ; ++i ){
		for( int j = 0 ; j != 3 ; ++j ){
			b[i][j] = bdb[i][j] ;
		}
	}
}
// -----------------------------------------------------------------------------

void yyy_short_base( double muu[3][3], int uSv[3][3], double mvv[3][3], int vSu[3][3] )
{

	double     sc_t[13] ;
	int        v_b[3][3], b_t[13] ;

	int v_t[26][3] = {
		{  1, 0, 0 },
		{  0, 1, 0 },
		{  0, 0, 1 },
		{  0, 1, 1 },
		{  0, 1,-1 },
		{  1, 0, 1 },
		{ -1, 0, 1 },
		{  1, 1, 0 },
		{  1,-1, 0 },
		{  1, 1, 1 },
		{ -1, 1, 1 },
		{  1,-1, 1 },
		{  1, 1,-1 },
		{  0, 0, 0 },
		{  0, 0, 0 },
		{  0, 0, 0 },
		{  0, 0, 0 },
		{  0, 0, 0 },
		{  0, 0, 0 },
		{  0, 0, 0 },
		{  0, 0, 0 },
		{  0, 0, 0 },
		{  0, 0, 0 },
		{  0, 0, 0 },
		{  0, 0, 0 },
		{  0, 0, 0 },
	} ;

	for( int i = 0 ; i != 3 ; ++i ){
		for( int j = 0 ; j != 3 ; ++j ){
			mvv[i][j] = muu[i][j]/ 144 ;
		}
	}
	for( int it = 13 ; it != 26 ; ++it ){
		for( int i = 0 ; i != 3 ; ++i ){
			v_t[it][i] = - v_t[it-13][i] ;
		}
	}
	yyy_transform_m( uSv, mvv ) ;
	int icy = 0 ;
	bool found = false ;
	while( ! found ){
		icy = icy + 1 ;
		if( icy >= 1000 )						throw "yyy_short_base_:a" ;
		double sc_max = 0 ;
		for( int it = 0 ; it != 13 ; ++it ){
			b_t[it] = 0 ;
			sc_t[it] = 0 ;
			for( int i = 0 ; i != 3 ; ++i ){
				for( int j = 0 ; j != 3 ; ++j ){
					sc_t[it] = sc_t[it] + v_t[it][j]* mvv[i][j]* v_t[it][i] ;
				}
			}
			sc_max = std::max( sc_max, sc_t[it] ) ;
		}

		found = true ;
		for( int ib = 0 ; ib != 3 ; ++ib ){
			double sc_min = sc_max* 2 + 1 ;
			int jt = - 1 ;
			for( int it = 0 ; it != 13 ; ++it ){
				if( ( b_t[it] == 0 ) & ( sc_t[it] < sc_min ) ){
					for( int i = 0 ; i != 3 ; ++i ){
						v_b[ib][i] = v_t[it][i] ;
					}
					bool ok = true ;
					if( ib == 2 ){
						int detb = 0 ;
						yyy_invert_int( v_b, vSu, detb ) ;
						int ilen = 9 ;
						yyy_divide_int( ilen, *vSu, detb, ok ) ;
					}
					if( ok ){
						sc_min = sc_t[it] ;
						jt = it ;
					}
				}
			}
			if( ( jt < 0 ) | ( jt > 12 ) )				throw "yyy_short_base_:b" ;
			b_t[jt] = ib + 1 ;
			for( int i = 0 ; i != 3 ; ++i ){
				v_b[ib][i] = v_t[jt][i] ;
			}
			found = found & ( jt < 3 ) ;
		}
		yyy_transform_b( v_b, uSv ) ;
		yyy_transform_m( v_b, mvv ) ;
	}
	int detx = 0 ;
	yyy_invert_int( uSv, vSu, detx ) ;
	if( ! found | ( detx == 0 ) )						throw "yyy_short_base_:c" ;
	int dety = detx/ 12 ;
	if( detx != dety* 12 )							throw "yyy_short_base_:d" ;
	for( int i = 0 ; i != 3 ; ++i ){
		for( int j = 0 ; j != 3 ; ++j ){
			detx = vSu[i][j]/ dety ;
			if( vSu[i][j] != detx* dety )				throw "yyy_short_base_:e" ;
			vSu[i][j] = detx ;
		}
	}
}
// -----------------------------------------------------------------------------

void yyy_generate_ops( int &mo, int &no, int vv_o[][3][3] )
{

	int v_t[26][3] = {
		{ 1, 0, 0 },
		{ 0, 1, 0 },
		{ 0, 0, 1 },
		{ 0, 1, 1 },
		{ 0, 1,-1 },
		{ 1, 0, 1 },
		{-1, 0, 1 },
		{ 1, 1, 0 },
		{ 1,-1, 0 },
		{ 1, 1, 1 },
		{-1, 1, 1 },
		{ 1,-1, 1 },
		{ 1, 1,-1 },
		{ 0, 0, 0 },
		{ 0, 0, 0 },
		{ 0, 0, 0 },
		{ 0, 0, 0 },
		{ 0, 0, 0 },
		{ 0, 0, 0 },
		{ 0, 0, 0 },
		{ 0, 0, 0 },
		{ 0, 0, 0 },
		{ 0, 0, 0 },
		{ 0, 0, 0 },
		{ 0, 0, 0 },
		{ 0, 0, 0 }
	} ;

	int  w[3][3], ww[3][3], www[3][3], wr[3][3] ;

	for( int it = 0 ; it != 13 ; ++it ){
		for( int i = 0 ; i != 3 ; ++i ){
			v_t[it+13][i] = - v_t[it][i] ;
		}
	}

	no = 0 ;
	for( int it = 0 ; it != 26 ; ++it ){
		for( int jt = 0 ; jt != 26 ; ++jt ){
			for( int kt = 0 ; kt != 26 ; ++kt ){
				for( int i = 0 ; i != 3 ; ++i ){
					w[0][i] = v_t[it][i] ;
					w[1][i] = v_t[jt][i] ;
					w[2][i] = v_t[kt][i] ;
				}

				for( int i = 0 ; i != 3 ; ++i ){
					for( int j = 0 ; j != 3 ; ++j ){
						ww[i][j] = w[i][j] ;
					}
				}
				int det = 0 ;
				yyy_invert_int( w, wr, det ) ;
				int tr = 0 ;
				int trr = 0 ;
				for( int i = 0 ; i != 3 ; ++i ){
					tr = tr + w[i][i] ;
					trr = trr + wr[i][i] ;
				}
				bool found = ( tr == trr ) & ( det == 1 ) ;
				int cou = 0 ;
				while( found & ( cou < tr + 3 ) ){
					cou = cou + 1 ;
					for( int i = 0 ; i != 3 ; ++i ){
						for( int j = 0 ; j != 3 ; ++j ){
							www[i][j] = 0 ;
							for( int k = 0 ; k != 3 ; ++k ){
								www[i][j] = www[i][j] + ww[k][j]* w[i][k] ;
							}
						}
					}
					for( int i = 0 ; i != 3 ; ++i ){
						for( int j = 0 ; j != 3 ; ++j ){
							ww[i][j] = www[i][j] ;
						}
					}
					for( int i = 0 ; i != 3 ; ++i ){
						for( int j = 0 ; j != 3 ; ++j ){
							found = found & ( std::abs( www[i][j] ) <= 1 ) ;
						}
					}
				}
				if( found ){
					for( int i = 0 ; i != 3 ; ++i ){
						for( int j = 0 ; j != 3 ; ++j ){
							vv_o[no][i][j] = w[i][j] ;
						}
					}
					no = no + 1 ;
					if( no > mo )				throw "yyy_generate_ops_:a" ;
				}
			}
		}
	}
	if( no != mo )								throw "yyy_generate_ops_:b" ;
}
// -----------------------------------------------------------------------------

void yyy_eigv3( double w[3][3], double a[3][3], double r[3], double b[3][3] )
{

	double     m[3][3], ab[2][2], abc, maxab ;

	int        c3[5] = { 0, 1, 2, 0, 1 } ;
	double     dp, op, tp, cp, sp ;
	double     dm, om, tm, cm, sm ;
	double     cosa, sina, testa ;
	double     cosb, sinb, testb ;
	double     cosx, sinx, tanx, m2, m3 ;
	double     tmax, tplus ;

	for( int i1 = 0 ; i1 != 3 ; ++i1 ){
		for( int j1 = 0 ; j1 != 3 ; ++j1 ){
			m[i1][j1] = w[i1][j1] ;
			a[i1][j1] = 0 ;
			b[i1][j1] = 0 ;
		}
		a[i1][i1] = 1 ;
		b[i1][i1] = 1 ;
	}
	int i1 = 0 ;
	int icy = 0 ;
	bool more[3] = { true, true, true } ;
	while( more[0] | more[1] | more[2] ){
		if( ++icy > 40 )						throw "yyy_eigv3_:a" ;

		int i2 = c3[i1+1] ;
		int i3 = c3[i1+2] ;
		dp = m[i2][i2] + m[i3][i3] ;
		dm = m[i2][i2] - m[i3][i3] ;
		op = m[i2][i3] + m[i3][i2] ;
		om = m[i2][i3] - m[i3][i2] ;
		tp = sqrt( pow( om, 2 ) + pow( dp, 2 ) ) ;
		tm = sqrt( pow( op, 2 ) + pow( dm, 2 ) ) ;

		// handling zeroes:

		tmax = std::max( tp, tm ) ;

		tplus = tmax + dp ;
		if( tplus == tmax ) dp = 0 ;

		tplus = tmax + dm ;
		if( tplus == tmax ) dm = 0 ;

		tplus = tmax + op ;
		if( tplus == tmax ) op = 0 ;

		tplus = tmax + om ;
		if( tplus == tmax ) om = 0 ;

		tplus = tmax + tp ;
		if( tplus == tmax ) tp = 0 ;

		tplus = tmax + tm ;
		if( tplus == tmax ) tm = 0 ;

		// ordering:

		if( i1 != 1 ) tm = - tm ;

		// sines, cosines:

		if( tm == 0 ){
			cp = 1 ;
			sp = 0 ;
		}
		else{
			cp = dm/ tm ;
			sp = op/ tm ;
		}
		if( tp == 0 ){
			cm = 1 ;
			sm = 0 ;
		}
		else{
			cm = dp/ tp ;
			sm = om/ tp ;
		}

		ab[0][0] = ( cm + cp )/ 2 ;
		ab[1][1] = ( cm - cp )/ 2 ;
		ab[0][1] = ( sp + sm )/ 2 ;
		ab[1][0] = ( sp - sm )/ 2 ;
		m[i2][i2] = ( tp + tm )/ 2 ;
		m[i3][i3] = ( tp - tm )/ 2 ;
		m[i2][i3] = 0 ;
		m[i3][i2] = 0 ;

		int ja = - 1 ;
		int jb = - 1 ;
		maxab = - 1.0 ;
		for( int ib = 0 ; ib != 2 ; ++ib ){
			for( int ia = 0 ; ia != 2 ; ++ia ){
				abc = std::fabs( ab[ib][ia] ) ;
				if( maxab < abc ){
					ja = ia ;
					jb = ib ;
					maxab = abc ;
				}
			}
		}

		if( ( ja == - 1 ) | ( jb == - 1 ) )				throw "yyy_eigv3_:b" ;

		tanx = ab[jb][1-ja]/ ab[jb][ja] ;
		cosx = 1/ sqrt( 1 + pow( tanx, 2 ) ) ;
		sinx = tanx* cosx ;
		if( ab[jb][ja] < 0 ){
			cosx = - cosx ;
			sinx = - sinx ;
		}

		testa = 1 + tanx ;
		if( ja == 0 ){
			cosa = cosx ;
			sina = sinx ;
		}
		else{
			cosa = sinx ;
			sina = cosx ;
		}

		tanx = ab[1-jb][ja]/ ab[jb][ja] ;
		cosx = 1/ sqrt( 1 + pow( tanx, 2 ) ) ;
		sinx = tanx* cosx ;

		testb = 1 + tanx ;
		if( jb == 0 ){
			cosb = cosx ;
			sinb = sinx ;
		}
		else{
			cosb = sinx ;
			sinb = cosx ;
		}

		m2 = m[i1][i2] ;
		m3 = m[i1][i3] ;
		m[i1][i2] = cosa* m2 + sina* m3 ;
		m[i1][i3] = cosa* m3 - sina* m2 ;
		m2 = m[i2][i1] ;
		m3 = m[i3][i1] ;
		m[i2][i1] = cosb* m2 + sinb* m3 ;
		m[i3][i1] = cosb* m3 - sinb* m2 ;

		for( int j1 = 0 ; j1 != 3 ; ++j1 ){
			m2 = a[i2][j1] ;
			m3 = a[i3][j1] ;
			a[i2][j1] = cosa* m2 + sina* m3 ;
			a[i3][j1] = cosa* m3 - sina* m2 ;
			m2 = b[j1][i2] ;
			m3 = b[j1][i3] ;
			b[j1][i2] = cosb* m2 + sinb* m3 ;
			b[j1][i3] = cosb* m3 - sinb* m2 ;
		}
		more[i1] = ( testa != 1 ) | ( testb != 1 ) ;
		i1 = c3[i1+1] ;
	}
	for( int i1 = 0 ; i1 != 3 ; ++i1 ){
		r[i1] = m[i1][i1] ;
	}
//	for( int i1 = 0 ; i1 != 3 ; ++i1 ){
//		for( int j1 = 0 ; j1 != 3 ; ++j1 ){
//			m[i1][j1] = 0 ;
//			for( int k1 = 0 ; k1 != 3 ; ++k1 ){
//				m[i1][j1] = m[i1][j1] + a[k1][j1]* r[k1]* b[i1][k1] ;
//			}
//		}
//	}
}
// -----------------------------------------------------------------------------

void yyy_score_ops( double mvv[3][3], int &no, int vv_o[][3][3], double sc_o[] )
{

	double     a[3][3], r[3], b[3][3] ;

	double     of[3][3], fo[3][3] ;
	double     qo[3][3], oqo[3][3] ;

	double     pi = 4.0* std::atan( 1.0 ) ;

	yyy_eigv3( mvv, a, r, b ) ;

	for( int i = 0 ; i != 3 ; ++i ){
		r[i] = sqrt( r[i] ) ;
	}
	for( int i = 0 ; i != 3 ; ++i ){
		for( int j = 0 ; j != 3 ; ++j ){
			of[i][j] = 0 ;
			for( int k = 0 ; k != 3 ; ++k ){
				of[i][j] = of[i][j] + a[k][j]* r[k]* b[i][k] ;
			}
		}
	}
	for( int i = 0 ; i != 3 ; ++i ){
		r[i] = 1/ r[i] ;
	}
	for( int i = 0 ; i != 3 ; ++i ){
		for( int j = 0 ; j != 3 ; ++j ){
			fo[i][j] = 0 ;
			for( int k = 0 ; k != 3 ; ++k ){
				fo[i][j] = fo[i][j] + a[k][j]* r[k]* b[i][k] ;
			}
		}
	}
	for( int io = 0 ; io != no ; ++io ){
		for( int i = 0 ; i != 3 ; ++i ){
			for( int j = 0 ; j != 3 ; ++j ){
				qo[i][j] = 0 ;
				for( int k = 0 ; k != 3 ; ++k ){
					qo[i][j] = qo[i][j] + vv_o[io][k][j]* fo[i][k] ;
				}
			}
		}
		for( int i = 0 ; i != 3 ; ++i ){
			for( int j = 0 ; j != 3 ; ++j ){
				oqo[i][j] = 0 ;
				for( int k = 0 ; k != 3 ; ++k ){
					oqo[i][j] = oqo[i][j] + of[k][j]* qo[i][k] ;
				}
			}
		}
		yyy_eigv3( oqo, a, r, b ) ;
		sc_o[io] = 0 ;
		for( int i = 0 ; i != 3 ; ++i ){
			sc_o[io] = sc_o[io] + pow( ( 1 - r[i] ), 2 ) ;
		}

//-		score as mean square deviation of bond lengths in the structure with exact twin lattice
//-		sc_o[io] = sqrt( sc_o[io]/ 3 ) ;

//-		score as obliquity angle (in degrees)
		sc_o[io] = 180* atan( sqrt( sc_o[io]/ 2 ) )/ pi ;
	}
}
// -----------------------------------------------------------------------------

void yyy_sort_key( int &nc, int c[], double rc[] )
{

	for( int ic = 0 ; ic != nc ; ++ic ){
		c[ic] = ic ;
		c[ic+nc] = ic ;                           // for nc = 1
	}
	int jc = nc* 2 ;
	int kc = 1 ;
	while( kc < nc ){
		int ia = 0 ;
		int ib = kc ;
		int ic = nc ;
		while( ic < jc ){
			int ja = std::min( nc, ia + kc ) ;
			int jb = std::min( nc, ib + kc ) ;
			while( ( ia < ja ) & ( ib < jb ) ){
				if( rc[c[ib]] - rc[c[ia]] >= 0 ){
					c[ic] = c[ia] ;
					ia = ia + 1 ;
				}
				else{
					c[ic] = c[ib] ;
					ib = ib + 1 ;
				}
				ic = ic + 1 ;
			}
			while( ia < ja ){
				c[ic] = c[ia] ;
				ia = ia + 1 ;
				ic = ic + 1 ;
			}
			while( ib < jb ){
				c[ic] = c[ib] ;
				ib = ib + 1 ;
				ic = ic + 1 ;
			}
			ia = ia + kc ;
			ib = ib + kc ;
		}
		for( int ic = 0 ; ic != nc ; ++ic ){
			c[ic] = c[ic+nc] ;
		}
		kc = kc* 2 ;
	}
}
// -----------------------------------------------------------------------------

void yyy_multiply_ops( int a[3][3], int b[3][3], int c[3][3] )
{

	for( int i = 0 ; i != 3 ; ++i ){
		for( int j = 0 ; j != 3 ; ++j ){
			c[i][j] = 0 ;
			for( int k = 0 ; k != 3 ; ++k ){
				c[i][j] = c[i][j] + a[k][j]* b[i][k] ;
			}
		}
	}
}
// -----------------------------------------------------------------------------

void yyy_the_same_ops( int a[3][3], int b[3][3], bool &found )
{

	found = false ;
	for( int i = 0 ; i != 3 ; ++i ){
		for( int j = 0 ; j != 3 ; ++j ){
			if( a[i][j] != b[i][j] )				return ;
		}
	}
	found = true ;
}
// -----------------------------------------------------------------------------

// Multiplication table is not collected as elements of group
// are collected in order different from sorting order required

void yyy_collect_group(
	int &ivb, double &sc_tol,
	int &no, int h_s[], int o_s[], int vv_o[24][3][3], double sc_o[],
	int &nh, int s_h[], int vv_h[24][3][3], double sc_h[24] )
{

	int vv_x[3][3] ;

	if( ivb >= 3 ){
		for( int i = 0 ; i != 32 ; ++i ) std::cout << '-' ;
		std::cout << std::endl ;
	}

	int ns = no ;
	for( int is = 0 ; is != no ; ++is ){
		h_s[is] = 0 ;
	}

	nh = 0 ;
	int is = 0 ;
	bool found = true ;
	while( found ){
		int ih = nh ;
		int jh = nh ;
		int kh = 0 ;
		h_s[is] = ih + 1 ;
		s_h[ih] = is ;
		while( found & ( jh <= ih ) ){
			yyy_multiply_ops(
				vv_o[o_s[s_h[kh]]], vv_o[o_s[s_h[jh]]], vv_x
			) ;

			int ts = - 1 ;
			found = false ;
			while( ! found & ( ts < no - 1 ) ){
				ts = ts + 1 ;
				yyy_the_same_ops( vv_x, vv_o[o_s[ts]], found ) ;
			}

			if( found & ( h_s[ts] == 0 ) ){
				ih = ih + 1 ;
				s_h[ih] = ts ;
				h_s[ts] = ih + 1 ;
			}

			if( ivb >= 3 ){
				printf( "%4d%2s%3d%1s", ts + 1, "(", h_s[ts], ")" ) ;
				std::cout << std::endl ;
			}

			int th = kh ;
			kh = jh ;
			jh = th ;
			if( kh < jh ){
				kh = kh + 1 ;
			}
			else if( kh == jh ){
				kh = 0 ;
				jh = jh + 1 ;
			}
		}
		if( found ){
			nh = ih + 1 ;
		}
		else{
			for( int jh = nh ; jh <= ih ; jh++ ){
				h_s[s_h[jh]] = 0 ;
			}
		}

		if( ivb >= 3 ){
			for( int i = 0 ; i != 32 ; ++i ) std::cout << '-' ;
			std::cout << std::endl ;
			for( int jh = 0 ; jh <= ih ; jh++ ) printf( "%4d", s_h[jh] + 1 ) ;
			for( int jh = 0 ; jh + 1 <= nh ; jh++ ) printf( "%4d", s_h[jh] + 1 ) ;
			std::cout << std::endl ;
			for( int i = 0 ; i != 32 ; ++i ) std::cout << '-' ;
			std::cout << std::endl ;
		}

		found = false ;
		while( ! found & ( is < ns - 1 ) ){
			is = is + 1 ;
			found = h_s[is] == 0 ;
		}
		found = found & ( sc_o[o_s[is]] <= sc_tol ) ;
	}

	if( nh > 24 )								throw "yyy_collect_group_:a" ;
	int ih = - 1 ;
	for( int is = 0 ; is != no ; ++is ){
		if( h_s[is] != 0 ){
			ih = ih + 1 ;
			if( ih + 1 > nh )					throw "yyy_collect_group_:b" ;
			s_h[ih] = is ;
			sc_h[ih] = sc_o[o_s[is]] ;

			for( int i = 0 ; i != 3 ; ++i ){
				for( int j = 0 ; j != 3 ; ++j ){
					vv_h[ih][i][j] = vv_o[o_s[is]][i][j] ;
				}
			}
		}
	}
	if( ih + 1 != nh )							throw "yyy_collect_group_:c" ;

	if( ivb >= 3 ){
		for( int jh = 0 ; jh + 1 <= nh ; jh++ ) printf( "%4d", s_h[jh] + 1 ) ;
		std::cout << std::endl ;
		for( int i = 0 ; i != 32 ; ++i ) std::cout << '-' ;
		std::cout << std::endl ;
	}
}
// -----------------------------------------------------------------------------

// Produces a group with generators vv_o(io) such as sc_o(io) < sc_tol.
// An operator vv_o(jo) with sc_o(jo) < sc_tol is not present in the group
// if it is not consistent with an operator vv_o(io) and vv_o(io) is present
// in the group and sc_o(io) < sc_o(jo).

void yyy_cell_group(
	int &ivb, double &sc_tol, double mvv[3][3],
	int &nh, int vv_h[24][3][3], double sc_h[24] )
{

	int mo = 504 ;
	int no ;
	int vv_o[mo][3][3] ;
	int o_s[mo*2] ;
	double sc_o[mo] ;

	int h_s[mo] ;
	int s_h[mo] ;

	yyy_generate_ops( mo, no, vv_o ) ;
	yyy_score_ops( mvv, no, vv_o, sc_o ) ;

	sc_o[0] = - 1 ;
	yyy_sort_key( no, o_s, sc_o ) ;
	sc_o[0] = 0 ;

	yyy_collect_group(
		ivb, sc_tol,
		no, h_s, o_s, vv_o, sc_o,
		nh, s_h, vv_h, sc_h
	) ;
}
// -----------------------------------------------------------------------------

// Collects tables of multiplication and reciprocal elements
// thus checking integrity of a group vv_h(,,)

void yyy_test_group( int &nh, int vv_h[24][3][3], int h_hh[24][24], int h_h[24], int tr_h[24] )
{

	int vv_x[3][3] ;
	for( int ih = 0 ; ih != nh ; ++ih ){
		for( int jh = 0 ; jh != nh ; ++jh ){
			yyy_multiply_ops( vv_h[jh], vv_h[ih], vv_x ) ;
			int kh = - 1 ;
			bool found = false ;
			while( ! found & ( kh + 1 < nh ) ){
				kh = kh + 1 ;
				yyy_the_same_ops( vv_x, vv_h[kh], found ) ;
			}
			if( ! found )						throw "yyy_test_group_:a" ;
			h_hh[ih][jh] = kh ;
		}
	}

	bool found = false ;
	if( h_hh[0][0] != 0 )							throw "yyy_test_group_:b" ;
	for( int ih = 0 ; ih != nh ; ++ih ){
		int jh = - 1 ;
		while( ! found & ( jh + 1 < nh ) ){
			jh = jh + 1 ;
			found = h_hh[ih][jh] == 0 ;
		}
		if( ! found )							throw "yyy_test_group_:c" ;
		h_h[ih] = jh ;
		found = false ;
		while( jh + 1 < nh ){
			jh = jh + 1 ;
			if( h_hh[ih][jh] == 0 )					throw "yyy_test_group_:d" ;
		}
	}

	for( int ih = 0 ; ih != nh ; ++ih ){
		if( h_h[h_h[ih]] != ih )					throw "yyy_test_group_:e" ;
	}

	for( int ih = 0 ; ih != nh ; ++ih ){
		tr_h[ih] = 0 ;
		for( int i = 0 ; i != 3 ; ++i ){
			tr_h[ih] = tr_h[ih] + vv_h[ih][i][i] ;
		}
	}
}
// -----------------------------------------------------------------------------

void yyy_write_group( char ch[], int &nh, int vv_h[24][3][3], int h_hh[24][24], int h_h[24], int tr_h[24] )
{

	std::cout << "--------------------------------------" << std::endl ;
	std::cout << ch << std::endl ;
	std::cout << std::endl ;
	std::cout << "    no    tr" << std::endl ;
	for( int ih = 0 ; ih != nh ; ++ih ){
		std::cout << std::endl ;
		for( int j = 0 ; j != 3 ; ++j ){
			if( j == 0 ){
				printf( "%6d%6d", ih + 1, tr_h[ih] ) ;
			}
			else{
				std::cout << "            " ;
			}
			for( int i = 0 ; i != 3 ; ++i ){
				printf( "%4d", vv_h[ih][i][j] ) ;
			}
			std::cout << std::endl ;
		}
	}
	std::cout << std::endl ;
	std::cout << "      " ;
	for( int ih = 0 ; ih != nh ; ++ih ){
		printf( "%3d", ih + 1 ) ;
	}
	std::cout << std::endl ;
	std::cout << std::endl ;
	for( int jh = 0 ; jh != nh ; ++jh ){
		printf( "%3d   ", jh + 1 ) ;
		for( int ih = 0 ; ih != nh ; ++ih ){
			printf( "%3d", h_hh[ih][jh] + 1 ) ;
		}
		std::cout << std::endl ;
	}
	std::cout << std::endl ;
	std::cout << "      " ;
	for( int ih = 0 ; ih != nh ; ++ih ){
		printf( "%3d", ih + 1 ) ;
	}
	std::cout << std::endl ;
	std::cout << std::endl ;
	std::cout << "      " ;
	for( int ih = 0 ; ih != nh ; ++ih ){
		printf( "%3d", h_h[ih] + 1 ) ;
	}
	std::cout << std::endl ;
	std::cout << std::endl ;
}
// -----------------------------------------------------------------------------

void yyy_transform_u_to_v( int uSv[3][3], int vSu[3][3], int &np, int uu_p[24][3][3], int vv_p[24][3][3] )
{

	int vu[3][3] ;

	for( int ip = 0 ; ip != np ; ++ip ){
		yyy_multiply_ops( vSu, uu_p[ip], vu ) ;
		yyy_multiply_ops( vu, uSv, vv_p[ip] ) ;
	}
}
// -----------------------------------------------------------------------------

// Collects tables of multiplication and reciprocal elements
// thus checking integrity of a group vv_h(,,)

void yyy_test_group12( int &nh, int uu_h[24][3][3], int h_hh[24][24], int tr_h[24] )
{

	int uu_x[3][3] ;
	int uu_y[3][3] ;

	for( int ih = 0 ; ih != nh ; ++ih ){
		for( int jh = 0 ; jh != nh ; ++jh ){
			yyy_multiply_ops( uu_h[jh], uu_h[ih], uu_x ) ;
			int kh = h_hh[ih][jh] ;
			for( int i = 0 ; i != 3 ; ++i ){
				for( int j = 0 ; j != 3 ; ++j ){
					uu_y[i][j] = uu_h[kh][i][j]* 12 ;
				}
			}
			bool ok ;
			yyy_the_same_ops( uu_x, uu_y, ok ) ;
			if( ! ok )						throw "yyy_test_group12_:a" ;
		}
	}

	for( int ih = 0 ; ih != nh ; ++ih ){
		int tr_x = 0 ;
		for( int i = 0 ; i != 3 ; ++i ){
			tr_x = tr_x + uu_h[ih][i][i] ;
		}
		if( tr_x != tr_h[ih]* 12 )					throw "yyy_test_group12_:b" ;
	}
}
// -----------------------------------------------------------------------------

void yyy_map_subgroup( int &ng, int uu_g[][3][3], int &nh, int uu_h[24][3][3], int p_h[24], int &np, int h_p[24] )
{

	int uu_x[3][3], uu_y[3][3] ;

	for( int ih = 0 ; ih != nh ; ++ih ){
		p_h[ih] = - 1 ;
	}
	int ip = - 1 ;
	for( int ig = 0 ; ig != ng ; ++ig ){
		for( int i = 0 ; i != 3 ; ++i ){
			for( int j = 0 ; j != 3 ; ++j ){
				uu_x[i][j] = uu_g[ig][i][j]* 12 ;
			}
		}
		int d_x ;
		yyy_invert_int( uu_x, uu_y, d_x ) ;
		if( d_x < 0 ){
			for( int i = 0 ; i != 3 ; ++i ){
				for( int j = 0 ; j != 3 ; ++j ){
					uu_x[i][j] = - uu_x[i][j] ;
				}
			}
		}
		int ih = - 1 ;
		bool ok = false ;
		while( ! ok & ( ih + 1 < nh ) ){
			ih = ih + 1 ;
			yyy_the_same_ops( uu_x, uu_h[ih], ok ) ;
		}
		if( ! ok )							throw "yyy_map_subgroup_:a" ;
		if( p_h[ih] == - 1 ){
			ip = ip + 1 ;
			h_p[ip] = ih ;
			p_h[ih] = ip ;
		}
	}
	np = ip + 1 ;
}
// -----------------------------------------------------------------------------

// input:   nh, h_hh, tr_h, sc_h
//          sc_eps
//          np, h_pc(1:np,1)

// output:  p_h
//          nc, h_c
//          h_pc(1:np,2:nc), nc2

void yyy_right_classes(
	int &nh, int h_hh[24][24], int p_h[24], int tr_h[24], double sc_h[24],
	double &sc_eps,
	int &nc, int h_c[24],
	int &np, int h_pc[24][24], int &nc2 )
{

	double     sc_min, sc_max, sc_pc ;
	int        tr_min, tr_max, tr_pc ;
	int        ip_min, ip_max ;

	for( int ip = 0 ; ip != np ; ++ip ){
		int ih = h_pc[0][ip] ;
		p_h[ih] = ip ;
	}

	int ic = 0 ;
	for( int ih = 0 ; ih != nh ; ++ih ){
		if( p_h[ih] == - 1 ){
			ic = ic + 1 ;
			for( int ip = 0 ; ip != np ; ++ip ){
				int jh = h_hh[ih][h_pc[0][ip]] ;
				if( p_h[jh] != - 1 )				throw "yyy_right_classes_:a" ;
				p_h[jh] = ip ;
				h_pc[ic][ip] = jh ;
			}
		}
	}

	nc = ic + 1 ;
	if( np* nc != nh )							throw "yyy_right_classes_:b" ;
	if( h_pc[0][0] != 0 )							throw "yyy_right_classes_:c" ;
	h_c[0] = 0 ;
	nc2 = 0 ;
	for( int ic = 1 ; ic != nc ; ++ic ){
		tr_min = + 4 ;
		tr_max = - 2 ;
		ip_min = - 1 ;
		ip_max = - 1 ;
		for( int ip = 0 ; ip != np ; ++ip ){
			tr_pc = tr_h[h_pc[ic][ip]] ;
			if( tr_max <= tr_pc ){
				tr_max = tr_pc ;
				ip_max = ip ;
			}
			if( tr_min >= tr_pc ){
				tr_min = tr_pc ;
				ip_min = ip ;
			}
		}
		if( ( tr_min < -1 ) | ( tr_max > 3 ) )				throw "yyy_right_classes_:d" ;
		h_c[ic] = h_pc[ic][ip_min] ;
		if( tr_min == -1 ) nc2 = nc2 + 1 ;
	}

	for( int ic = 0 ; ic != nc ; ++ic ){
		sc_max = 0 ;
		for( int ip = 0 ; ip != np ; ++ip ){
			sc_pc = sc_h[h_pc[ic][ip]] ;
			sc_max = std::max( sc_max, sc_pc ) ;
		}
		sc_min = sc_max ;
		for( int ip = 0 ; ip != np ; ++ip ){
			sc_pc = sc_h[h_pc[ic][ip]] ;
			sc_min = std::min( sc_min, sc_pc ) ;
		}
		if( sc_max - sc_min > sc_eps )					throw "yyy_right_classes_:e" ;
	}
}
// -----------------------------------------------------------------------------

void yyy_write_classes( int &np, int &nc, int h_c[24], int h_pc[24][24], int tr_h[24], double sc_h[24], int &ivb )
{

	if( ivb >= 2 ){
		for( int i = 0 ; i != 78 ; ++i ) std::cout << "-" ;
		std::cout << std::endl << "CG-op vs PG-op\\\\class, " ;
		std::cout << "and selected CG-op and score*1000 vs class:" << std::endl ;

		std::cout << std::endl << "      " ;
		for( int ic = 0 ; ic != nc ; ++ic ){
			printf( "%3d", ic + 1 ) ;
		}
		std::cout << std::endl << std::endl ;
		for( int ip = 0 ; ip != np ; ++ip ){
			printf( "%3d   ", ip + 1 ) ;
			for( int ic = 0 ; ic != nc ; ++ic ){
				printf( "%3d", h_pc[ic][ip] + 1 ) ;
			}
			std::cout << std::endl ;
		}
		std::cout << std::endl << "sel   " ;
		for( int ic = 0 ; ic != nc ; ++ic ){
			printf( "%3d", h_c[ic] + 1 ) ;
		}
		std::cout << std::endl << "sco   " ;
		for( int ic = 0 ; ic != nc ; ++ic ){
//-			printf( "%3d", int( 0.5 + 1000* sc_h[h_c[ic]] ) ) ;
			printf( "%3d", int( 0.5 + 14* sc_h[h_c[ic]] ) ) ;
		}
		std::cout << std::endl ;
	}

	if( ivb >= 1 ){
		std::cout << std::endl ;
		for( int i = 0 ; i != 78 ; ++i ) std::cout << "-" ;
		std::cout << std::endl << "tr(CG-op) vs PG-op\\\\class, " ;
//-		std::cout << "and tr(selected CG-op) and score*1000 vs class:" << std::endl ;
		std::cout << "and tr(selected CG-op) and score*14 vs class:" << std::endl ;

		std::cout << std::endl << "      " ;
		for( int ic = 0 ; ic != nc ; ++ic ){
			printf( "%3d", ic + 1 ) ;
		}
		std::cout << std::endl << std::endl ;
		for( int ip = 0 ; ip != np ; ++ip ){
			printf( "%3d   ", ip + 1 ) ;
			for( int ic = 0 ; ic != nc ; ++ic ){
				printf( "%3d", tr_h[h_pc[ic][ip]] ) ;
			}
			std::cout << std::endl ;
		}
		std::cout << std::endl << "sel   " ;
		for( int ic = 0 ; ic != nc ; ++ic ){
			printf( "%3d", tr_h[h_c[ic]] ) ;
		}
		std::cout << std::endl << "sco   " ;
		for( int ic = 0 ; ic != nc ; ++ic ){
//-			printf( "%3d", int( 0.5 + 1000* sc_h[h_c[ic]] ) ) ;
			printf( "%3d", int( 0.5 + 14* sc_h[h_c[ic]] ) ) ;
		}
		std::cout << std::endl ;
	}

	if( ivb >= 2 ){
		std::cout << std::endl ;
		for( int i = 0 ; i != 78 ; ++i ) std::cout << "-" ;
		std::cout << std::endl ;
	}
}
// -----------------------------------------------------------------------------

//    cell = cell
//    sc_tol = limit acceptable value of: (.05)
//   ( coordinate error on the surface of a sphere )/ ( radius of the sphere )
//         to estimate error, compared are:
//         1) operator that fits translations, but is not orthogonal
//         2) its best fit to an orthogonal operator
//    ng = No of symops
//    uu_g = rotation matrixes ( h(j)* uu_g(j,i)* x(i) !!! )
//    x'(j) = uu_g(j,1)* x(1) + uu_g(j,2)* x(2) + uu_g(j,3)* x(3)
//    u_g = symops' translations (*12)
//    lc = dimension of array for storage of c ( = 48)
//    nc = total number of non-equivalent twin operators
//    nc - 1
//    uu_c = matrices of the above operators ( h(j)* uu_c(j,i)* x(i) !!! ) (*12)
//    trace =-1 --> 2; =0 --> 3; =1 --> 4 ; =2 --> 6 ; =3 --> [1]
//    sc_c = scores for operators e ( in the above sense of errors on the
//    surface of the sphere), i.e. sc_c < sc_tol except those operators
//    from the factor group, that are generated from the low score
//    generators to close the factor group.

//    ivb == 0 : no output,
//    ivb == 1 : a table of traces of cell-group operators (vs SG-ops and classses)
//    ivb == 2 : same and operators of cell-group etc.
//    ivb == 3 : more

//    ivb > 0 : print errors
//    ierr == 0 : OK
//    ierr == 1 : error

//-   13/01/2012 ------
//-   new_score (obliquity) = 70.175* old_score
//-   (see commented lines at the end of "subroutine yyy_score_ops")
//-   e.g. score 0.03 approximately corresponds to obliquity 2 degrees

void yyy_cell2tg(
	double cell[6], double &sc_tol,
	int &ng, int uu_g[][3][3], int u_g[][3],
	int &lc, int &nc, int &nc2, int uu_c[][3][3], double sc_c[],
	int &ivb, int &ierr )
{

//-	double     sc_eps = 1.0e-04 ;
	double     sc_eps = 7.0e-03 ;

	double     muu[3][3], mvv[3][3] ;
	double     celt[6] ;
	int        uSv[3][3], vSu[3][3] ;

	int        nh ;
	int        vv_h[24][3][3] ;
	int        uu_h[24][3][3] ;
	int        h_hh[24][24] ;
	int        h_h[24] ;
	int        tr_h[24] ;
	int        p_h[24] ;
	double     sc_h[24] ;

	int        np ;
	int        h_c[24] ;
	int        h_pc[24][24] ;

	try{
		yyy_cell2met( cell, muu ) ;
		yyy_met2cell( celt, muu ) ;
		yyy_find_base( ng, uu_g, u_g, uSv ) ;
		yyy_short_base( muu, uSv, mvv, vSu ) ;

		yyy_cell_group( ivb, sc_tol, mvv, nh, vv_h, sc_h ) ;
		yyy_test_group( nh, vv_h, h_hh, h_h, tr_h ) ;

		if( ivb >= 2 ){
			char ch[] = { "vv_h" } ;
			yyy_write_group( ch, nh, vv_h, h_hh, h_h, tr_h ) ;
		}

		yyy_transform_u_to_v( vSu, uSv, nh, vv_h, uu_h ) ;
		yyy_test_group12( nh, uu_h, h_hh, tr_h ) ;

		if( ivb >= 2 ){
			char ch[] = { "uu_h" } ;
			yyy_write_group( ch, nh, uu_h, h_hh, h_h, tr_h ) ;
		}

		yyy_map_subgroup( ng, uu_g, nh, uu_h, p_h, np, *h_pc ) ;
		yyy_right_classes(
			nh, h_hh, p_h, tr_h, sc_h, sc_eps,
			nc, h_c,
			np, h_pc, nc2
		) ;

		if( ivb >= 1 ){
			yyy_write_classes( np, nc, h_c, h_pc, tr_h, sc_h, ivb ) ;
		}

		if( nc > lc )							throw "yyy_cell2tg_:a" ;
		for( int ic = 0 ; ic != nc ; ++ic ){
			sc_c[ic] = sc_h[h_c[ic]] ;
			for( int i = 0 ; i != 3 ; ++i ){
				for( int j = 0 ; j != 3 ; ++j ){
					uu_c[ic][i][j] = uu_h[h_c[ic]][i][j] ;
				}
			}
		}
		ierr = 0 ;
	}
	catch( const char *msg ){
		ierr = 1 ;
		if( ivb > 0 ){
			std::cerr << msg << std::endl ;
			std::cout << msg << std::endl ;
		}
	}
	catch( ... ){
		ierr = 1 ;
		if( ivb > 0 ){
			std::cerr << "yyy_cell2tg_:b" << std::endl ;
			std::cout << "yyy_cell2tg_:b" << std::endl ;
		}
	}
}
// -----------------------------------------------------------------------------

