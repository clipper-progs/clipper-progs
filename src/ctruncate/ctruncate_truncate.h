//
//     CTRUNCATE
//     Copyright (C) 2006-2011 Norman Stein, Charles Ballard
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//


#ifndef __CTRUNCATE_TRUNCATE_H
#define __CTRUNCATE_TRUNCATE_H

#include "alt_hkl_datatypes.h"
#include "clipper/clipper.h"
#include "clipper/clipper-ccp4.h"

namespace ctruncate {
	
	/*int truncate(clipper::HKL_data<clipper::data32::I_sigI>& isig, clipper::HKL_data<clipper::data32::I_sigI>& jsig,
				 clipper::HKL_data<clipper::data32::F_sigF>& fsig, clipper::ResolutionFn& Sigma, float scalef, 
				 clipper::Resolution& reso, int& nrej, bool debug);
	int truncate(clipper::HKL_data<clipper::data32::J_sigJ_ano>& isig, clipper::HKL_data<clipper::data32::J_sigJ_ano>& jsig, 
				 clipper::HKL_data<clipper::data32::G_sigG_ano>& fsig, clipper::ResolutionFn& Sigma, float scalef, 
				 clipper::Resolution& reso, int& nrej, bool debug);*/
	template<class T1, class T2, class T3, class T4> int truncate(const clipper::HKL_data<clipper::datatypes::I_sigI<T1> >& isig, clipper::HKL_data<clipper::datatypes::I_sigI<T2> >& jsig,
				 clipper::HKL_data<clipper::datatypes::F_sigF<T3> >& fsig, const clipper::HKL_data<clipper::datatypes::I_sigI<T4> >& Sigma, float scalef, 
				 clipper::Resolution& reso, int& nrej, bool debug);
	template<class T1, class T2, class T3, class T4> int truncate(const clipper::HKL_data<clipper::datatypes::I_sigI_ano<T1> >& isig, clipper::HKL_data<clipper::datatypes::I_sigI_ano<T2> >& jsig, 
				 clipper::HKL_data<clipper::datatypes::F_sigF_ano<T3> >& fsig, const clipper::HKL_data<clipper::datatypes::I_sigI<T4> >& Sigma, float scalef, 
				 clipper::Resolution& reso, int& nrej, bool debug);

	template<class T1, class T2, class T3> int truncate(clipper::HKL_data<clipper::datatypes::I_sigI<T1> >& isig, clipper::HKL_data<clipper::datatypes::I_sigI<T2> >& jsig,
				 clipper::HKL_data<clipper::datatypes::F_sigF<T3> >& fsig, float scalef, 
				 clipper::Resolution& reso, int& nrej, bool debug);
	template<class T1, class T2, class T3> int truncate(clipper::HKL_data<clipper::datatypes::I_sigI_ano<T1> >& isig, clipper::HKL_data<clipper::datatypes::I_sigI_ano<T2> >& jsig, 
				 clipper::HKL_data<clipper::datatypes::F_sigF_ano<T3> >& fsig, float scalef, 
				 clipper::Resolution& reso, int& nrej, bool debug);
	
	template<class T1, class T2, class T3> int truncate_sivia(clipper::HKL_data<clipper::datatypes::I_sigI<T1> >& isig, clipper::HKL_data<clipper::datatypes::I_sigI<T2> >& jsig,
				 clipper::HKL_data<clipper::datatypes::F_sigF<T3> >& fsig, float scalef, 
				 clipper::Resolution& reso, int& nrej, bool debug);
	template<class T1, class T2, class T3> int truncate_sivia(clipper::HKL_data<clipper::datatypes::I_sigI_ano<T1> >& isig, clipper::HKL_data<clipper::datatypes::I_sigI_ano<T2> >& jsig, 
				 clipper::HKL_data<clipper::datatypes::F_sigF_ano<T3> >& fsig, float scalef, 
				 clipper::Resolution& reso, int& nrej, bool debug);
	
	int truncate_centric(float I, float sigma, float S, float& J, float& sigJ, float& F, float& sigF, int& nrej, bool debug);
	int truncate_acentric(float I, float sigma, float S, float& J, float& sigJ, float& F, float& sigF, int& nrej, bool debug);
	int truncate_flat(float I, float sigma, float& J, float& sigJ, float& F, float& sigF, int& nrej, bool debug);
	int truncate_sivia_calc(float I, float sigma, float& J, float& sigJ, float& F, float& sigF, int& nrej, bool debug);
	
	
	// Wilson truncate: takes I's as input
	
	 template<class T1, class T2, class T3, class T4> int truncate(const clipper::HKL_data<clipper::datatypes::I_sigI<T1> >& isig, clipper::HKL_data<clipper::datatypes::I_sigI<T2> >& jsig,
											 clipper::HKL_data<clipper::datatypes::F_sigF<T3> >& fsig, const clipper::HKL_data<clipper::datatypes::I_sigI<T4> >& Sigma, float scalef, 
											 clipper::Resolution& reso, int& nrej, bool debug)
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		float J, sigJ, F, sigF;
		int iflag;
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() && ih.invresolsq() <=  reso.invresolsq_limit() ) {
				float I = isig[ih].I();
				float sigma = isig[ih].sigI();
				float S = Sigma[ih].I();
				clipper::HKL hkl = ih.hkl();
				float weight(ih.hkl_class().epsilon() );
				float sqwt = sqrt(weight);
				
				I /= weight;
				sigma /= weight;
				
				// handle acentric and centric reflections separately
				if ( ih.hkl_class().centric() ) iflag = truncate_centric(I, sigma, S, J, sigJ, F, sigF, nrej, debug);
				else iflag = truncate_acentric(I, sigma, S, J, sigJ, F, sigF, nrej, debug);					if (iflag) {
					jsig[ih].I() = J*weight;
					jsig[ih].sigI() = sigJ*weight;
					fsig[ih].f() = F*scalef*sqwt;
					fsig[ih].sigf() = sigF*scalef*sqwt;
				} else {
					jsig[ih].I() = clipper::Util::nan();
					jsig[ih].sigI() = clipper::Util::nan();
					fsig[ih].f() = clipper::Util::nan();
					fsig[ih].sigf() = clipper::Util::nan();
				}
			} else {
				jsig[ih].I() = clipper::Util::nan();
				jsig[ih].sigI() = clipper::Util::nan();
				fsig[ih].f() = clipper::Util::nan();
				fsig[ih].sigf() = clipper::Util::nan();
			}
		}
		return(1);
	}
	
	// Wilson truncate: takes anomalous I's as input. 
	
	template<class T1, class T2, class T3, class T4> int truncate(const clipper::HKL_data<clipper::datatypes::I_sigI_ano<T1> >& isig, clipper::HKL_data<clipper::datatypes::I_sigI_ano<T2> >& jsig, 
																  clipper::HKL_data<clipper::datatypes::F_sigF_ano<T3> >& fsig, const clipper::HKL_data<clipper::datatypes::I_sigI<T4> >& Sigma, float scalef, 
																  clipper::Resolution& reso, int& nrej, bool debug)
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		float J, sigJ, F, sigF;
		int iflag;
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( ih.invresolsq() <=  reso.invresolsq_limit() ) {
				if ( !clipper::Util::is_nan(isig[ih].I_pl() ) ) {
					float I = isig[ih].I_pl();
					float sigma = isig[ih].sigI_pl();
					float S = Sigma[ih].I();
					clipper::HKL hkl = ih.hkl();
					float weight(ih.hkl_class().epsilon() );
					float sqwt = sqrt(weight);
					
					I /= weight;
					sigma /= weight;
					
					// handle acentric and centric reflections separately
					if ( ih.hkl_class().centric() ) iflag = truncate_centric(I, sigma, S, J, sigJ, F, sigF, nrej, debug);
					else iflag = truncate_acentric(I, sigma, S, J, sigJ, F, sigF, nrej, debug);		 if (iflag) {
						jsig[ih].I_pl() = J*weight;
						jsig[ih].sigI_pl() = sigJ*weight;
						fsig[ih].f_pl() = F*scalef*sqwt;
						fsig[ih].sigf_pl() = sigF*scalef*sqwt;
					} else {
						jsig[ih].I_pl() = clipper::Util::nan();
						jsig[ih].sigI_pl() = clipper::Util::nan();
						fsig[ih].f_pl() = clipper::Util::nan();
						fsig[ih].sigf_pl() = clipper::Util::nan();
					}
				} else {
					jsig[ih].I_pl() = clipper::Util::nan();
					jsig[ih].sigI_pl() = clipper::Util::nan();
					fsig[ih].f_pl() = clipper::Util::nan();
					fsig[ih].sigf_pl() = clipper::Util::nan();
				}
				
				if ( !clipper::Util::is_nan(isig[ih].I_mi() ) ) {
					float I = isig[ih].I_mi();
					float sigma = isig[ih].sigI_mi();
					float S = Sigma[ih].I();
					clipper::HKL hkl = ih.hkl();
					float weight(ih.hkl_class().epsilon() );
					float sqwt = sqrt(weight);
					
					I /= weight;
					sigma /= weight;
					
					// handle acentric and centric reflections separately
					if ( ih.hkl_class().centric() ) iflag = truncate_centric(I, sigma, S, J, sigJ, F, sigF, nrej, debug);
					else iflag = truncate_acentric(I, sigma, S, J, sigJ, F, sigF, nrej, debug);	
					if (iflag) {
						jsig[ih].I_mi() = J*weight;
						fsig[ih].f_mi() = F*scalef*sqwt;
						fsig[ih].sigf_mi() = sigF*scalef*sqwt;
					} else {
						jsig[ih].I_mi() = clipper::Util::nan();
						jsig[ih].sigI_mi() = clipper::Util::nan();
						fsig[ih].f_mi() = clipper::Util::nan();
						fsig[ih].sigf_mi() = clipper::Util::nan();
					}
				} else {
					jsig[ih].I_mi() = clipper::Util::nan();
					jsig[ih].sigI_mi() = clipper::Util::nan();
					fsig[ih].f_mi() = clipper::Util::nan();
					fsig[ih].sigf_mi() = clipper::Util::nan();
				}
			} else {
				jsig[ih].I_pl() = clipper::Util::nan();
				jsig[ih].sigI_pl() = clipper::Util::nan();
				fsig[ih].f_pl() = clipper::Util::nan();
				fsig[ih].sigf_pl() = clipper::Util::nan();
				jsig[ih].I_mi() = clipper::Util::nan();
				jsig[ih].sigI_mi() = clipper::Util::nan();
				fsig[ih].f_mi() = clipper::Util::nan();
				fsig[ih].sigf_mi() = clipper::Util::nan();
			}
		}
		return(1);
	}
	
	// flat prior based on means
	
	template<class T1, class T2, class T3> int truncate(clipper::HKL_data<clipper::datatypes::I_sigI<T1> >& isig, clipper::HKL_data<clipper::datatypes::I_sigI<T2> >& jsig,
														clipper::HKL_data<clipper::datatypes::F_sigF<T3> >& fsig, float scalef, 
														clipper::Resolution& reso, int& nrej, bool debug)
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		float J, sigJ, F, sigF;
		int iflag;
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() && ih.invresolsq() <=  reso.invresolsq_limit() ) {
				float I = isig[ih].I();
				float sigma = isig[ih].sigI();
				clipper::HKL hkl = ih.hkl();
				float weight(ih.hkl_class().epsilon() );
				float sqwt = sqrt(weight);
				
				I /= weight;
				sigma /= weight;
				
				iflag = truncate_flat(I, sigma, J, sigJ, F, sigF, nrej, debug);
				if (iflag) {
					jsig[ih].I() = J*weight;
					jsig[ih].sigI() = sigJ*weight;
					fsig[ih].f() = F*scalef*sqwt;
					fsig[ih].sigf() = sigF*scalef*sqwt;
				} else {
					jsig[ih].I() = clipper::Util::nan();
					jsig[ih].sigI() = clipper::Util::nan();
					fsig[ih].f() = clipper::Util::nan();
					fsig[ih].sigf() = clipper::Util::nan();
				}
			} else {
				jsig[ih].I() = clipper::Util::nan();
				jsig[ih].sigI() = clipper::Util::nan();
				fsig[ih].f() = clipper::Util::nan();
				fsig[ih].sigf() = clipper::Util::nan();
			}
		}
		return(1);
	}
	
	// flat prior based on means, input as anomalous Is
	
	template<class T1, class T2, class T3> int truncate(clipper::HKL_data<clipper::datatypes::I_sigI_ano<T1> >& isig, clipper::HKL_data<clipper::datatypes::I_sigI_ano<T2> >& jsig, 
														clipper::HKL_data<clipper::datatypes::F_sigF_ano<T3> >& fsig, float scalef, 
														clipper::Resolution& reso, int& nrej, bool debug)
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		float J, sigJ, F, sigF;
		float invr2 = reso.invresolsq_limit();
		int iflag;
		
        for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
            if ( ih.invresolsq() <=  reso.invresolsq_limit() ) {
                if ( !clipper::Util::is_nan(isig[ih].I_pl() ) ) {
                    float I = isig[ih].I_pl();
                    float sigma = isig[ih].sigI_pl();
                    clipper::HKL hkl = ih.hkl();
                    float weight(ih.hkl_class().epsilon() );
                    float sqwt = sqrt(weight);
                    
                    I /= weight;
                    sigma /= weight;
                    
                    iflag = truncate_flat(I, sigma, J, sigJ, F, sigF, nrej, debug);
                    if (iflag) {
                        jsig[ih].I_pl() = J*weight;
                        jsig[ih].sigI_pl() = sigJ*weight;
                        fsig[ih].f_pl() = F*scalef*sqwt;
                        fsig[ih].sigf_pl() = sigF*scalef*sqwt;
                    } else {
                        jsig[ih].I_pl() = clipper::Util::nan();
                        jsig[ih].sigI_pl() = clipper::Util::nan();
                        fsig[ih].f_pl() = clipper::Util::nan();
                        fsig[ih].sigf_pl() = clipper::Util::nan();
                    }
                } else {
                    jsig[ih].I_pl() = clipper::Util::nan();
                    jsig[ih].sigI_pl() = clipper::Util::nan();
                    fsig[ih].f_pl() = clipper::Util::nan();
                    fsig[ih].sigf_pl() = clipper::Util::nan();
                }
                
                if ( !clipper::Util::is_nan(isig[ih].I_mi() ) ) {
                    float I = isig[ih].I_mi();
                    float sigma = isig[ih].sigI_mi();
                    clipper::HKL hkl = ih.hkl();
                    float weight(ih.hkl_class().epsilon() );
                    float sqwt = sqrt(weight);
                    
                    I /= weight;
                    sigma /= weight;
                    
                    iflag = truncate_flat(I, sigma, J, sigJ, F, sigF, nrej, debug);
                    if (iflag) {
                        jsig[ih].I_mi() = J*weight;
                        jsig[ih].sigI_mi() = sigJ*weight;
                        fsig[ih].f_mi() = F*scalef*sqwt;
                        fsig[ih].sigf_mi() = sigF*scalef*sqwt;
                    } else {
                        jsig[ih].I_mi() = clipper::Util::nan();
                        jsig[ih].sigI_mi() = clipper::Util::nan();
                        fsig[ih].f_mi() = clipper::Util::nan();
                        fsig[ih].sigf_mi() = clipper::Util::nan();
                    }
                } else {
                    jsig[ih].I_mi() = clipper::Util::nan();
                    jsig[ih].sigI_mi() = clipper::Util::nan();
                    fsig[ih].f_mi() = clipper::Util::nan();
                    fsig[ih].sigf_mi() = clipper::Util::nan();
                }
            } else {
                jsig[ih].I_pl() = clipper::Util::nan();
                jsig[ih].sigI_pl() = clipper::Util::nan();
                fsig[ih].f_pl() = clipper::Util::nan();
                fsig[ih].sigf_pl() = clipper::Util::nan();
                jsig[ih].I_mi() = clipper::Util::nan();
                jsig[ih].sigI_mi() = clipper::Util::nan();
                fsig[ih].f_mi() = clipper::Util::nan();
                fsig[ih].sigf_mi() = clipper::Util::nan();
            }
        }
		return(1);
	}
	
	// flat prior of Sivia, based on arithmetic solution of second order approximation
	
	template<class T1, class T2, class T3> int truncate_sivia(clipper::HKL_data<clipper::datatypes::I_sigI<T1> >& isig, clipper::HKL_data<clipper::datatypes::I_sigI<T2> >& jsig,
															  clipper::HKL_data<clipper::datatypes::F_sigF<T3> >& fsig, float scalef, 
															  clipper::Resolution& reso, int& nrej, bool debug)
	{
		typedef clipper::HKL_data_base::HKL_reference_index HRI;
		float J, sigJ, F, sigF;
		int iflag;
		
		for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
			if ( !isig[ih].missing() && ih.invresolsq() <=  reso.invresolsq_limit() ) {
				float I = isig[ih].I();
				float sigma = isig[ih].sigI();
				clipper::HKL hkl = ih.hkl();
				float weight(ih.hkl_class().epsilon() );
				float sqwt = sqrt(weight);
				
				I /= weight;
				sigma /= weight;
				
				iflag = truncate_sivia_calc(I, sigma, J, sigJ, F, sigF, nrej, debug);
				if (iflag) {
					jsig[ih].I() = J*weight;
					jsig[ih].sigI() = sigJ*weight;
					fsig[ih].f() = F*scalef*sqwt;
					fsig[ih].sigf() = sigF*scalef*sqwt;
				} else {
					jsig[ih].I() = clipper::Util::nan();
					jsig[ih].sigI() = clipper::Util::nan();
					fsig[ih].f() = clipper::Util::nan();
					fsig[ih].sigf() = clipper::Util::nan();
				}
			} else {
				jsig[ih].I() = clipper::Util::nan();
				jsig[ih].sigI() = clipper::Util::nan();
				fsig[ih].f() = clipper::Util::nan();
				fsig[ih].sigf() = clipper::Util::nan();
			}
		}
		return(1);
	}
	
	// flat prior of Sivia, based on arithmetic solution of second order approximation, anomalous input
	
    template<class T1, class T2, class T3> int truncate_sivia(clipper::HKL_data<clipper::datatypes::I_sigI_ano<T1> >& isig, clipper::HKL_data<clipper::datatypes::I_sigI_ano<T2> >& jsig,
                                                              clipper::HKL_data<clipper::datatypes::F_sigF_ano<T3> >& fsig, float scalef,
                                                              clipper::Resolution& reso, int& nrej, bool debug)
    {
        typedef clipper::HKL_data_base::HKL_reference_index HRI;
        float J, sigJ, F, sigF;
        float invr2 = reso.invresolsq_limit();
        int iflag;
        
        for ( HRI ih = isig.first(); !ih.last(); ih.next() ) {
            if ( ih.invresolsq() <=  reso.invresolsq_limit() ) {
                if ( !clipper::Util::is_nan(isig[ih].I_pl() ) ) {
                    float I = isig[ih].I_pl();
                    float sigma = isig[ih].sigI_pl();
                    clipper::HKL hkl = ih.hkl();
                    float weight(ih.hkl_class().epsilon() );
                    float sqwt = sqrt(weight);
                    
                    I /= weight;
                    sigma /= weight;
                    
                    iflag = truncate_sivia_calc(I, sigma, J, sigJ, F, sigF, nrej, debug);
                    if (iflag) {
                        jsig[ih].I_pl() = J*weight;
                        jsig[ih].sigI_pl() = sigJ*weight;
                        fsig[ih].f_pl() = F*scalef*sqwt;
                        fsig[ih].sigf_pl() = sigF*scalef*sqwt;
                    } else {
                        jsig[ih].I_pl() = clipper::Util::nan();
                        jsig[ih].sigI_pl() = clipper::Util::nan();
                        fsig[ih].f_pl() = clipper::Util::nan();
                        fsig[ih].sigf_pl() = clipper::Util::nan();
                    }
                } else {
                    jsig[ih].I_pl() = clipper::Util::nan();
                    jsig[ih].sigI_pl() = clipper::Util::nan();
                    fsig[ih].f_pl() = clipper::Util::nan();
                    fsig[ih].sigf_pl() = clipper::Util::nan();
                }
                
                if ( !clipper::Util::is_nan(isig[ih].I_mi() ) ) {
                    float I = isig[ih].I_mi();
                    float sigma = isig[ih].sigI_mi();
                    clipper::HKL hkl = ih.hkl();
                    float weight(ih.hkl_class().epsilon() );
                    float sqwt = sqrt(weight);
                    
                    I /= weight;
                    sigma /= weight;
                    
                    iflag = truncate_sivia_calc(I, sigma, J, sigJ, F, sigF, nrej, debug);
                    if (iflag) {
                        jsig[ih].I_mi() = J*weight;
                        jsig[ih].sigI_mi() = sigJ*weight;
                        fsig[ih].f_mi() = F*scalef*sqwt;
                        fsig[ih].sigf_mi() = sigF*scalef*sqwt;
                    } else {
                        jsig[ih].I_mi() = clipper::Util::nan();
                        jsig[ih].sigI_mi() = clipper::Util::nan();
                        fsig[ih].f_mi() = clipper::Util::nan();
                        fsig[ih].sigf_mi() = clipper::Util::nan();
                    }
                } else {
                    jsig[ih].I_mi() = clipper::Util::nan();
                    jsig[ih].sigI_mi() = clipper::Util::nan();
                    fsig[ih].f_mi() = clipper::Util::nan();
                    fsig[ih].sigf_mi() = clipper::Util::nan();
                }
            } else {
                jsig[ih].I_pl() = clipper::Util::nan();
                jsig[ih].sigI_pl() = clipper::Util::nan();
                fsig[ih].f_pl() = clipper::Util::nan();
                fsig[ih].sigf_pl() = clipper::Util::nan();
                jsig[ih].I_mi() = clipper::Util::nan();
                jsig[ih].sigI_mi() = clipper::Util::nan();
                fsig[ih].f_mi() = clipper::Util::nan();
                fsig[ih].sigf_mi() = clipper::Util::nan();
            }
        }
        return(1);
    }
	
}

#endif
