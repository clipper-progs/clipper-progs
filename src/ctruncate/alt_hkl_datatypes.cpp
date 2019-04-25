/*! \file alt_hkl_datatypes.cpp
    alternative datatypes for mtz output for the clipper libraries
*/


#include "alt_hkl_datatypes.h"

namespace clipper {

namespace datatypes {

// compile template types

template class J_sigJ_ano<ftype32>;
template class G_sigG_ano<ftype32>;

template class J_sigJ_ano<ftype64>;
template class G_sigG_ano<ftype64>;

} // namespace datatypes

} // namespace clipper
