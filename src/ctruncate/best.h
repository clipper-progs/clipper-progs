//
//     CTRUNCATE best routine
//     Copyright (C) Airlie McCoy
//     Originally best.h in phaser
//
//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Application.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//

#include "clipper/core/clipper_util.h"
namespace ctruncate {

    class Best {
    public:
        //contained within best range
        static bool contains(const clipper::ftype);
        //return value
        static clipper::ftype value(const clipper::ftype);
        //maximum invresolsq
        static clipper::ftype invresolsq_max();
        //minimum invresolsq
        static clipper::ftype invresolsq_min();
        
    private:
        static const double k = 8.563714792513845e-06;         // magic scale factor to put best data on scale for average electron squared
        static const double offset = 0.009; // low resolution shell at 10.5409 A resolution
        static const double step = 0.004091973;
    };

    float BEST(float);

}
