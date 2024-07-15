#ifndef CONSTANT_HPP
#define CONSTANT_HPP

#if __cplusplus < 201402L
    #pragma message("Error : Requires C++14 Support !")
#endif
#ifdef __CUDA_ARCH__
    #define NAN_REAL CUDART_NAN
    #define INF_REAL CUDART_IN
    #include <math_constants.h>
#else
    template <typename Real>
    const Real INF_REAL= std::numeric_limits<Real>::infinity();

    template <typename Real>
    const Real NAN_REAL= std::numeric_limits<Real>::quiet_NaN();
#endif

#ifdef __MINGW32__
    #warning "Using MinGW compiler"
#elif _MSC_VER
    #warning "Using Visual C++ compiler"
#endif
    
#include <limits>

namespace constant {

#ifndef PI
    template <typename Real>
    const Real PI = Real(3.1415);
#endif
#ifndef GRAVITY
    template <typename Real>  
    const Real GRAVITY = Real(9.81) ;
#endif
#ifndef TIME_STEP
    template <typename Real>
    const Real TIME_STEP = Real(0.001) ;
#endif
#ifndef EPSILON
    template <typename Real>
    const Real EPSILON = std::numeric_limits<Real>::epsilon();
#endif
#ifndef DIFF_STEP
    template <typename Real>
    const Real DIFF_STEP = Real(0.001);
#endif





}; // namespace constant
#endif // CONSTANTS_HPP