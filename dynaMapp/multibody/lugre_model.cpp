#include "lugre_model.hpp"
 
namespace friction {
  
template <typename Real, typename Int>
Lugre<Real, Int>::Lugre(Real sigma0_, Real sigma1_, Real sigma2_, 
                        Real vs_, Real coulomb_, Real stiction_)
{
    this->coulomb = coulomb_;
    this->stiction = stiction_;
    this->vs = vs_;
    this->sigma0 = sigma0_;
    this->sigma1 = sigma1_;
    this->sigma2 = sigma2_;
    this->initialized = false;
}

template <typename Real, typename Int>
void Lugre<Real, Int>::init()
{
    this->t0     = 0.0;
    this->z0     = 0.0;
    this->z_dot0 = 0.0;
    this->z      = z0;
    this->t      = t0;
    this->z_dot  = z_dot0;
    this->initialized = true;
}

template <typename Real, typename Int>
void Lugre<Real, Int>::init(Real t0, Real z0, Real z_dot0)
{
    this->t0 = t0;
    this->z0 = z0;
    this->z_dot0 = z_dot0;

    this->t = t0;
    this->z = z0;
    this->z_dot = z_dot0;

    this->initialized = true;
}
template <typename Real, typename Int>
Real Lugre<Real, Int>::friction_torque(Real velocity, Real t) noexcept{   
    if (!this->initialized)
    {
        throw std::runtime_error("Lugre model is not initialized!");
    }
    Real torque ;
    while (this->t<=t)
    {   
        torque= this->sigma0*this->z+this->sigma1*z_dot+sigma2*velocity;
        Lugre::update_state(velocity);
        this->t += constant::TIME_STEP<Real>;
    }
    return  torque;
}

template <typename Real, typename Int>
inline void Lugre<Real, Int>::update_state(Real velocity) noexcept
{   
    if (!this->initialized)
    {
        throw std::runtime_error("Lugre model is not initialized!");
    }
    Real g_v = Lugre<Real, Int>::stribek_force(velocity);
    this->z_dot=velocity-this->sigma0*(abs(this->velocity))/g_v*this->z;
    this->z = this->z + constant::TIME_STEP<Real> * z_dot;
}

template <typename Real, typename Int>
inline Real Lugre<Real, Int>::stribek_force(Real velocity)
{
    Real r = velocity/this->vs;
    return this->coulomb_force+(this->stiction_force-this->coulomb_force)
        *exp(-pow(r,2));
}

template <typename Real, typename Int>
inline Real Lugre<Real, Int>::steady_torque(Real velocity)
{

}

template <typename Real, typename Int>
inline Real Lugre<Real, Int>::steady_state(Real velocity)
{

}


template<typename Real, typename Int>
LugreTemp<Real, Int>::LugreTemp(Real Te, Real sigma0, Real sigma0e, Real sigma1,
                                Real sigma1e,Real sigma2, Real sigma2e,Real vs,
                                Real coulomb, Real stiction, Real gradT)
{
    this->Te=Te;
    this->stiction=stiction;
    this->gradT=gradT;
    this->coulomb=coulomb;
    this->vs=vs;
    this->sigma0 = sigma0;
    this->sigma1 = sigma1;
    this->sigma2 = sigma2;
    this->initialized = false ;
}

template <typename Real, typename Int>
inline void LugreTemp<Real, Int>::init()
{
    this->t0     = 0.0;
    this->z0     = 0.0;
    this->z_dot0 = 0.0;
    this->z      = z0;
    this->t      = t0;
    this->z_dot  = z_dot0;
    this->initialized = true;
}

template <typename Real, typename Int>
inline void LugreTemp<Real, Int>::init(Real Ti, Real t0, Real z0, Real z_dot0)
{
    this->t0     = t0;
    this->z0     = z0;
    this->z_dot0 = z_dot0;

    this->t     = t0;
    this->z     = z0;
    this->z_dot = z_dot0;
    this->Ti    = Ti;

    this->initialized = true;

}

template <typename Real, typename Int>
Real LugreTemp<Real, Int>::friction_torque(Real velocity, Real t)
{
    return Real();
}

template <typename Real, typename Int>
inline Real LugreTemp<Real, Int>::steady_state(Real velocity)
{
    return Real();
}

}; //namespace friction 