#ifndef LUGRE_HPP
#define LUGRE_HPP

 
#include <stdexcept>
#include <cmath>

#include "utils/constant.hpp"

namespace friction { 

/**
 * @headerfile  LuGre friction model class declaration. 
 * @copyright (c) 2024  Dynamium
 * @date   22-March-2024
 *  
 * @brief this
 * @ref A New Model For Control of Systems with Friction - 1995
 * 
*/
template<typename Real, typename Int>
class Lugre
{

public:
    //! @brief default constructor 
    Lugre(Real sigma0, Real sigma1,Real sigma2,
        Real vs, Real coulomb, Real stiction);

    //! @brief intilize the model with null values 
    inline void init();

    //! @brief initlize the model with given values
    //! @param t0 - Intial time value.
    //! @param z0 - Initial deflection state value. 
    //! @param z_dot0 - Initial deflection derivative value.
    inline void init(Real t0, Real z0, Real z_dot0);

    //! @brief  compute the friction torque at time date t
    //! @param  t - time date
    //! @return friction torque computed  
    Real friction_torque(Real velocity, Real t);

    //! @brief  update the avarage deflection state z and it's 
    //! first derivative by one step time.
    inline void  update_state(Real velocity);

    //! @brief compute the stribeck function value.
    //! @return stribeck force 
    inline Real stribek(Real velocity);

    //! @brief compute the torque value at steady state
    inline Real steady_torque(Real velocity);

    //! @brief compute the avarage deflection z at steady state 
    inline Real steady_state(Real velocity);
   
private:
    // Kineatic transition velocity 
    Real vs ;
    //  coulomb force 
    Real coulomb;
    //  Stiction friction force
    Real stiction;
    // Initial and current time
    Real t0, t;
    // Initial and current state 
    Real z0, z;
    // Initial and current state first derivative
    Real z_dot0, z_dot; 
    // Is the model Intilized ?
    bool initialized; 
    // Model static parameters
    Real sigma0, sigma1, sigma2;
};

template<typename Real, typename Int>
class LugreTemp
{

public:
    LugreTemp(Real Te, Real sigma0, Real sigma0e, Real sigma1,
              Real sigma1e,Real sigma2, Real sigma2e,Real vs,
              Real coulomb, Real stiction, Real gradT );
    
    //! @brief intilize the model with null values.
    inline void init();

    //! @brief initlize the time-depend model variables
    //! with given values
    //! @param t0 - Intial time value.
    //! @param z0 - Initial deflection state value. 
    //! @param z_dot0 - Initial deflection derivative value.
    //! @param Ti intial model temperature.
    inline void init(Real Ti, Real t0, Real z0, Real z_dot0);

    //! @brief  compute the friction torque at time date t
    //! @param  t - time date
    //! @return friction torque computed  
    Real friction_torque(Real velocity, Real t);

    //! @brief  update the avarage deflection state z and it's 
    //! first derivative by one step time.
    inline void  update_state(Real velocity);

    //! @brief compute the stribeck function value.
    //! @return stribeck force 
    inline Real stribek(Real velocity);

    //! @brief compute the torque value at steady state
    inline Real steady_torque(Real velocity);

    //! @brief compute the avarage deflection z at steady state 
    inline Real steady_state(Real velocity);

private:
    // initial temperature 
    Real Ti;
    // ambiant temperature.
    Real Te;
    // temperature gradient coefficent
    Real gradT;
    // Kineatic transition velocity 
    Real vs ;
    //  coulomb force 
    Real coulomb;
    //  Stiction friction force
    Real stiction;
    // Initial and current time
    Real t0, t;
    // Initial and current state 
    Real z0, z;
    // Initial and current state first derivative
    Real z_dot0, z_dot; 
    // Is the model Intilized ?
    bool initialized; 
    // Model static parameters
    Real sigma0, sigma0e, sigma1, sigma1e,sigma2,sigma2e;
};


 
} // namespace friction
 
#endif // LUGRE_HPP