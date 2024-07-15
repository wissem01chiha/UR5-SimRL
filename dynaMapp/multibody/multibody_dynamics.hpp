#ifndef MULTIBODY_DYNAMICS_HPP
#define MULTIBODY_DYNAMICS_HPP

#include <stddef.h>
#include <stdexcept>
#include <vector>
#include "../utils/constant.hpp"
#include "multibody_kinemtics.hpp"
#include <Eigen/Dense>

using namespace Eigen;

namespace dynamics{

/**
 * @brief Provide multiple implmentattion of computation of 
 * joints torques vector, in differents cases and with differents
 * alogrithms.
 */
template<typename Real, typename Int, size_t ndof>
class Torque
{
public:

    /**
     * @brief class constructor 
     */
    Torque();

    /**
     * @return the torque vetcor applied to the multibody
     * joints using the Recursive Newton-Euler (RNE) formulation.
     * @note  extenal forces and torques are taken into account.
     * @warning the function dosent handle prismatic multibody joint
     * types. i-e by default the boolean joint_type is set to "1". 
     */
    Vector<Real,ndof> general_torque(const std::vector<Real> q,
                                     const std::vector<Real> qp,
                                     const std::vector<Real> qpp) const;
    
    /**
     * @return the torque vector applied to multibody joints,
     * when the joint "j" velocity only is null. 
     */
    Vector<Real,ndof> torque_no_velocity(const Int j) const;

    /**
     * @return the torque vector applied to multibody joints, when
     *  velocities is 0, based on Recursive Newton-Euler formulation.
     */
    Vector<Real,ndof> torque_no_velocity(const std::vector<Real> q,
                                        const std::vector<Real> qpp) const;

    /**
     * @return the torque vector applied to multibody joints, when
     * accelerartions is 0, based on Recursive Newton-Euler formulation.
     */
    Vector<Real,ndof> torque_no_acceleration() const;
    
    /** 
     * @return the torque  vector due to centrifugal and Corriolis effects
     * based on Recursive Newton-Euler formulation.
     */ 
    Vector<Real,ndof> torque_centrifugal();
    
    /**
     * @return the torque vector of the linearized dynamic model formulation
     * [1] 
     */
    Vector<Real,ndof> linear_torque() const;

    /**
     * @return the gravitational force applied to each 
     * multibody joint.
     */
    Vector<Real,ndof> gravity_torque(const std::vector<Real> q) const;
    
    /**
    * @brief we compute the multibody inertia matrix using 
    * the computed torques at null joints velocity.
    * @details  we use this formulation to compute the vector b
    * of gravity and corcolis factor and solve the equation: 
    * 
    *           M(q) q_ddot = tau - b
    * 
    * for the evaluation of the matrix D we can use the formulationn
    * given or use elment-de robotique formulation ! 
    * we use the fact that when joint velocities are null 
    * the acceleration vector qpp=ej
    * 
    * @cite 
    * [1] Efficient Dynamic Computer Simulation of Robotic Mechanisms 
    * M.W.Walker, D.E.Orin ASME-1982
    */
    Matrix<Real,ndof,ndof> inertia(const std::vector<Real> q);

    /**
     * @brief compute and return the corcolis matrix.
     * @details  
     */ 
    Matrix<Real,ndof,ndof> cororlis(const std::vector<Real> q,
                                    const std::vector<Real> qp);

private:

    Real alpha[ndof];
    Real m[ndof];
    Vector<Real,3> z0;
    Vector<Real,3> g;
    Vector<Real,3> w[ndof], wp[ndof];
    Vector<Real,3> v[ndof], vp[ndof];
    Vector<Real,3> vc[ndof];
    Matrix<Real,3,3> Ic[ndof];
    Vector<Real,3> N[ndof],F[ndof],f[ndof],n[ndof];
    Vector<Real,3> r[ndof];
     
}; // class Torque




}; // namespace dynamics
#include "multibody_dynamics.hpp"
#endif // MULTIBODY_DYNAMICS_HPP