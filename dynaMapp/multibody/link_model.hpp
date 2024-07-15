#ifndef LINK_MODEL_HPP
#define LINK_MODEL_HPP

#include <stdexcept>
#include <cmath>
#include "Eigen/Dense"
#include "../utils/constant.hpp"
#include "../multibody/lugre_model.hpp"

using namespace Eigen;

namespace manipulator {
     
/**
 * @brief Genral class template for multibody link.
 * @cite
 * [1] Understanding and Modeling the Behavior of a Harmonic Drive Gear 
 * Transmission Timothy D.Tuttlie - MIT Artificial Intelligence Lab 1993.
 * @cite
 * [2] Modeling, Identification, and Compensation of Friction in Harmonic 
 * Drives, P.S. Gandhi, F.H. Ghorbel, J. Dabney, IEEE 2002. 
*/
template <typename Real, typename Int>
class Link
{
public:
     /**
      * @brief class constructor with minimal given paramters.
      * @param joint_type = 1 for revolute and 0 for prismatric.
      * @param a, d, q, alpha link DH paramters. 
      * @param mx, my, mz  link center of mass coordinates respect to 
      * joint local frame 
      * @param I  link inertia tensor matrix respect to joint local frame.
      * @param q_max, q_min link joint position limits.
      * @param tau_min, tau_max link joint actuator torque limits bounds.
      * @param qp_min, qp_max link joint velocity limits.
      * @param joint_offset link joint offset value at initial time.
      * @param friction_model = "lugre" , "coulomb", "stribeck", "maxwell".
      * @param temperature 
      * @param damping link joint damping coefficent (usully very small value).
      * @param static if "True" the link associated joint is considered fixed.
      * all parameter values will be "reset".
     */
    Link(
        const Real     q,
        const Real     d, 
        const Real     a, 
        const Real     alpha, 
        const Real     mass, 
        const Vector3d r,
        const Matrix3d I, 
        const Int      ratio, 
        const Real     couloumb,
        const Real     stribek, 
        const Real     damping, 
        const bool     static_ = false, 
        const Real     mx,
        const Real     my, 
        const Real     mz,
        const Real     q_max, 
        const Real     q_min, 
        const Real     joint_offset,
        const Real     Jm, 
        const Real     temperature);

    //! @brief initlize model with default values.
    void init();      

    //! @brief this functions is used to reset the rotation matrix 
    //! @brief and the vector p is the link joint is prismatic                         
    inline void transform(const Real& q_);

    //! @return link DH parmter alpha
    Real get_alpha() const { return this->alpha; }  

    //! @return  link DH paramter d
    Real get_d() const { return this->d; } 

    //! @return link joint genralized variable 
    Real get_q() const;

    //! @return minimun bound of joint variable. 
    Real get_q_min() const {return this->q_min;}
    Real get_q_max() const {return this>q_max};

    //! @return joint reduction system ratio
    Real get_ratio() const {return this->ratio; }

    //! @return joint stiffness  
    Real get_stiffness(){ return this->stiffness; }

    //! @return link joint tempertaure 
    //! @note if the variable "temperature_depandancy" is set 
    //! to false then it will return the initial temperature given 
    //! at class constructor
    Real get_temperature(){return this->temperature;}

    //! @return joint offset  
    Real get_joint_offset(){return this->joint_offset;}

    //! @return joint typ "revolute" 0 "prismatic "1"
    Int get_joint_type(){return this->joint_type;}

    //! @return link mass 
    Real get_mass(){return this->mass;}

    //! @return joint global inertia 
    Matrix<Real,3,3> get_inertia(){return this->Jm;}

    //! @return true is the link is considred static.
    bool is_static(){return this->static_;}

    //! @return the joint veclocity at output
    Real get_qm(){this->qm;}

    //! @return the joint velocity at motor shaft
    Real get_qpm(){this->qpm;}

    //! @return 
    Real get_delta_q(){this->delta_q;}

    //! @return time date 
    Real get_time(){this->t};
    
    //! @brief set the link joint generalised variable
    void set_joint_position(Real q);

    //! @brief  DC-Brushed motor constants defalu value = 1.02
    void set_kv(Real kv);

    //! @brief DC-Brushed motor constants
    void set_kt(Real kt);

    //! @brief set manully the joint temperature,
    //! if the "temperature_depandancy" is set to false 
    //! this will not ovverite the temperature value
    void set_temperature(Real temperature_);

    //! @brief set of the sensor stiffness value.
    void set_sensor_stiffness(Real value);

    //! @brief set of harmonic drive reducer stiffness coefficents.
    void set_transmision_stiffness(Real a0, Real a1, Real a2);
    //! @brief set of 
    void set_transmission_stiffness(Real value);

    //! @brief enable temperture depandent computations
    void set_temperature_depandancy(const bool status);

    //! @brief set the link center of mass coordinates vector 
    void set_mc(const Real mx, const Real my, const Real mz);

    //! @brief set the joint global inertia 
    void set_inertia(const Real Jm);

    //! @brief set the link inertia tensor
    void set_inertia_tensor(const Matrix3d I_);

    //! @brief set the link satuts to static 
    void set_static(const bool status);

    //! @brief update time state by one time step 
    void update_time();

    //! @brief set time to t0 
    void reset_time();

    //! @brief update the motor shaft angular spped (RPM) 
    //! given the voltage and current values
    void update_qpm(Real current, Real voltage);

    //! @brief update the angular position at motor shaft
    //! from state k to state k+1 given the qpm at state k 
    void update_qm(Real current, Real voltage);    

    //! @brief compute the motor output torque 
    Real compute_motor_torque(Real current, Real volatge);

    //! @brief torque in joint general contact phenomenas 
    //! at time date t.
    Real contact_friction_torque();    

    //! @brief compute the DC-brushed motor friction torque
    //! @todo add a guard to assert that the friction torque 
    //! do not excced  a ..% value 
    Real motor_friction(Real current, Real volatge); 

    //! @brief estimates the joint torsinal angle
    void update_delta_q(Real current, Real volatge);
                    
    //! @return the link joint overhall stiffness  
    Real compute_stiffness(Real current, Real voltage);

    Matrix3d R ;                ///< Orientation matrix of actual link respect to previous link.
    Vector3d p;                 ///< end of link local position vector
    Vector3d mc;                ///< vector of principal links inertia moments
    Matrix4d T;                 ///< homogenus transformation matrix of the link 
    Real qp,                    ///< Joint output mesured velocity
         qpm,                   ///< joint veolcity at motor shaft 
         qpp,                   ///< joint output mesured accelartion.
         qppm;                  ///< joint accelration at motor shaft.

private:
    Int joint_type;             ///< Joint type.
    Real q,                     ///< q DH parameter.
         qm,                    ///< Motor shaft position
         d,                     ///< d DH parameter.
         a,                     ///< a DH parameter.
         alpha,                 ///< alpha DH parameter.
         q_min,                 ///< Min joint angle.
         q_max,                 ///< Max joint angle.
         delta_q,               ///< the joint  
         joint_offset;          ///< Offset in joint angle (rotoide and prismatic).
    Vector3d r,                 ///< Position of center of link mass respect to global coordinate system.
             p;                 ///< Position vector of actual link respect to previous link.
    Real m,                     ///< Mass of the link.
        Jm,                     ///< Joint Inertia.
        ratio,                  ///< Joint Harmonic drive Gear Ratio. kinova gen 3 = 100
        coulomb,                ///< Coulomb friction coefficient.
        stribek,                ///< Stribek friction coefficent.
        sigma0,                 ///< Lugre friction model paramters
        sigma1,
        sigma2,
        vs,
        temperature,            ///< Joint temperature 
        damping,                ///< Link joint damping 
        mass;                   ///< Link mass
   Matrix3d I;                  ///< Inertia matrix respect center of mass and link coordinate system orientation.
   bool static_;                ///< true if the joint is to be considered locked - ignored for inverse kinematics 
   Real c0, c1, c2;             ///< friction coefficents in the DC-Brushedd motor, consntes given by datasheets
   Real sensor_stiffness,       ///< joint torque sensor stiffness, constant given by datasheets
        a0, a1, a2;             ///< harmonic drive stifness parameters values
   Real kv, kt;                 ///< Joint Motor Constants 
   bool temperature_depandancy; ///< when "true" the calculation of friction torques , and other stuff is 
                                ///< taking into account the evolution of tehieir parmters respect to temperature 
                                ///< defualt value is false !
   Real tau_m_max;              ///< max motor torque output delivred defualt is set to ....from datsheetet
   Real t0 = 0, t;              ///< time date, by default t0 = 0 
};

 

} // namespace manipulator

#include "link_model.cpp"
#endif // LINK_HPP