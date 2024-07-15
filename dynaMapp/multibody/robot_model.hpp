#ifndef ROBOT_MODEL_HPP
#define ROBOT_MODEL_HPP

#include "Eigen/Dense"
#include "link.hpp"

using namespace Eigen;

namespace manipulator {

/**  
 * @brief a template interface for reprsenting general structure  
 * robot manipulators model with 
 * @details  class conati implmentation of  the folowing dynamics 
 *  models : linerized dynamic model
 *  RNE recursive newtom euler alogrithm 
 * 
 * [1] "Elment de robotique" G.clement- Laval University -2021
 * Walker and Orin 
 * 
 */
template<typename Real,typename Int, size_t ndof>
class Robot 
{   
    
public:

    //! @brief default class constructor
    //! initilize each link paramter based on 
    //! given args 
    Robot(const Int      joints_type[ndof], 
        const Real       q[ndof],
        const Real       d, 
        const Real       a, 
        const Real       alpha, 
        const Real       mass[ndof], 
        const Vector3d   r[ndof],
        const Matrix3d   I[ndof], 
        const Int        ratio[ndof], 
        const Real       couloumb[ndof],
        const Real       stribek[ndof], 
        const Real       damping[ndof], 
        const bool       static_[ndof] , 
        const Real       mx[ndof],
        const Real       my[ndof], 
        const Real       mz[ndof],
        const Real       q_max[ndof], 
        const Real       q_min[ndof], 
        const Real       joint_offset[ndof],
        const Real       Jm[ndof], 
        const Real       temperature[ndof],
        const Real       temperature_dependancy[ndof]);

    //! @brief intilize the model time depend paramters
    //! set all initial values to zero.
    inline void init();

    //! @brief return the robot joint "i" variable 
    Real get_joint_position(const Int& i);

    // return the robot joint "i" velocity 
    Real get_joint_velocity(const Int& i);

    //! @brief return the robot joint "i" acceleration
    Real get_joint_acceleration(const Int& i);

    //! @return the robot joints variable vector 
    VectorXd get_joint_position() const;

    //! @brief return the robot joints velocity vector
    VectorXd get_joint_velocity();

    //! @brief return the robot joints acceleration vector 
    VectorXd get_joint_acceleration();

    //! @brief set the robot joint "i" variable value     
    void set_joint_poistion(const Int& i, const Real& qi_);

    //! @brief set the joint "i" velocity value
    void set_joint_velocity(const Int& i, const Real& qpi_);

    //! @brief set the robot joint "i" acceleration value 
    void set_joint_acceleration(const Int& i, const Real& qppi_);

    //! @brief set full joint variable vector
    void set_joint_poistion(const Real(&qpp_)[ndof]);

    //! @brief set full joint velocity vector 
    void set_joint_velocity(const Real(&qpp_)[ndof]);

    //! @brief set full joint  accelartion vector 
    void set_joint_acceleration(const Real (&qpp_)[ndof]);

    //! @brief set the joint "i" static
    //! @note if already static , will skip  
    void set_joint_status(const Int i);

    //! @brief compute the forward kinematic form
    //! until link "j" starting from base link
    //! @return the homogenus transformation matrix 
    //! until link "j".
    Matrix4d forward_kinematic(const int j) const;

    //! @brief compute the forward kinematic chain form 
    //! until the end-effector link 
    Matrix4d forward_kinematic();

    //! @brief compute the robot jacobian matrix 
    MatrixXd jacobian(const Real(&qpp_)[ndof]);
    
    //! @brief compute the inverse kinmatics for a given position
    VectorXd inverse_kinematic();

    //! @brief Joint torque based on Recursive Newton-Euler (RNE) 
    //! formulation. extenal forces and torques are taken into 
    //! account.   
    VectorXd torque(const Real(&q_)[ndof], 
                    const Real(&qp_)[ndof],
                    const Real(&qpp_)[ndof],
                    const Real(&Fext_)[ndof],
                    const Real(&Next_)[ndof]);

    //! @brief Joint torque, without contact force, based on Recursive 
    //! Newton-Euler formulation (RNE).
    VectorXd torque(const Real(&q_)[ndof], 
                    const Real(&qp_)[ndof],
                    const Real(&qpp_)[ndof]);

    //!@brief  Joint torques. when joint velocity is 0, based on Recursive
    //! Newton-Euler formulation.(RNE).
    VectorXd torque_novelocity(const Real(&q_)[ndof]);

   
   
    //! this function use the torque_nonvelocity function 
    MatrixXd inertia_matrix(Real(&q_)[ndof]);

    //! @brief compute the robot inertia matrix, iusing
    //! the ful recursive RNE described in Walker. and Orin.
    //! paper
    MatrixXd inertia_matrix_RNE(Real(&q_)[ndof]);

    //! @brief  compute the coriolis matrix term, regrouping
    //! centrefugal and coriolis effects, using inertia_matrix
    //! coefficents
    MatrixXd coriolis_matrix(Real(&q_)[ndof],Real(&qp_)[ndof]);

    //! @brief Joints acceleration with contact forces. 
    VectorXd acceleration(const Real(&q_)[ndof], 
                        const Real(&qp_)[ndof],
                        const Real(&tau_)[ndof],
                        const Real(&Next_)[ndof]);
    
    //! @return Joints acceleration without contact forces.
    VectorXd acceleration(const Real(&q_)[ndof], 
                        const Real(&qp_)[ndof]);

    //! @brief Joint torque due to gravity based on Recursive 
    //! Newton-Euler formulation (RNE).
    VectorXd torque_gravity();

    //!  @brief oint torque due to centrifugal and Corriolis based on 
    //! Recursive Newton-Euler formulation.
    VectorXd torque_centrifugal(const Real(&qp_)[ndof]);

    //! @brief compute the joints friction torque applied by 
    //! the transmission system and contact pehonema all included 
    //! @todo temperature depandancy ! 
    //! this version  implment no temp depandzncy !  
    VectorXd friction_torque(const Real(&qp_)[ndof]);

    //! @brief Computes the 
    VectorXd delta_torques();

    //! @brief compute the differentifel of the torque 
    //! vector  
    void delta_torque(const Real(&q_)[ndof], 
                    const Real(&qp_)[ndof],
                    const Real(&qpp_)[ndof], 
                    const Real(&dq_)[ndof],
                    const Real(&dqp_)[ndof], 
                    const Real(&dqpp_)[ndof],
                    const Real(&torque_)[ndof],
                    const Real(&dtorque_)[ndof]);

    

    //! @brief compute the the term S(q,qdot,qddot)delta_q
    VectorXd dq_torque();

    //! @brief retun the diagonal matrix of joints stiffness 
    MatrixXd stiffness_matrix(const Real (&q_)[ndof]); 

    //! @brief return the Ng
    MatrixXd ratio_matrix();

    //! @brief compute the state space dynamic model paramters
    //! of the robot
    MatrixXd state_matrix_A();

    //! @brief state matrix C
    MatrixXd state_matrix_C();

    //! @return state matrix B
    MatrixXd state_matrix_B();

    //! @return state matrix D 
    MatrixXd state_matrix_D();   

private:
    // store the robot links data
    Link<Real, Int> links[ndof];
    // global gravity vector
    Vector3d gravity_vector;
    // global vertical vector 
    Vector3d z0;
    // internal computing vectors 
    Vector3d w[ndof],  wp[ndof], 
            vp[ndof],  a[ndof], 
            f[ndof],   f_nv[ndof],
            n[ndof],   n_nv[ndof], dn[ndof],
            F[ndof],   dF[ndof],  
            dp[ndof],  p[ndof], pp[ndof],
            dw[ndof],  dwp[ndof], 
            dvp[ndof], da[ndof], df[ndof],  
            N[ndof],   dN[ndof] ;
    MatrixXd R;
    // is the model intilizied ?
    bool intilizied = false;
};



}; // namespace manipulator

#endif  // ROBOT_MODEL_HPP