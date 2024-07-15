#include "link.hpp"

namespace manipulator{ 
 
template <typename Real, typename Int>
inline Link<Real, Int>::Link(const Int       joint_type, 
                            const Real       q,
                            const Real       d, 
                            const Real       a, 
                            const Real       alpha, 
                            const Real       mass, 
                            const Vector3d   r,
                            const Matrix3d   I, 
                            const Int        ratio, 
                            const Real       couloumb,
                            const Real       stribek, 
                            const Real       damping, 
                            const bool       static_ = false, 
                            const Real       mx,
                            const Real       my, 
                            const Real       mz,
                            const Real       q_max, 
                            const Real       q_min, 
                            const Real       joint_offset,
                            const Real       Jm, 
                            const Real       temperature)
{
    this->joint_type = joint_type;
    this->kv = 125;  
    this->kt = 0.076;
    this->a0 = 0.02;
    this->a1 = 0.02;
    this->a2 = 0.02;
    this->c0 = 0.118;
    this->c1 = 0.0187;
    this->c2 = 0.0101;
    this->qp = 0.0 ;
    this->qpm = 0.0 ;
    this->qpp = 0.0;
    this->qppm = 0.0;
    this->sensor_stiffness = 266;
    this->vs = 0.14;

    this->damping = damping;
    this->alpha   = alpha;
    this->a       = a;
    this->d       = d;
    this->ratio   = ratio;
    this->coulomb = couloumb;
    this->static_ = static_;
    this->stribek = stribek;
    this->mass    = mass; 
    this->q_max = q_max;
    this->q_min = q_min;
    this->joint_offset =joint_offset;
    this->r=r;
    this->I=I;
    this->sigma0 = 0.015;
    this->sigma1 = 0.15;
    this->sigma2 = 0.12;

    this->t= 0;
    this->R.setZero(3,3);
    this->p.setZero(3,1);
    this->T.setZero(4,4);

    if (joint_type == 0){
	    this->q += joint_offset;
    }
    else{
	    this->d += joint_offset;
    }
    Real ct, st, ca, sa;
    ct = cos(this->q);
    st = sin(this->q);
    ca = cos(this->alpha);
    sa = sin(this->alpha);
    this->R(1,1) = ct;
    this->R(2,1) = st;
    this->R(3,1) = 0.0;
    this->R(1,2) = -ca*st;
    this->R(2,2) = ca*ct;
    this->R(3,2) = sa;
    this->R(1,3) = sa*st;
    this->R(2,3) = -sa*ct;
    this->R(3,3) = ca;

    this->p(1) = a*ct;
    this->p(2) = a*st;
    this->p(3) = d;
    
    this->mc(1) = mx;
    this->mc(2) = my;
    this->mc(3) = mz;
    this->temperature_depandancy = false ;
    this->tau_m_max = 17.00;
    this->T<< R, p 
              0, 1; 
}

template <typename Real, typename Int>
Link<Real, Int>::~Link(){}

template <typename Real, typename Int>
inline void Link<Real, Int>::transform(const Real& q_){
    if(this->joint_type == 0){ // revolute joint  
        Real ct, st, ca, sa;
        this->q = this->q + this->joint_offset;
        ct = cos(this->q);
        st = sin(this->q);
        ca = this->R(3,3);
        sa = this->R(3,2);

        this->R(1,1) = ct;
        this->R(2,1) = st;
        this->R(1,2) = -ca*st;
        this->R(2,2) = ca*ct;
        this->R(1,3) = sa*st;
        this->R(2,3) = -sa*ct;
        this->p(1) = this->a*ct;
        this->p(2) = this->a*st;
        }
      else // prismatic joint
        this->p(3) = this->q + this->joint_offset;
        this->d = this->q + this->joint_offset;
}

template <typename Real, typename Int>
Real Link<Real, Int>::get_q() const {
    if (this->joint_type == 0){
        return (this->q - this->joint_offset);
    }else{
        return (this->d - this->joint_offset);
    } 
}

template <typename Real, typename Int>
void  Link<Real, Int>::set_q(Real q_){
    if (this->joint_type == 0){
        this->q = q_ + this->joint_offset;
    }else{
        this->d=q_ + this->joint_offset;
    }
}

template <typename Real, typename Int>
void Link<Real, Int>::set_kv(Real kv_) { this->kv = kv_; }

template <typename Real, typename Int>
void Link<Real, Int>::set_kt(Real kt_){this->kt = kt_;}

template <typename Real, typename Int>
void  Link<Real, Int>::set_temperature(Real temperature_){
    this->temperature = temperature_;}

template <typename Real, typename Int>
void Link<Real, Int>::set_sensor_stiffness(Real value){
    this->sensor_stiffness = value;}

template <typename Real, typename Int>
void Link<Real, Int>::set_transmision_stiffness(Real a0, Real a1, Real a2){
    this->a0 = a0; this->a1 = a1;  this->a2 = a2;
}

template <typename Real, typename Int>
void  Link<Real, Int>::set_temperature_depandancy(const bool status){
    this->temperature_depandancy = status;
}

template <typename Real, typename Int>
void Link<Real, Int>::set_mc(const Real mx, const Real my, const Real mz){
    this->mc(0) = mx;
    this->mc(1) = my;
    this->mc(2) = mz;
}

template <typename Real, typename Int>
void Link<Real, Int>::set_inertia(const Real Jm_){
    this->Jm=Jm_;
}

template <typename Real, typename Int>
void Link<Real, Int>::set_inertia_tensor(const Matrix3d I_){
    this->I=I_;
}

template <typename Real, typename Int>
void Link<Real, Int>::set_static(const bool status){
    this->static_ =status;
}

template <typename Real, typename Int>
void Link<Real, Int>::update_time(){this->t += constant::TIME_STEP; }

template <typename Real, typename Int>
void Link<Real, Int>::reset_time(){ this->t =0}

template <typename Real, typename Int>
void  Link<Real, Int>::update_qpm(Real current, Real voltage){ 
    update_time<Real,Int>();
    this->qpm = this->kv * voltage; 
}

template <typename Real, typename Int>
void Link<Real, Int>::update_qm(Real current, Real voltage){   
    update_qpm<Real,Int>(current,voltage);
    this->qm += this->qpm * constant::TIME_STEP;
}

template <typename Real, typename Int>
Real Link<Real, Int>::motor_friction(Real current, Real voltage){
    update_qm(current, voltage);
    return this->c0+this->c1*this->qm*this->qm+this->c2*this->qm*this->qm*this->qm;
}

template <typename Real, typename Int>
Real  Link<Real, Int>::compute_motor_torque(Real current, Real voltage) noexcept {   
    try{
        Real torque =  this->kt * current;
         if (torque > this->tau_m_max) {
            throw std::runtime_error("Computed torque exceeds maximum allowed value");
        }
        return torque;
    }
    catch(const std::runtime_error& e){
        std::cerr << "Warning: " << e.what() << std::endl;
        return this->tau_m_max;
    }
}

template <typename Real, typename Int>
Real Link<Real, Int>::contact_friction() noexcept {   
    if (!this->temperature_depandancy)
    {   
        friction::LuGre<Real,Int> lugre(this->sigma0,this->sigma1,this->sigma2,
        this->vs,this->coulomb, this->stribek);
        lugre.init<Real, Int>();
        return lugre.friction_torque<Real, Int>(this->qp,this->t);
    }else{}
}

template <typename Real, typename Int>
void Link<Real, Int>::update_delta_q(Real current , Real voltage) noexcept {   
    try{
        update_qm<Real,Int>(current,voltage);
        if (this->ratio <= 0){
             throw std::runtime_error("Non Positive Transmission Ratio Value"); 
        }
        this->delta_q = this->qm-this->q * 1/this->ratio;
    }
    catch(const std::exception& e){
        std::cerr << "Warning: " << e.what() << std::endl;
    }
}

template <typename Real, typename Int>
Real Link<Real, Int>::compute_stiffness(Real current, Real voltage){
    update_delta_q<Real,Int>(current,voltage);
    return a0+this->delta_q*a0+a1*this->delta_q+a2*this->delta_q*this->delta_q;
}
 
}; // namespace manipulator 