#include "robot.hpp"

using namespace manipulator;
 
template <typename Real, typename Int, size_t ndof>
Robot<Real, Int, ndof>::Robot(const Int      joints_type[ndof], 
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
                            const Real       temperature_dependancy[ndof])
{ 
    this->gravity_vector << 0, 0, constant::GRAVITY<Real>;
    this->z0 << 0 , 0 , Real(1);
    for (Int i = 0; i < ndof; i++){

        this->links[i].set_q(q[i]);
        this->links[i].set_d(d[i]);
        this->links[i].set_alpha(alpha[i]);
        this->links[i].set_mass(mass[i]);
        this->links[i].set_inertia_tensor(I[i]);
        this->links[i].set_ratio(ratio[i]);
        this->links[i].set_coulomb(couloumb[i]);
        this->links[i].set_stribek(stribek[i]);
        this->links[i].set_mc(mx[i], my[i], mz[i]);
        this->links[i].set_joint_offset(joint_offset[i]);
        this->links[i].set_inertia(Jm[i]);
        this->links[i].set_temperature_depandancy(temperature_depandancy[i]);
        this->links[i].set_temperature(temperature[i]);

    }
    
    this->intilizied = true ;  
}

template <typename Real, typename Int, size_t ndof>
inline void Robot<Real, Int, ndof>::init(){}

template <typename Real, typename Int, size_t ndof>
Real Robot<Real, Int, ndof>::get_q(const Int &i){
    return this->links[i].get_q(i);
}

template <typename Real, typename Int, size_t ndof>
Real Robot<Real, Int, ndof>::get_qp(const Int &i){
    return this->links[i].get_qp(i);
}

template <typename Real, typename Int, size_t ndof>
Real Robot<Real, Int, ndof>::get_qpp(const Int &i){
    return this->links[i].get_qpp(i);
}

template <typename Real, typename Int, size_t ndof>
VectorXd Robot<Real, Int, ndof>::get_q() const {
    VectorXd result = VectorXd::Zero(ndof);
    for (size_t i = 0; i < ndof; i++){
        result[i]=this->links[i].get_q();
    }
    return result;
}

template <typename Real, typename Int, size_t ndof>
VectorXd Robot<Real, Int, ndof>::get_qp(){
    VectorXd result = VectorXd::Zero(ndof);
    for (size_t i = 0; i < ndof; i++){
        result[i]=this->links[i].get_qp();
    }
    return result;
}

template <typename Real, typename Int, size_t ndof>
VectorXd Robot<Real, Int, ndof>::get_qpp(){
    VectorXd result = VectorXd::Zero(ndof);
    for (size_t i = 0; i < ndof; i++){
        result[i]=this->links[i].get_qpp();
    }
    return result;
}

template <typename Real, typename Int, size_t ndof>
void  Robot<Real, Int, ndof>::set_q(const Int &i, const Real &qi_){
    this->links[i].set_q(qi_); }

template <typename Real, typename Int, size_t ndof>
void  Robot<Real, Int, ndof>::set_qp(const Int &i, const Real &qpi_){
    this->links[i].set_qp(qpi_); }

template <typename Real, typename Int, size_t ndof>
void Robot<Real, Int, ndof>::set_qpp(const Int &i, const Real &qppi_){
    this->links[i].set_qpp(qppi_); }

template <typename Real, typename Int, size_t ndof>
void  Robot<Real, Int, ndof>::set_q(const Real (&q_)[ndof]){
    for (Int i = 0; i < ndof; i++){
        set_q(i,q_[i]);
    } 
}

template <typename Real, typename Int, size_t ndof>
void Robot<Real, Int, ndof>::set_qp(const Real (&qp_)[ndof]){
    for (Int i = 0; i < ndof; i++){
        set_qp(i, qp_[i]);
    }
}

template <typename Real, typename Int, size_t ndof>
void  Robot<Real, Int, ndof>::set_qpp(const Real (&qpp_)[ndof]){
    for (Int i = 0; i < ndof; i++){
        set_qpp(i, qpp_[i]);
    }
}

template <typename Real, typename Int, size_t ndof>
void Robot<Real, Int, ndof>::set_static(const Int i){
     if (!this->links[i].is_static()){
        this->links[i].set_static(true);
     }
}

template <typename Real, typename Int, size_t ndof>
Matrix4d Robot<Real, Int, ndof>::forward_kinematic(const int j) const{
    Matrix4d Thomo = Matrix4d::Identity(4,4) ;
    for (Int i = 0; i < j; i++){
            Thomo = Thomo *  this->links[i].T; 
    }
    return Thomo;  
}

template <typename Real, typename Int, size_t ndof>
Matrix4d Robot<Real, Int, ndof>::forward_kinematic(){
    return forward_kinematic(ndof);
}

template <typename Real, typename Int, size_t ndof>
MatrixXd Robot<Real, Int, ndof>::jacobian(const Real (&q_)[ndof]){
    MatrixXd J = MatrixXd::Zero(6,ndof);
    Vector3d r[ndof], z[ndof];
    for (Int i = 0; i < ndof; i++){
        if (i == 0){
            z[i] << 0, 0, 1; 
            }else{
                z[i]= this->links[i].R*z[i-1];
            }
    }
    for (Int i = ndof; i>=0; i--){
        if (i == ndof){
            r[i]= this->links[i].p;
        }else{
            r[i]=this->links[i].p+ this->links[i].R*r[i+1];
        }
    }
    for (Int i = 0; i < ndof; i++){
        if (this->links[i].get_joint_type()==0){
            J.col(i) << z[i].cross(r[i]),z[i];
        }else{
            J.col(i) << z[i],0;
        }
    }
    return J; 
}

template <typename Real, typename Int, size_t ndof>
VectorXd Robot<Real, Int, ndof>::inverse_kinematic()
{
    return VectorXd();
}

template <typename Real, typename Int, size_t ndof>
VectorXd Robot<Real, Int, ndof>::torque(const Real (&q_)[ndof], 
                                        const Real (&qp_)[ndof], 
                                        const Real (&qpp_)[ndof], 
                                        const Real (&Fext_)[ndof], 
                                        const Real (&Next_)[ndof])
{
   VectorXd ltorque = VectorXd::Zero(ndof,1); 
   Matrix3d Rt, temp;
   set_q(q_);
   set_qp(qp_);
   set_qpp(qpp_);
   this->vp[0] << this->gravity_vector;
   for (Int i = 0; i < ndof; i++)
   {
        Rt = this->links[i].R;
        if (this->links[i].get_joint_type()==0)
        {
            this->w[i] = Rt*(w[i-1] + z0*qp(i));
            this->wp[i] = Rt*(wp[i-1] + z0*qpp(i)
                    + this->w[i-1].cross(this->z0 * get_qp(i)));
            this->vp[i] = this->wp[i].cross(this->link[i].p)
                 + this->w[i].cross(this->w[i].cross(this->links[i].p))
                 + Rt * (this->vp[i-1]);
        }else{
            this->w[i] = Rt * this->w[i-1];
            this->wp[i] = Rt*wp[i-1];
         vp[i] = Rt*(vp[i-1] + z0*qpp(i))
                 + 2.0*CrossProduct(w[i],Rt*z0*qp(i))
                 + CrossProduct(wp[i],p[i])
                 + CrossProduct(w[i],CrossProduct(w[i],p[i]));
        }
        

    
   }
   

}

template <typename Real, typename Int, size_t ndof>
VectorXd Robot<Real, Int, ndof>::torque(const Real (&q_)[ndof],
                                        const Real (&qp_)[ndof],
                                        const Real (&qpp_)[ndof])
{
     
}

template <typename Real, typename Int, size_t ndof>
VectorXd Robot<Real, Int, ndof>::torque_novelocity(const Real (&qpp_)[ndof])
{
    VectorXd ltorque = VectorXd::Zero(ndof,1);
    MatrixXd Rt, temp;
    VectorXd vp = VectorXd::Zero(ndof);
    for (Int i = 0; i < ndof; i++){
        Rt = this->links[i].R;
        if (this->links[i].joint_type == 0){
            this->wp[i] = Rt *(this->wp[i-1] + this->z0*get_qpp(i));
            this->vp[i] = this->wp[i].cross(this->p[i])+ Rt * this->vp[i];
        }else{
            this->wp[i] = Rt * this->wp[i-1];
            this->vp[i] = Rt * (this->vp[i-1] + this->z0*get_qpp(i))
                 + wp[i].cross(this->p[i]);
        }
        this->a[i] = this->wp[i].cross(this->links[i].p) + this->vp[i];  
    }
    for (Int i = ndof-1; i => 0; i--){
        this->F[i] = this->a[i] * this->links[i].get_mass();
        this->N[i] = this->links[i].I * wp[i];
        if (i == dof){
            this->f_nv[i]=F[i];
            this->n_nv[i]=this->p[i].cross(this->f_nv[i])+this->links[i].p.cross(this->F[i]) 
            +this->N[i]; 
        }else{
            this->f_nv[i]=this->links[i+1].R*this->f_nv[i+1]+this->F[i];
            this->n_nv[i]=this->links[i+1].R*this->n_nv[i+1]+this->p[i].cross(this->f_nv[i])
            + this->links[i].r.cross(this->F[i]) + this->N[i];
        }
        if(links[i].get_joint_type() == 0){ 
            temp = ((this->z0 * this->links[i].R) * this->n_nv[i]);
        }
        else{ 
            temp = ((this->z0 * this->links[i].R) * this->f_nv[i]);
        }
        ltorque(i) = temp(1,1)+
        this->links[i].Jm*this->links[i].get_ratio()*links[i].get_ratio()*get_qpp(i);
    }
    return ltorque;
}

template <typename Real, typename Int, size_t ndof>
inline MatrixXd Robot<Real, Int, ndof>::inertia_matrix(Real (&q_)[ndof]){
    MatrixXd M = MatrixXd::Zero(ndof,ndof);
    Real e[ndof] ;
    set_q(q_);
    for (Int i = 0; i < ndof; i++){
        for (Int j = 0; j < ndof; j++){
            e(j)=(i==j ? : Real(1.0) : Real(0.0));
        }
        M.col(i) << torque_novelocity(e);
    }
    return M;
}

template <typename Real, typename Int, size_t ndof>
MatrixXd Robot<Real, Int, ndof>::inertia_matrix_RNE(Real (&q_)[ndof]){
    return MatrixXd();
}

template <typename Real, typename Int, size_t ndof>
MatrixXd Robot<Real, Int, ndof>::coriolis_matrix(Real (&q_)[ndof], 
                                                Real (&qp_)[ndof]){
    
    
}

template <typename Real, typename Int, size_t ndof>
VectorXd Robot<Real, Int, ndof>::acceleration(const Real (&q_)[ndof], 
                                              const Real (&qp_)[ndof], 
                                              const Real (&tau_)[ndof], 
                                              const Real (&Next_)[ndof])
{
    return VectorXd();
}

template <typename Real, typename Int, size_t ndof>
VectorXd Robot<Real, Int, ndof>::acceleration(const Real (&q_)[ndof], 
                                            const Real (&qp_)[ndof] ){
    return VectorXd();
}

template <typename Real, typename Int, size_t ndof>
VectorXd Robot<Real, Int, ndof>::torque_gravity(){
    VectorXd ltorque = VectorXd::Zero(ndof);
    MatrixXd Rt, temp;
    this->vp[0] = (Real)constant::GRAVITY;
    for (Int i = 0; i < ndof; i++){
        Rt = this->links[i].R;
        if(this->links[i].get_joint_type() == 0){ 
            this->vp[i] = Rt*(this->vp[i-1]);
        }
        else{ 
            this->vp[i] = Rt * this->vp[i-1];
            this->a[i] = vp[i];
        }
    }
    for(Int i = ndof-1; i => 0; i--) {
      this->F[i] =  this->a[i] * this->links[i].get_mass();
      if(i == ndof) {
         this->f[i] = this->F[i];
         this->n[i] = this->p[i].cross(this->f[i])
                + this->links[i].p.cross(this->F[i]);
      } else {
         this->f[i] = this->links[i+1].R * this->f[i+1] + this->F[i];
         this->n[i] = this->links[i+1].R * this->n[i+1] + this->p[i].cross(this->f[i])
                +this->links[i].p.cross(this->F[i]);
      }
      if(this->links[i].get_joint_type() == 0){ 
            temp = ((this->z0.t()*links[i].R)*n[i]);
      }
      else{ 
         temp = ((this->z0 * this->links[i].R) * this->f[i]);}
      ltorque(i) = temp(1,1);
   }
    return ltorque;
}

template <typename Real, typename Int, size_t ndof>
MatrixXd Robot<Real, Int, ndof>::stiffness_matrix(const Real (&q_)[ndof])
{
     
}

template <typename Real, typename Int, size_t ndof>
MatrixXd Robot<Real, Int, ndof>::ratio_matrix()
{
    MatrixXd result = MatrixXd::Zero(ndof,ndof);
    for (Int i = 0; i < ndof; i++){
        result(i,i)= this->links[i].get_ratio(); 
    }
    return result;
}

template <typename Real, typename Int, size_t ndof>
MatrixXd Robot<Real, Int, ndof>::state_matrix_A()
{
    return MatrixXd();
}

template <typename Real, typename Int, size_t ndof>
MatrixXd  Robot<Real, Int, ndof>::state_matrix_C()
{
    return MatrixXd();
}

template <typename Real, typename Int, size_t ndof>
Eigen::MatrixXd Robot<Real, Int, ndof>::state_matrix_B()
{
    return MatrixXd();
}

template <typename Real, typename Int, size_t ndof>
Eigen::MatrixXd Robot<Real, Int, ndof>::state_matrix_D()
{
    return Eigen::MatrixXd();
}
