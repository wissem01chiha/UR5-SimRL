#include "multibody_dynamics.hpp"

namespace dynamics{

template<typename Real, typename Int, size_t ndof>
Torque<Real,Int,ndof>::Torque()
{
    z0 << Real(0), Real(0), Real(1);
    g << Real(0), Real(0), constant::GRAVITY<Real>;
}

template <typename Real, typename Int, size_t ndof>
Vector<Real, ndof> Torque<Real, Int, ndof>::general_torque(const std::vector<Real> q,
                                                           const std::vector<Real> qp, 
                                                           const std::vector<Real> qpp) const
{
    Vector<Real, ndof> ltorque;
    
    if (q.size() != ndof)
    {
        throw std::runtime_error("Error: Dimension mismatch in q. Expected " + 
        std::to_string(ndof) + "elements.");
    }
    if (qp.size() != ndof)
    {
        throw std::runtime_error("Error: Dimension mismatch in qp. Expected " + 
        std::to_string(ndof) + "elements."); 
    }
    if (qpp.size() != ndof)
    {
        throw std::runtime_error("Error: Dimension mismatch in qpp. Expected " + 
        std::to_string(ndof) + "elements."); 
    }
    w[0] = wp[0] << Real(0), Real(0), Real(0);
    vp[0] = -g;
    for (size_t i = 0; i < ndof; i++)
    {
        Matrix<Real,3,3> R = kinematics::orthogonal_rotation(q.at(i), alpha[i]);
        w[i]  = R.transpose()*(w[i-1]+ z0 * qp.at(i));
        wp[i] = R.transpose()*(wp[i-1]+ z0 * qpp.at(i) + w[i-1].cross(z0*qp.at(i)));
        vp[i] = R.transpose()*vp[i-1];
    }
    for (size_t i = ndof-1; i => 0; i--)
    {
        Matrix<Real,3,3> R = kinematics::orthogonal_rotation(q.at(i), alpha[i]);
        vc[i] = vp[i] + wp[i].cross(r[i]) + w[i].cross(w[i].cross(r[i]));
        F[i] = m[i] * vc[i];
        N[i] = Ic[i] * wp[i]+ wp[i].cross(Ic[i] * w[i]);
        f[i] = R *f[i+1] + F[i];
        n[i] = R * n[i+1] + p[i] * f[i] + N[i]+ r[i].cross(F[i]);
        ltorque[i]= n[i].transpose() * ( R.transpose() * z0);
    }
    return ltorque;
}

template <typename Real, typename Int, size_t ndof>
Vector<Real, ndof> Torque<Real, Int, ndof>::torque_no_velocity(const std::vector<Real> q,
                                                               const std::vector<Real> qpp) const
{
   std::vector<Real> qp(ndof,Real(0));
   return general_torque(q,qp,qpp);
}

template <typename Real, typename Int, size_t ndof>
Vector<Real, ndof> Torque<Real, Int, ndof>::gravity_torque(const std::vector<Real> q) const 
{
    
}

template <typename Real,typename Int, size_t ndof>
Matrix<Real, ndof, ndof> Torque<Real,Int,ndof>::inertia(const std::vector<Real> q)
{
    Matrix<Real, ndof, ndof> M;
    Vector<Real,ndof> torque;
   for(size_t i = 1; i <= ndof; i++)
   {
      for(size_t j = 1; j <= ndof; j++)
      {
        torque(j) = (i == j ? Real(1.0): Real(0.0));
      }
      torque = torque_novelocity(q,torque);
      M.col(i) = torque;
   }
    return M;
}

template <typename Real, typename Int, size_t ndof>
Matrix<Real,ndof,ndof> Torque<Real,Int,ndof>::cororlis(const std::vector<Real> q, 
                                                       const std::vector<Real> qp)
{
    return Matrix<Real, ndof, ndof>();
}

}; //namespace dynamics