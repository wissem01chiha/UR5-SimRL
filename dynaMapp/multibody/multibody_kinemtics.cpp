#include "multibody_kinemtics.hpp"

namespace kinematics{

template <typename Real>
Matrix<Real,3,3> orthogonal_rotation(const Real theta, const Real alpha)
{
    Matrix<Real,3,3> R;
    R << cos(theta), -sin(alpha)*cos(theta), sin(alpha)*cos(theta),
         sin(theta), cos(alpha)*sin(theta),  -sin(alpha)*cos(theta),
         Real(0)           sin(alpha),            cos(alpha);
         return R;
}

template <typename Real>
Vector<Real, 3> frame_position(const Real alpha, const Real d, const Real a)
{
     Vector<Real, 3> p;
     p << a, d* sin(alpha), d*cos(alpha);
     return p;
}

template <typename Real>
Matrix<Real,3,3> diff_orthogonal_rotation(const Real theta, const Real alpha)
{
    Matrix<Real,3,3> Q;
    Q << Real(0), Real(-1), Real(0),
         Real(1), Real(0),  Real(0),
         Real(0), Real(0),  Real(0);
    return Q * orthogonal_rotation(theta, alpha);
}

template <typename Real, size_t N>
Vector<Real, N> kinematic_partial_derivative()
{
    return Vector<Real, N>();
}

template <typename Real>
Matrix<Real,4,4> homogenus_transformation(const Real theta, 
                                          const Real alpha, 
                                          const Real a, 
                                          const Real d)
{
     Matrix<Real, 3, 3> R = orthogonal_rotation(theta,alpha) ;
     Vector<Real, 3> p = frame_position(alpha,d,a);
     Matrix<Real, 4, 4> T;
     T << R, p,
         Real(0), Real(1);
    return T;
}

template <typename Real, typename Int>
Matrix<Real, 4, 4> forward_kinematics(const std::vector<Real> theta, 
                                      const std::vector<Real> alpha, 
                                      const std::vector<Real> a, 
                                      const std::vector<Real> d, 
                                      const Int j)
{
    Matrix<Real,4,4> Thomo;
    Thomo.setIdentity();
    for (size_t i = 0; i <= j; i++)
    {
        Thomo = Thomo * homogenus_transformation(theta,alpha,a,d); 
    }
    retrun Thomo;
}

}; //namespace kinematics