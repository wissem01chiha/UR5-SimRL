#ifndef MULTIBODY_KINEMTICS_HPP
#define MULTIBODY_KINEMTICS_HPP
#include <vector>
#include <cmath>
#include "Eigen/Dense"
using namespace Eigen;

namespace kinematics{

/**
 * @brief orthogonal rotation matrix which transforms a 
 * vector in the (i-1)st coordinate frame to a coordinate 
 * frame which is parallel to the i-th coordinate frame.
 */
template<typename Real>
Matrix<Real,3,3> orthogonal_rotation(const Real theta, const Real alpha);

/**
 * @brief position of the i-th coordinate frame with 
 * respect to the (i— 1)st coordinate frame.
*/
template<typename Real>
Vector<Real,3> frame_position(const Real alpha, const Real d, const Real a);

/**
 * @brief differential form of the orthogonal rotation matrix
 * respect to joint variable.
 */
template<typename Real>
Matrix<Real,3,3> diff_orthogonal_rotation(const Real theta, const Real alpha);

/**
 * @return the partial derivative of the forward kinematics.
 * @tparam N - number of the multibody links 
 */
template<typename Real, size_t N>
Vector<Real,N> kinematic_partial_derivative();

/**
 * @return the homgenus transformation matrix betwwen 
 * link "i" and link "i-1"
 */
template<typename Real>
Matrix<Real,4,4> homogenus_transformation(const Real theta,
                                          const Real alpha,
                                          const Real a,
                                          const Real d);
/**
 * @brief compute the forward kinematic form 
 * until link "j" starting from base link.
 * @return the homogenus transformation matrix 
 * until link "j".
 */
template<typename Real, typename Int>
Matrix<Real,4,4> forward_kinematics(const std::vector<Real> theta,
                                    const std::vector<Real> alpha,
                                    const std::vector<Real> a,
                                    const std::vector<Real> d,
                                    const Int j) ;

/**
 * @return 
 */
MatrixXd jacobien(const int n);



}; // namespace kinematics
#include "multibody_kinemtics.cpp"
#endif // MULTIBODY_KINEMTICS_HPP