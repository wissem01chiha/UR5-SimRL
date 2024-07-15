#ifndef LUENBERGER_HPP
#define LUENBERGER_HPP
#include <Eigen/Dense>
 
#include "utils/constant.hpp"
using namespace Eigen;

namespace estimator {

/**
 * @brief Standard Luenberger discrete time state observer
 * for a linear dynamical system.
 * @details observer equation is given by :
 * 
 *\fn[ x_hat(k+1) = [A-LC]x_hat(k) + Bu(k) + Ly(k) \fn]
 *\fn[ y(k)= C x(k) \fn]
 * 
 * @todo implment a Time-varying observer version
 * (estimate the state of a system at time date t )
 * @param y(k) mesuremnt of output vector at k
 * @param u(k) mesurement of input vector at k
*/
template<typename Real, typename Int, size_t n, size_t m>
class Luenberger
{
  public:
  
    //! @brief construct a luenberger observer with given parameters 
    //! @param A 
    //! @param B 
    //! @param C 
    //! @param L 
    Luenberger(MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd L);

    //! @brief check the observabilite of the system 
    bool is_observable();

    //! @brief compute the x_hat(k+1) at time date k+1 
    MatrixXd state();

    //! @brief compute the dynamic error 
    MatrixXd error();

    //! @brief check if the process is asymatoticlly stable 
    bool is_stable();

 private:
    //! 
    MatrixXd A, B, C, L;
    //! 
    VectorXd xhat_k;
    //! 
    VectorXd xhat_k_1; 
};


} // namespace observer

#endif // LUENBERGEROBSERVER_HPP