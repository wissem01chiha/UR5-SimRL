#ifndef KALMAN_HPP
#define KALMAN_HPP

#include <stdexcept>

#include "Eigen/Dense"
#include "utils/constant.hpp"
using namespace Eigen;

namespace estimator {   

/**
 * @brief Kalman template class for disreate kalman filters.
 * [1] An Introduction to the Kalman Filter - G.Welch and G.Bishop
 */
template<typename Real, typename Int, size_t m, size_t n> 
class Kalman {

public:
  //! @brief default constructor  
  //! allocate memory 
  Kalman();

  //! @brief constructor with given parmeters 
  Kalman(const MatrixXd& A,const MatrixXd& C, const MatrixXd& Q,
        const MatrixXd& R, const MatrixXd& P);

  //! @brief Initialize the filter with initial states as zero.
  inline void init();

  //! @brief Initialize the filter with a guess for initial states.
  //! @param  t0 - intial real_ instant value
  //! @param  x0 - initial guess state vector
  inline void init(Real & t0, const VectorXd& x0);

  //! @brief Update the estimated state based on measured values. The
  //! time step is assumed to remain constant.
  //! @param y measured value at time t
  inline  void update(const VectorXd& y);

  //! @brief Update the estimated state based on measured values,
  //! using the given time step and dynamics matrix.
  //! @param  y - Mesured values vector
  //! @param  A - System dynamics matrix
  inline void update(const VectorXd& y, const MatrixXd& A);

  //! @brief  the current state vector.
  VectorXd get_state();

  //! @return  the current time date.
  Real get_time();

private:

  // Matrices for computation
  MatrixXd A, C, Q, R, P, K, P0;
  // Initial and current real_
  Real t0, t;
  // Is the filter initialized?
  bool initialized;
  // n-size identity matrix
  MatrixXd I;
  // Estimated system states
  VectorXd x_hat, x_hat_new;
};

/**
 * @brief A template base class for extended kalman discrete 
 * filters
 */
template<typename Real, typename Int>
class ExtendedKalman
{
  public:
  Extendedkalman( );
 
  private:
    /* data */
};
  

  
};  // namespace filter
 

#endif  