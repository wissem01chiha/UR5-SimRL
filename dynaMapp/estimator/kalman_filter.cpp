#include "kalman_filter.hpp"

namespace estimator { 

template<typename Real, typename Int, size_t m, size_t n> 
Kalman<Real, Int, m, n>::Kalman(){
  this->A.setZero(m,n);
  this->Q.setZero(n,n);
  this->C.setZero(m,n);
  this->R.setZero(m,m);
  this->P.setZero(n,n);
}

template<typename Real, typename Int, size_t m, size_t n> 
Kalman<Real,Int,m,n>::Kalman(const MatrixXd &A, const MatrixXd &C,
                             const MatrixXd &Q, const MatrixXd &R,
                             const MatrixXd &P){
  this->A = A;
  this->C = C;
  this->Q = Q;
  this->R = R;
  this->P0 = P;
  this->initialized = false;
  I.setIdentity();
  I.resize(n,n);
  x_hat.resize(n);
  x_hat_new.resize(n);
};

template<typename Real, typename Int, size_t m, size_t n> 
void Kalman<Real,Int,m,n>::init(Real& t0, const VectorXd& x0){
  this->x_hat = x0;
  this->P = P0;
  this->t0 = t0;
  this->t = t0;
  this->initialized = true;
}

template<typename Real, typename Int, size_t m, size_t n> 
inline void  Kalman<Real,Int,m,n>::init() {
  this->x_hat.setZero();
  this->P = P0;
  this->t0 = 0;
  this->t = t0;
  this->initialized = true;
}

template<typename Real, typename Int, size_t m, size_t n> 
void Kalman<Real,Int,m,n>::update(const VectorXd& y) noexcept {
  if(!this->initialized){ 
    throw std::runtime_error("Filter is not initialized!");
  }
  this->x_hat_new = A * this->x_hat;
  this->P=this->A*this->P*this->A.transpose() + this->Q;
  this->K=this->P*this->C.transpose()*(this->C*this->P*this->C.transpose()+this->R).inverse();
  this->x_hat_new += this->K * (y - this->C*this->x_hat_new);
  this->P = (this->I - this->K*this->C)*this->P;
  this->x_hat = this->x_hat_new;
  this->t += constant::TIME_STEP;
}

template<typename Real, typename Int, size_t m, size_t n>
void  Kalman<Real,Int,m,n>::update(const VectorXd& y, const MatrixXd& A){
  this->A = A;
  update(y);
}

template<typename Real, typename Int, size_t m, size_t n>
VectorXd  Kalman<Real,Int,m,n>::get_state(){ return this->x_hat;}

template<typename Real, typename Int, size_t m, size_t n>
Real  Kalman<Real,Int,m,n>::get_time(){ return this->t; };

}; // namespace filter 
 
 
  
 