#ifndef RICCATI_SOLVER_HPP
#define RICCATI_SOLVER_HPP

#include "Eigen/Dense"

using namespace Eigen;
 

/**
 * @brief Computes the unique solution X of the discrete
 * algebric riccati equation :
 * 
 *      AᵀXA − X − AᵀXB(BᵀXB + R)⁻¹BᵀXA + Q = 0 
 * 
 * @throws std::exception if Q is not symmetric positive semidefinite.
 * @throws std::exception if R is not symmetric positive definite.
 * @throws std::exception if (A, B) isn't a stabilizable pair.
 * @throws std::exception if (A, C) isn't a detectable pair where Q = CᵀC.
 */
MatrixXd steady_riccati_equation(
    const Ref<const MatrixXd>& A, 
    const Ref<const MatrixXd>& B, 
    const Ref<const MatrixXd>& Q,
    const Ref<const MatrixXd>& R);

/**
 * @brief Computes numericlly a solution of the SDRE Problem 
 * given by :
 *  
 *      Ⅱ(x)A(x)+A(x)ᵀⅡ(x)−Ⅱ(x)BR⁻¹BᵀⅡ(x)+ Q = 0
 * 
 * we split A(x)= A₀ + ε ∆A(x) with a power serie expension 
 * for of  Ⅱ(x) matrix. power series terms are computed 
 * recursivlly.
 * @note the term ∆A(x), is computed by an other function
 * not yet implemented.
 * @return the matrix Ⅱ(x) evaluted at state x. 
 */
MatrixXd state_depand_riccati_equation(const Ref<const VectorXd>& state);

/*
[1] E. K.-W. Chu, H.-Y. Fan, W.-W. Lin & C.-S. Wang "Structure-Preserving
    Algorithms for Periodic Discrete-Time Algebraic Riccati Equations",
    International Journal of Control, 77:8, 767-788, 2004.
    DOI: 10.1080/00207170410001714988

*/



#endif // RICCATI_SOLVER_HPP