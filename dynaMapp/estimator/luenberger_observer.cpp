#include "luenberger_observer.hpp"

namespace estimator{ 

template<typename Real, typename Int, size_t n, size_t m>
Luenberger<Real, Int, n, m>::Luenberger(MatrixXd A, 
                                        MatrixXd  B, 
                                        MatrixXd  C, 
                                        MatrixXd  L){

}

template<typename Real, typename Int, size_t n, size_t m>
bool  Luenberger<Real, Int, n, m>::is_observable()
{
    return false;
}

template<typename Real, typename Int, size_t n, size_t m>
MatrixXd  Luenberger<Real, Int, n, m>::state()
{
return MatrixXd();
}
template<typename Real, typename Int, size_t n, size_t m>
MatrixXd  Luenberger<Real, Int, n, m>::error()
{
return MatrixXd();
}

template<typename Real, typename Int, size_t n, size_t m>
bool  Luenberger<Real, Int, n, m>::is_stable()
{
return false;
}

}; // namespace estimator
