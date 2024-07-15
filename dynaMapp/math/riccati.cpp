#include "riccati.hpp"

namespace math{ 
namespace riccati{ 

MatrixXd steady_riccati_equation(
    const Ref<const MatrixXd>& A, 
    const Ref<const MatrixXd>& B, 
    const Ref<const MatrixXd>& Q,
    const Ref<const MatrixXd>& R
){
    return MatrixXd();
}
}; // namespace riccati
}; // namespace math 