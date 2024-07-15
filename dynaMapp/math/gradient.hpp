#ifndef GRADIENT_HPP
#define GRADIENT_HPP

#include "Eigen/Dense"
using namespace Eigen;

namespace math { 

namespace gradient {

/**
 * @return the gradient of scalar function " f "
 * with respect to some or all variables "x".
 * @note we assume that the function f satisfy the 
 * conditions that the derivative exsits 
*/
template<typename fun, typename... Args, typename... Vars>
void gradient(const fun& f, const std::tuple<Args...>& args,const std::tuple<Vars...>& vars);



}; // namespace gradient

}; // namespace math















#endif 