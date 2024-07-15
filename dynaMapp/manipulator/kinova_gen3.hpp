#ifndef KINOVA_GEN3_HPP
#define KINOVA_GEN3_HPP
#include "../multibody/robot_model.hpp"

namespace manipulator{  
    
/**
 * @class Kinova 
 * 
 * @brief  ase class for kinova gen3 manipuators
 * model
 * 
 * @ref  
 */
template<typename Real, typename Int, size_t ndof>
class kinova : public Robot<Real,Int,ndof>
{

public:
    kinova(/* args */);
private:
    /* data */
};




}; //namespace manipulator
#endif //  KINOVA_GEN3_HPP


