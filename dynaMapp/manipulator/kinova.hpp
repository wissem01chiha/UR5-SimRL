#ifndef KINOVA_HPP
#define KINOVA_HPP
#include "robot.hpp"

namespace manipulator{  
    
/**
 * @class Kinova 
 * @copyright 
 * @date 
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
#endif


