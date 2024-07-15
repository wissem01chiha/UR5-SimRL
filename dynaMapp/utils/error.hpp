#ifndef ERROR_HPP
#define ERROR_HPP

#include <exception>
#include <string>

namespace error { 

/**
 * @class Base class for all DynaMapp exceptions
 * @brief saves an error message that then can be read when catching the 
 * error in the calling code, by calling the what() function.
 * @note for each type of error we crate a custom class that override this
 * class and implment handler specific functions 
 */
class Error : public std::exception
{
public:
    Error(std::string text);
    virtual const char* what() const noexcept;
    virtual void handle() const;
   
protected:
  std::string text_data;
};

 



}; //namespace error

#endif // ERROR_HPP