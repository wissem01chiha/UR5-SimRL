#ifndef PLOT_HPP
#define PLOT_HPP

#if defined(__WIN32__) || defined(_WIN32)    
    #define GNUPLOT "gnuplot.exe"
    #include <windows.h>
#else   
    #define GNUPLOT "gnuplot"
#endif
#ifndef NCURVESMAX
    #define NCURVESMAX  10 
#endif

#include <vector>
#include <stdio.h>
#include <stdexcept>

/**
 * @class Plot2d
 * @brief  
 */
class Plot2d
{

public:
    //! @brief default class contructor. 
    Plot2d();
    //! @brief member functions.
    void set_title(const std::string & title);
    //! @brief 
    void set_x_label(const std::string & xlabel);
    //! @brief 
    void set_y_label(const std::string & ylabel);
    //! @brief 
    void set_size(int x,int y);
    //! @brief save a gnupolt object to a file.
    bool save(const std::string & filename);
    //!  @brief creates a gnuplot graphic.
    void gnuplot();
    //! @brief show the figure. 
    void show();
private:
   
   std::string  title;         ///< Graph title.
   std::string  x_label;       ///< Graph x axis.
   std::string  y_label;       ///< Graph y axis.
   std::string  gnu_command;   ///< GNU plot command.
 };

Plot2d::Plot2d()
{
}

inline void Plot2d::gnuplot()
{
}

#endif