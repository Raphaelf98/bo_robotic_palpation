#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_
#include<stdlib.h>
struct TumorModelParameters
{
    //Triangle Parameters
    double triangle_low=1;
    double triangle_high=2;
    double triangle_radius=0.1;
    double triangle_x_trans=0.5;
    double triangle_y_trans=0.5;
    
    double triangle_epsilon=0.1;
    //Rectanlge Parameters
    double rectangle_low=1;
    double rectangle_high=2;
    double rectangle_radius=0.15;
    double rectangle_x_trans=0.5;
    double rectangle_y_trans=0.5;
    double rectangle_epsilon=0.2;
    //Circle Parameters
    double circle_low=1;
    double circle_high=2;
    double circle_radius=0.1;
    double circle_x_trans=0.5;
    double circle_y_trans=0.5;
    double circle_epsilon=0.1;
    //Two Circle Parameters
    double two_circles_low=1;
    double two_circles_high=2;
    double two_circles_radius_1=0.05;
    double two_circles_radius_2=0.1;
    double two_circles_x_trans_1=0.1;
    double two_circles_x_trans_2=0.2;
    double two_circles_y_trans_1=0.7;
    double two_circles_y_trans_2=0.8;
    double two_circles_epsilon=0.1;
};
struct ContourParamters{
    size_t n_exploration_directions=10;
    size_t c_points=100;
    double means_shift_bandwidth=0.05;
    size_t lim_steps=1000;
    double threshold_multiplier=3.0;
};

#endif