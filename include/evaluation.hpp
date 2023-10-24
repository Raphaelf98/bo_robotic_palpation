#ifndef _EVALUATION_HPP_
#define _EVALUATION_HPP_
#include<boost/math/quadrature/trapezoidal.hpp>
#include "stdafx.h"
#include "interpolation.h"
#include "contour.hpp"

class evaluate
{

private: 
    Contour *contour_;
    alglib::spline1dinterpolant* spline_1_,* spline_2_;
    
    /* data */
public:
    evaluate(Contour *contour);
    bool computeIntegrand();
    bool computeArea();
    bool computeFalseNegative();
    bool computeFalsePositive();
    bool computeSpecificity();
    bool computeSensitiviy();
    ~evaluate();
};
#endif