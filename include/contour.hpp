#ifndef _CONTOUR_HPP_
#define _CONTOUR_HPP_

#include "tumorModel.hpp"
#include "meanShift.hpp"
#include "interpolation.h"
#include "stdafx.h"
inline std::vector<double> linspace(double min,double max,int n){
    std::vector<double> a;
    if(n<1){n=1;}
    a.resize(n);
    for(int i=0;i<n;++i){a[i]=min+(max-min)*i/(n-1);}
    return a;
};

class Contour
{
private:
    // Bayesian Optimization variables
    bayesopt::BayesOptBase* bopt_model_;
    size_t n_exploration_directions_;
    size_t c_points_;
    std::vector<std::vector<double>> c_;
    // Mean Shift variables
    MeanShift mean_shift_;
    std::vector<std::vector<double>> ms_data_; 
    double bandwidth_;
    int samples_;
    std::vector<Point> clusters_;
    // Contour Parameters
    size_t lim_steps_;
    double stiffness_threshold_;
    size_t n_directions_;
    std::vector<std::vector<Point>> contours_;
    // Approximation Parameters
    alglib::spline1dinterpolant spline_1_, spline_2_;
    size_t n_samples_;
public:
    Contour(bayesopt::BayesOptBase* bopt_model, size_t n_exploration_directions);
    ~Contour();
    void runGaussianProcess();
    void computeCluster();
    void exploreContour();
    void approximateTContour();
    std::vector<alglib::spline1dinterpolant*> getSplineInterpolant();
};

#endif