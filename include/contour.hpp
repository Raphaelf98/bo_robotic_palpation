#ifndef _CONTOUR_HPP_
#define _CONTOUR_HPP_

#include "tumorModel.hpp"
#include "meanShift.hpp"
#include "interpolation.h"
#include "stdafx.h"
#include "param_loader.hpp"
#include "parameters.hpp"
#include<numeric>

#include"helper.hpp"

typedef std::unique_ptr<alglib::spline1dinterpolant> SplineInterpolant1d_ptr;
typedef std::pair<std::shared_ptr<alglib::spline1dinterpolant>,std::shared_ptr<alglib::spline1dinterpolant>> SplineInterpolant_ptr_pair;
typedef std::vector<SplineInterpolant_ptr_pair> SplineInterpolant_ptr_pair_vec;

#define SPLINE_SAMPLES 1000;
/*
This class handles the methods to perform 
    1. Approximation through a gaussian process(GP) of an underlying tumor model
    2. Computation of maximum stiffness regions on posterior of GP in order to find tumor centroids
    3  Computation of points attributed to tumor contour by exploration of Centroids 
    4  Spline approximation of tumor contour
*/
class Contour
{
private:
    
    // Bayesian Optimization variables
    bayesopt::BayesOptBase* bopt_model_;
    std::vector<std::vector<double>> c_;// Holds posterior
    // Mean Shift variables
    MeanShift mean_shift_;
    std::vector<std::vector<double>> ms_data_; 
    double bandwidth_;
    int samples_;
    std::vector<Point> clusters_;
    // Contour Parameters
    size_t n_exploration_directions_;
    size_t c_points_;
    size_t lim_steps_;
    size_t n_directions_;
    double multiplier_;
    //Containers to store approximation data
    std::vector<std::vector<Point>> contours_;
    alglib::spline1dinterpolant spline_1_, spline_2_;
    SplineInterpolant_ptr_pair_vec spline_contours_;
    size_t n_samples_; //number of conotur points to perfrom approximation on
    //k-means variables
    std::vector<double> y_values_;
    size_t total_number_of_iterations_; 
    std::vector<double> tumor_stiffness_vec_;
    double tumor_stiffness_std_;
    double tumor_stiffness_mean_;
    double tumor_stiffness_guess_low_;
    double tumor_stiffness_guess_high_;
    //Variables
    double threshold_;
    double stiffness_threshold_;

    //File handling
    std::string experiment_path_; 
    bool contourPoint_(double &stiffness);
    void labelData_();
    void computeStiffnessThreshold_();
public:
    size_t number_of_step_runs;
    
    //Contour(bayesopt::BayesOptBase* bopt_model, size_t n_exploration_directions);
    Contour(bayesopt::BayesOptBase* bopt_model,ContourParamters cp, std::string experiment_path);
    Contour(){}
    ~Contour();

    //Methods to perform external compuation used by DisplayHeatMap2D class
    size_t getCPoints();
    bayesopt::vectord getLastSample();
    bayesopt::ProbabilityDistribution* getPredictionGaussianProcess(const vectord &q);
    void getInitialSamples( std::vector<double> &samples_x, std::vector<double> &samples_y );
    double evaluateGaussianProcess(const bayesopt::vectord &q);
    double evaluateCriteriaGaussianProcess(const bayesopt::vectord &q);
    void prepareGaussianProcess();
    void stepRunGaussianProcess();
    std::vector<Point> getClusters();//plot clusters
    std::vector<Point> getContourPoints();   //plot contour points
    SplineInterpolant_ptr_pair_vec getSplineInterpolant();

    //function only relevant for evaluation
    std::string getResultsPath();

    //Methods to perform internal compuation
    void runGaussianProcess();
    void computeCluster();
    void exploreContour();
    void approximateContour();
    void computeStiffnessThreshold();

    bool writePosterior();
    void printParameters(const bayesopt::Parameters& par);
};

#endif