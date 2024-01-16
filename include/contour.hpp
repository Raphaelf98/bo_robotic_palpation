#ifndef _CONTOUR_HPP_
#define _CONTOUR_HPP_

#include "tumorModel.hpp"
#include "meanShift.hpp"
#include "interpolation.h"
#include "stdafx.h"
//#include "fileParser.hpp"
#include "param_loader.hpp"
#include "parameters.hpp"
#include<numeric>

#include"helper.hpp"

typedef std::unique_ptr<alglib::spline1dinterpolant> SplineInterpolant1d_ptr;
typedef std::pair<std::shared_ptr<alglib::spline1dinterpolant>,std::shared_ptr<alglib::spline1dinterpolant>> SplineInterpolant_ptr_pair;
typedef std::vector<SplineInterpolant_ptr_pair> SplineInterpolant_ptr_pair_vec;

class Contour
{
private:
    //File handling
    std::string experiment_path_; 
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
    SplineInterpolant_ptr_pair_vec spline_contours_;
    size_t n_samples_;
    //k-means variables
    std::vector<double> y_values_;
    size_t total_number_of_iterations_;
    std::vector<double> tumor_stiffness_vec_;
    double tumor_stiffness_std_;
    double tumor_stiffness_mean_;
    //PARAMETER
    double multiplier_;
    double threshold_;
    bool contourPoint_(double &stiffness);
    void labelData_();
    void computeStiffnessThreshold_();
public:
    size_t number_of_step_runs;
    
    Contour(bayesopt::BayesOptBase* bopt_model, size_t n_exploration_directions);
    Contour(bayesopt::BayesOptBase* bopt_model,ContourParamters cp, std::string experiment_path);
    Contour(){}
    ~Contour();

    size_t getCPoints();
    bayesopt::vectord getLastSample();
    bayesopt::ProbabilityDistribution* getPredictionGaussianProcess(const vectord &q);
    void getInitialSamples( std::vector<double> &samples_x, std::vector<double> &samples_y );
    double evaluateGaussianProcess(const bayesopt::vectord &q);
    double evaluateCriteriaGaussianProcess(const bayesopt::vectord &q);
    void prepareGaussianProcess();
    void stepRunGaussianProcess();

    void runGaussianProcess();
    void computeCluster();
    void exploreContour();
    void approximateContour();
    void computeStiffnessThreshold();

    //plot clusters
    std::vector<Point> getClusters();

    //plot contour points
    std::vector<Point> getContourPoints();
    
    SplineInterpolant_ptr_pair_vec getSplineInterpolant();
    //function only relevant for evalu
    std::string getResultsPath();
    bool writePosterior();
    void printParameters(const bayesopt::Parameters& par);
};

#endif