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

#define SPLINE_SAMPLES 1000

/**
 * @brief This class handles the methods to perform 
 *  1. Approximation through a gaussian process(GP) of an underlying tumor model
 *  2. Computation of maximum stiffness regions on posterior of GP in order to find tumor centroids
 *  3.  Computation of points attributed to tumor contour by exploration of Centroids 
 *  4.  Spline approximation of tumor contour
 *
 * 
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
    /**
    * @brief Contour class constructor.
    *
    * @param[in] bopt_model* BayesOptBase model pointer to run optimization on.
    * @param[in] cp Contour parameters used to define functionality of Contour class.
    * @param[in] experiment_path Current working directory path as string.
    * 
    * 
    */
    Contour(bayesopt::BayesOptBase* bopt_model,ContourParamters cp, std::string experiment_path);
    Contour(){}
    ~Contour();

    //Methods to perform external computation used by DisplayHeatMap2D class
    /**
    * @brief  Returns granularity of two dimensional gaussian process.
    *
    * 
    *
    * @retval Number of points in x and y direction.
    * 
    */
    size_t getCPoints();
    /**
    * @brief Get last sample from Gaussian process.
    *
    * 
    *
    * @retval Returns vector of last sample.
    * 
    */
    bayesopt::vectord getLastSample();
    /**
    * @brief Get current posterior distribution.
    *
    * @param[in] q point at which posterior distribution should be evaluated.
    *
    * @retval Distribution at point q.
    * 
    */
    bayesopt::ProbabilityDistribution* getPredictionGaussianProcess(const vectord &q);
    /**
    * @brief Get initial samples from initialization routine.
    *
    * @param[in] &samples_x Container for x-coordinates of samples.
    * @param[in] &samples_y Container for y-coordinates of samples.
    * 
    * 
    */
    void getInitialSamples( std::vector<double> &samples_x, std::vector<double> &samples_y );

    /**
    * @brief Evaluates ground truth at desired point.
    *
    * @param[in] q Coordinate at which ground truth should be evaluated.
    *
    * @retval returns stiffness value.
    * 
    */
    double evaluateGaussianProcess(const bayesopt::vectord &q);

    /**
    * @brief Get value of criterion (Expected Improvement) at desired point.
    *
    * @param[in] q Coordinate at which criterion should be evaluated.
    *
    * @retval returns value of acquisition function.
    * 
    */
    double evaluateCriteriaGaussianProcess(const bayesopt::vectord &q);

    /**
    * @brief Wrapper for BayesOpt initializeOptimizaition() function
    * 
    */
    void prepareGaussianProcess();

    /**
    * @brief Call of this function performs single update step of Bayesian optimization.
    */
    void stepRunGaussianProcess();
    /**
    * @brief Get clusters found by means-shift clustering algorithm.
    *
    *
    * @retval Vector of Points with x and y coordinate.
    * 
    */
    std::vector<Point> getClusters();//plot clusters
    /**
    * @brief Get contour points found by centroid exploration routine.
    *
    *
    * @retval Vector of Points of current Contour.
    * 
    */
    std::vector<Point> getContourPoints();   //plot contour points
    /**
    * @brief Get pair of spline interpolants in x and y.
    *
    *
    * @retval Returns pair of spline interpolants that can be used to query parametric contour.
    * 
    */
    SplineInterpolant_ptr_pair_vec getSplineInterpolant();

    //function only relevant for evaluation
    /**
    * @brief Returns path where results are stored.
    *
    *
    * @retval Absolute path.
    * 
    */
    std::string getResultsPath();

    //Methods to perform internal computation
    /**
    *  @brief Wrapper to perform single step optimization. 
    *
    * 
    *
    * 
    * 
    */
    void runGaussianProcess();
    /**
    * @brief Wrapper to call clustering routine.
    *
    *
    *
    * 
    * 
    */
    void computeCluster();
    /**
    * @brief Wrapper to call Centroid exploration routine.
    * 
    */
    void exploreContour();
    /**
    * @brief Perform approximation on Contour points.
    */
    void approximateContour();
    /**
    * @brief Computes stiffness threshold used to detect contour points.
    * 
    */
    void computeStiffnessThreshold();
   
    /**
    * @brief Log parameters used for optimization.
    *
    * @param[in] par BayesOpt paramters.
    *
    * 
    * 
    */
    void printParameters(const bayesopt::Parameters& par);
};

#endif