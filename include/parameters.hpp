#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_
#include<iostream>
#include"fileparser.hpp"
#include <string>


#define DATA_PATH "/data"
#define LOG_PATH "/log/"
#define CONFIG_PATH "/config/"
#define RESULTS_PATH "/results/"
#define PARAMS_PATH "/parameters/"
#define FILE_METRICS "metrics.txt"
#define FILE_BO_PARAMETERS "bo_parameters.txt"
#define FILE_CONTOUR_PARAMETERS "contour_parameters.txt"
#define FILE_TUMOR_MODEL_PARAMETERS "tumor_model_parameters.txt"
#define FILE_PARAMETERS_STORED "parameters_stored.txt"
#define FILE_SCATTERED_DATA "scattered_data.csv"
#define FILE_NORMALIZED_DATA "normalized_data.csv"
#define FILE_MS_ASSIGNMENTS "ms_assignments.csv"
#define FILE_MS_CENTROIDS "ms_centroids.csv"
#define FILE_CONTOUR_POINTS "contour_points_"
#define FILE_CONTOUR_APPRX "contour_"
#define FILE_EVAL_GROUND_TRUTH "groundTruth.csv"
#define FILE_POSTERIOR "posterior.txt"
#define FILE_STANDARDDEV "standard_deviation.txt"
#define FILE_GP_SAMPLES "gp_samples.csv"
#define FILE_GROUND_TRUTH_HEATMAP "groundTruthHeatMap.csv"
/**
 * @file
 * 
 */
/** 
 * @brief Tumor model parameters.
 *This struct holds parameters for triangle, rectangle, circle and two circle tumor models.
 * A default parameter setting is provided but parameters can be overwritten. 
 */
struct TumorModelParameters
{
    //Triangle Parameters
    double triangle_low=1;
    double triangle_high=2;
    double triangle_radius=0.1;
    double triangle_x_trans=0.5;
    double triangle_y_trans=0.5;
    double triangle_epsilon=0.1;
    double triangle_noise = 0.01;
    //Rectangle Parameters
    double rectangle_low=1;
    double rectangle_high=2;
    double rectangle_radius=0.15;
    double rectangle_x_trans=0.5;
    double rectangle_y_trans=0.5;
    double rectangle_epsilon=0.2;
    double rectangle_noise = 0.01;
    //Circle Parameters
    double circle_low=1;
    double circle_high=2;
    double circle_radius=0.1;
    double circle_x_trans=0.5;
    double circle_y_trans=0.5;
    double circle_epsilon=0.1;
    double circle_noise = 0.01;
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
    double two_circles_noise=0.01;
    /** 
 * @brief Read model from file and store them into TumorModelParameters object. Is called by loadModelParameters.
 *
 *@param[in] fp file parser object that hold information about parameter file.
 *@param[in] cp TumorModelParameters class object that parameters from file are written into.
 */
    void loadModel(bayesopt::utils::FileParser &fp, TumorModelParameters &cp);

      /** 
 * @brief Read model from file and store them into TumorModelParameters object. Calls loadModel routine. 
 * 
 *@param[in] filename absolute file path to parameter file.
 *@param[in] cp TumorModelParameters class object that parameters from file are written into.
 */
    bool loadModelParameters(std::string filename, TumorModelParameters &cp);
    /** 
 * @brief Prints current parameters.

 */
    void printParameters();
};
/** 
 * @brief Contour class parameters.
 *This struct holds parameters to control the behaviour of centroid computation and counter point search in Contour class.
 */
struct ContourParameters{
    size_t n_exploration_directions=10;
    size_t c_points=100;
    double means_shift_bandwidth=0.05;
    size_t lim_steps=1000;
    double threshold_multiplier=3.0;
    double tumor_stiffness_guess_low = 0.1;
    double tumor_stiffness_guess_high = 0.9;
    /** 
 * @brief Read contour parameters from file and store them into ContourParameters object. Is called by loadContour.
 *
 *@param[in] fp file parser object that hold information about parameter file.
 *@param[in] cp TumorModelParameters class object that parameters from file are written into.
 */
    void loadContour(bayesopt::utils::FileParser &fp, ContourParameters &cp);
     /** 
 * @brief Read contour parameters from file and store them into ContourParameters object. Calls loadContour routine. 
 * 
 *@param[in] filename absolute file path to parameter file.
 *@param[in] cp ContourParameters class object that parameters from file are written into.
 */
    bool loadContourParameters(std::string filename, ContourParameters &cp);
     /** 
 * @brief Prints current parameters.

 */
    void PrintParameters();

};

#endif