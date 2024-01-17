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
#define FILE_GROUND_TRUTH_HEATMAP "groundTruthHeatMap.csv"

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
    //Rectanlge Parameters
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
    void loadModel(bayesopt::utils::FileParser &fp, TumorModelParameters &cp);
    bool loadModelParameters(std::string filename, TumorModelParameters &cp);
    void printParameters();
};
struct ContourParamters{
    size_t n_exploration_directions=10;
    size_t c_points=100;
    double means_shift_bandwidth=0.05;
    size_t lim_steps=1000;
    double threshold_multiplier=3.0;
    void loadContour(bayesopt::utils::FileParser &fp, ContourParamters &cp);
    bool loadContourParameters(std::string filename, ContourParamters &cp);
    void PrintParameters();

};

#endif