
#ifndef _MEANSHIFT_HPP_
#define _MEANSHIFT_HPP_
#include<mlpack/methods/mean_shift/mean_shift.hpp>
#include<mlpack/methods/kmeans/kmeans.hpp>
#include<vector>
#include"bayesopt/bayesopt.hpp"
#include"fileparser.hpp"
#include"parameters.hpp"
#include "helper.hpp"
/**
 * @brief This struct holds x and y coordinates.
 *
 * 
 */
struct Point
{
    double x;
    double y;
    Point(): x(0), y(0) {}
    Point(double _x, double _y): x(_x), y(_y) {}
};
/**
* @brief This class serves as a wrapper for mlpack's implementation of the means shift clustering algorithm.
*      Source: https://github.com/mlpack/mlpack/blob/master/src/mlpack/methods/mean_shift/mean_shift.hpp
 * 
 */
class MeanShift
{
private:

    std::string experiment_path_;
    std::vector<std::vector<double>> data_;
    
    std::vector<Point> newPoints, newCenters;
    std::vector<Point> scattered_points_;
    std::vector<Point> centers;
    double bandwidth_;
    int samples_;
    arma::mat centroids_;
    std::string working_dir_path_;
    /**
    * @brief Function used to produce scattered data from two dimensional posterior distribution.
    * @param[in] probabilityOfPoint value between 0 and 1 from normalized posterior distribution.
    * @retval Point or no point.
    * 
    */
    bool pointNoPoint(float probabilityOfPoint);

public:
    /**
    * @brief MeanShift default class constructor.
    *
    * 
    */
    MeanShift();
    /**
    * @brief MeanShift class constructor.
    *
    * @param[in] data Posterior distribution to run means shift on.
    * @param[in] bandwidth Means shift bandwidth parameter.
    * @param[in] experiment_path Current working directory path as string.
    * 
    * 
    */
    MeanShift(std::vector<std::vector<double>> data, double bandwidth, int samples, std::string experiment_path);
    ~MeanShift();
    /**
    * @brief Print clusters.
    *
    * 
    */
    void printClusters();
     /**
    * @brief Save scattered data to CSV.
    *
    * @retval Returns true if successful.
    */
    bool saveDataToCSV();
     /**
    * @brief Save Centroids to CSV.
    * @param[in] name Filename
    * @param[in] points_ Centroids
    * @retval Returns true if successful.
    */
    bool savePointsToCSV(std::string name, std::vector<Point>& points_);
    /**
    * @brief Convert posterior distribution to scattered data representation. 
    * @param[in] points_ holds scattered points. 
    * 
    */
    void scatterData(std::vector<Point>& points_);
    /**
    * @brief Run the means-shift algorithm.
    * 
    */
    void meanshift_mlpack();
    /**
    * @brief Retrieve centroids.
    *  @retval Vector of centroid points
    */
    std::vector<std::vector<double>> getCentroids();
};

/**
* @brief This class serves as a wrapper for mlpack's implementation of the k-means clustering algorithm.
*      Source: https://github.com/mlpack/mlpack/blob/master/src/mlpack/methods/mean_shift/mean_shift.hpp
 * 
 */
class K_means{
    public:
    /**
    * @brief Contour class constructor.
    * @param[in] bopt_model* BayesOptBase model pointer to run optimization on.
    * @param[in] cp Contour parameters used to define functionality of Contour class.
    * @param[in] cp Contour parameters used to define functionality of Contour class.
    * 
    */
    K_means(const std::vector<double> &vals,const double low_stiffness,const double high_stifness);
    /**
    * @brief Contour class constructor.
    * 
    */
    void cluster();
        /**
    * @brief Contour class constructor.
    * @param[in] bopt_model* BayesOptBase model pointer to run optimization on.
    * @param[in] cp Contour parameters used to define functionality of Contour class.
    * 
    */
    std::vector<double> getCentroids();
        /**
    * @brief Contour class constructor.
    * @param[in] bopt_model* BayesOptBase model pointer to run optimization on.
    * @param[in] cp Contour parameters used to define functionality of Contour class.
    * 
    */
    std::vector<std::pair<double,size_t>> getAssignments();
    private:
    arma::mat data_;
    size_t clusters_;
    arma::Row<size_t> assignments_;
    mlpack::kmeans::KMeans<> kmeans_; 
    
    arma::mat centroids_;
};
#endif