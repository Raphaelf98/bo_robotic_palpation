
#ifndef _MEANSHIFT_HPP_
#define _MEANSHIFT_HPP_
#include<mlpack/methods/mean_shift/mean_shift.hpp>
#include<mlpack/methods/kmeans/kmeans.hpp>
#include<vector>
#include"bayesopt/bayesopt.hpp"
#include"fileparser.hpp"
#include"parameters.hpp"
#include "helper.hpp"

struct Point
{
    double x;
    double y;
    Point(): x(0), y(0) {}
    Point(double _x, double _y): x(_x), y(_y) {}
};
/*
This class serves as a wrapper for mlpack's implementation of the means shift clustering algorithm.
Source: https://github.com/mlpack/mlpack/blob/master/src/mlpack/methods/mean_shift/mean_shift.hpp
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
    bool pointNoPoint(float probabilityOfPoint);
public:
    MeanShift();
    MeanShift(std::vector<std::vector<double>> data, double bandwidth, int samples, std::string experiment_path);
    ~MeanShift();
    void printClusters();
    bool saveDataToCSV();
    bool savePointsToCSV(std::string name, std::vector<Point>& points_);
    void scatterData(std::vector<Point>& points_);
    void meanshift_mlpack();
    std::vector<std::vector<double>> getCentroids();
};
/*
This class serves as a wrapper for mlpack's implementation of the means shift clustering algorithm.
Source: https://github.com/mlpack/mlpack/blob/master/src/mlpack/methods/mean_shift/mean_shift.hpp
*/
class K_means{
    public:
    K_means(const std::vector<double> &vals);
    void cluster();
    std::vector<double> getCentroids();
    std::vector<std::pair<double,size_t>> getAssignments();
    private:
    arma::mat data_;
    size_t clusters_;
    arma::Row<size_t> assignments_;
    mlpack::kmeans::KMeans<> kmeans_;
    arma::mat centroids_;
};
#endif