
#ifndef _MEANSHIFT_HPP_
#define _MEANSHIFT_HPP_
#include<mlpack/methods/mean_shift/mean_shift.hpp>
#include<mlpack/methods/kmeans/kmeans.hpp>
#include<vector>
#include"bayesopt/bayesopt.hpp"
#include"fileParser.hpp"

struct Point{
    bool merged;
    int id;
    double x;
    double y;
   
    Point(): merged(false),x(0), y(0) {}
    Point(double _x, double _y): x(_x), y(_y),merged(false) {}
    
    
    double distance(const Point& other) const {
        return std::sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y));
    }
};
class MeanShift
{
private:
    std::vector<std::vector<double>> data_;
    
    std::vector<Point> newPoints, newCenters;
    std::vector<Point> scattered_points_;
    std::vector<Point> centers;
    double bandwidth_;
    int samples_;
    arma::mat centroids_;

    bool pointNoPoint(float probabilityOfPoint);
public:
    MeanShift();
    MeanShift(std::vector<std::vector<double>> data, double bandwidth, int samples);
    ~MeanShift();
    void cluster();
    std::vector<Point> mergeClusters(double threshold);
    Point shiftPoint(const Point& p, const std::vector<Point>& data);
    void printClusters();
    bool saveDataToCSV();
    bool savePointsToCSV(std::string name, std::vector<Point>& points_);
    void scatterData(std::vector<Point>& points_);
    void meanshift_mlpack();
    std::vector<std::vector<double>> getCentroids();
};
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