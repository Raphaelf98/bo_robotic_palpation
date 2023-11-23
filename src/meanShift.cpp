
#include "meanShift.hpp"
#include "prob_distribution.hpp"
#include<iostream>
#include <fstream>
#include <random>
#include<math.h>

MeanShift::MeanShift(){}
MeanShift::MeanShift(std::vector<std::vector<double>> data, double bandwidth, int samples)
{
    bandwidth_ = bandwidth;
    samples_ = samples;
    data_ = data;
    
    std::vector<Point> scattered_points;
    std::vector<double> max_row;
    std::vector<double> min_row;
    for (auto& row :data_)
    {
        double max = *max_element(row.begin(), row.end());
        double min = *min_element(row.begin(), row.end());
        
        max_row.push_back(max);
        min_row.push_back(min);

    }
    double max_value = *max_element(max_row.begin(), max_row.end());
    double min_value = *min_element(min_row.begin(), min_row.end());
    for(size_t i=0; i<data_.size(); ++i)
	{       
        
        for (size_t j=0; j<data[0].size(); ++j)
        {
            data_[i][j] = (data[i][j]-min_value) / max_value;
        }
       
	}

    saveDataToCSV();
    scatterData(scattered_points);
    savePointsToCSV("/home/raphael/robolab/displaygp/config/scattered_data.csv",scattered_points);
    scattered_points_ = scattered_points;
   
}
void MeanShift::meanshift_mlpack()
{   
    //Load your dataset.
    arma::mat data;
    mlpack::data::Load("/home/raphael/robolab/displaygp/config/scattered_data.csv", data, true);
    bool forceConvergence = true; // Flag whether to force each centroid seed  to converge regardless of maxIterations.
    // Create a Mean Shift object. 
    // Perform Mean Shift clustering.
    arma::Row<size_t> assignments;  // Cluster assignments.

    mlpack::meanshift::MeanShift<> meanShift;
    std::cout<<"Compute clusters ..."<<std::endl;
    meanShift.Cluster(data, assignments, centroids_, forceConvergence, true);


    // Save the results, if needed.
    mlpack::data::Save("assignments.csv", assignments, true);
    mlpack::data::Save("centroids.csv", centroids_, true); 
    centroids_.print();
    
    
  
}
std::vector<std::vector<double>> MeanShift::getCentroids()
{
    std::vector<std::vector<double>> centroids{centroids_.n_rows, std::vector<double>(centroids_.n_cols)};
    for (size_t i = 0; i < centroids_.n_rows; i++) 
    {
    for (size_t j = 0; j < centroids_.n_cols; j++) 
    {
        centroids[i][j] = centroids_(i, j);
        
    }
    
}
return centroids;
}
MeanShift::~MeanShift()
{
  
}
bool MeanShift::pointNoPoint(float probabilityOfPoint)
{
  return rand()%100 < (probabilityOfPoint * 100);
}

void MeanShift::scatterData(std::vector<Point>& scattered_points)
{       
    double x = 0.0;
    double y = 0.0;
    double increment_x = 1.0/(double) samples_;
    double increment_y = 1.0/(double) samples_;
    for(size_t i=0; i<data_.size(); ++i)
	{       
        
        for (size_t j=0; j<data_[0].size(); ++j)
        {
            x += increment_x;
            if(pointNoPoint(1.0-data_[i][j]))
            {
                Point p(x,y);
                scattered_points.push_back(p);
            }
        }
        x=0;
        y += increment_y;
       
	}



}

void MeanShift::printClusters()
{ 
    std::cout <<"Centers: "<< newCenters.size() << std::endl;
    for (size_t i = 0; i < newCenters.size(); ++i)
    {   
        std::cout <<"X: "<< newCenters[i].x << "Y: " << newCenters[i].y << std::endl;
    }
}
bool MeanShift::saveDataToCSV()
{   
  std::ofstream file("/home/raphael/robolab/displaygp/config/normalized_data.csv");

    if (!file.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return false;
    }

    for (const auto& row : data_) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i != row.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();
    return true;
}

bool MeanShift::savePointsToCSV(std::string name, std::vector<Point>& points_)
{   
    std::ofstream file(name);

    if (!file.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return false;
    }

    // Optionally, you can write headers
    file << "x,y\n";

    for (const auto& point : points_) {
        file << point.x << "," << point.y  << "\n";
    }

    file.close();
    return true;
}




void loadData(std::string filename, std::vector<std::vector<double>> &data)
{
    FileParser fp(filename);
    fp.open(1);
    fp.read_stdvecOfvec("posterior", data);
}

bool saveToCSV(const std::string& filename, const std::vector<std::vector<double>>& data) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return false;
    }

    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i != row.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();
    return true;
}

K_means::K_means(const std::vector<double> &vals){
    data_ = arma::mat(vals);
    std::cout<<"DIMENSIONS: "<< data_.size()<< std::endl;
    //data_.reshape(1,0);
    data_.reshape(1,data_.size());
    //pre-fill centroids_ with guesses for initial clusters
    std::vector<double> cs = {0.1,0.9,0.5};
    centroids_ = arma::mat(cs);


}
void K_means::cluster()
{
    kmeans_.Cluster(data_,2,assignments_,centroids_);
    std::cout<<"DATA"<<std::endl;
    data_.print();
    std::cout<<"ASSIGNMENTS"<<std::endl;
    assignments_.print();
}
std::vector<double> K_means::getCentroids()
{
    return std::vector<double>(centroids_.begin(), centroids_.end());
}
std::vector<std::pair<double,size_t>> K_means::getAssignments()
{   std::vector<std::pair<double,size_t>> v( data_.size());
    for(size_t i = 0; i < data_.size(); i++)
    {   
        v[i].first = data_(0,i);
        v[i].second = assignments_[i];
        
    }
    return v;
}
/*
int main(){
    std::vector<std::vector<double>> mat;
    loadData("/home/raphael/robolab/displaygp/config/posterior.txt", mat);
    saveToCSV("/home/raphael/robolab/displaygp/config/posterior.csv", mat);
    MeanShift ms(mat, 0.05, 100);
    ms.meanshift_mlpack();
    //ms.cluster();
    
    //std::vector<Point>  clusterCenter = ms.mergeClusters(2);
    //ms.printClusters();
    
    //ms.savePointsToCSV("../config/clusterCenters.csv" ,clusterCenter);
    return 0;
}
*/