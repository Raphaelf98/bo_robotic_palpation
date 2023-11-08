
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
    /*
    double x = 0.0;
    double y = 0.0;
    double increment_x = 1.0/(double) samples;
    double increment_y = 1.0/(double) samples;
    for(size_t i=0; i<samples; ++i)
	{       
        
          for (size_t j=0; j<samples; ++j)
          {
            data_->at(i).at(j) = data[i][j]; 
            Point point(x,y, data[i][j]);
            points.push_back(point);
            x += increment_x;
            std::cout<<"Point at: X:"<<point.x <<"  Y:"<<point.y<<" f:"<<point.f<< std::endl;
          }
          x=0;
        y += increment_y;
	}
    for(auto& row:*data_){
        std::cout<<"[";
        for(auto& col:row){
      std::cout<<col<< ", ";
    }
     std::cout<<"]"<<std::endl;
}
    */
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
    meanShift.Cluster(data, assignments, centroids_, forceConvergence);


    // Save the results, if needed.
    mlpack::data::Save("assignments.csv", assignments, true);
    mlpack::data::Save("centroids.csv", centroids_, true); 
    centroids_.print();
    
    
  
}
std::vector<std::vector<double>> MeanShift::getCentroids()
{
    std::vector<std::vector<double>> centroids{centroids_.n_rows, std::vector<double>(centroids_.n_cols)};
    for (size_t i = 0; i < centroids_.n_rows; i++) {
    for (size_t j = 0; j < centroids_.n_cols; j++) {
        centroids[i][j] = centroids_(i, j);
        
    }
    
}
return centroids;
}
MeanShift::~MeanShift()
{
  
    std::cout<<"data object removed from heap"<<std::endl;
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
Point MeanShift::shiftPoint(const Point& p, const std::vector<Point>& data) {
        Point shiftedPoint;
        double totalWeight = 0;
        
        for (const auto& point : data) {
            double distance = p.distance(point);

            if (distance < bandwidth_ ) {
                double weight = exp(-0.5*distance*distance);
                shiftedPoint.x += point.x * weight;
                shiftedPoint.y += point.y * weight;
                totalWeight += weight;
            }
        }
        shiftedPoint.x /= totalWeight;
        shiftedPoint.y /= totalWeight;
       
        //std::cout<<"new point is: x/y "<< shiftedPoint.x <<" / " <<shiftedPoint.y<<std::endl;
        return shiftedPoint;
    }
void MeanShift::cluster()
{
    std::cout <<"Cluster"<<  std::endl;
    std::vector<Point> oldPoints = scattered_points_;
    newPoints = std::vector<Point>(scattered_points_.size());
    double maxShiftDistance;

    do {
            std::cout <<"shift..."<<  std::endl;
            maxShiftDistance = 0;
            for (size_t i = 0; i < scattered_points_.size(); ++i) {
                newPoints[i] = shiftPoint(oldPoints[i], scattered_points_);
                double shiftDistance = oldPoints[i].distance(newPoints[i]);
                
                if (shiftDistance > maxShiftDistance) 
                {
                    maxShiftDistance = shiftDistance;
                }
            }
            oldPoints = newPoints;
        } while (maxShiftDistance > 0.004);  // You can adjust this threshold as needed
        savePointsToCSV("../config/centers.csv",newPoints);


}
std::vector<Point> MeanShift::mergeClusters(double threshold)
{   
    std::cout <<"Merge Clusters"<<  std::endl;
    centers= newPoints;
    
    for (size_t i = 0; i < centers.size(); ++i) {
        for (size_t j = i + 1; j < centers.size();) {   
            if (centers[i].distance(centers[j]) < threshold) {
                // Merge centers[j] into centers[i]
                centers[i].x = (centers[i].x + centers[j].x) / 2;
                centers[i].y = (centers[i].y + centers[j].y) / 2;
                
                // Remove centers[j]
                centers.erase(centers.begin() + j);
               
            } else {
                ++j; // Only increment if no merging was done
            }
        }
        newCenters = centers;
    }
    return newCenters;
    
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
    //data_.reshape(1,0);
    data_.resize(1,1);
    data_.print();
    //pre-fill centroids_ with guesses for initial clusters
    std::vector<double> cs = {0.5,0.6};
    centroids_ = arma::mat(cs);
    centroids_.resize(1,1);

}
void K_means::cluster(){
    kmeans_.Cluster(data_,2,assignments_,centroids_);
}
std::vector<double> K_means::getCentroids(){
    return std::vector<double>(centroids_.begin(), centroids_.end());
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