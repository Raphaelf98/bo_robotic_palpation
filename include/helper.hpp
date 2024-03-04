#ifndef __HELPER_HPP_
#define __HELPER_HPP_
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <string>
#include <filesystem>
#include<numeric>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include<fstream>

#include "stdafx.h"
#include "interpolation.h"
/**
 * @file
 * 
 */

std::string generateFilePath(const char* dir, const char* file);
bool createOrOverwriteDirectory(const std::string& path);
std::string createShapeDirectory(const std::string& parentDir, const std::string& shapeBaseName);
std::string generateExperimentFilePath(std::string experiment_dir, const char* dir, const char* file);
bool copyFileToDirectory(const std::string& sourceFilePath, const std::string& destinationDirectory, const std::string& destinationFilename);
bool saveMetricsToFile(int pairNumber,double specificity, double sensitivity, const std::string& filePath) ;
bool saveFileToCSV(const std::string& filename, const std::vector<std::vector<double>>& data);

/**
 * @brief Type conversion std::vector<boost::numeric::ublas::vector<double>> to std::vector<std::vector<double>>

*/
inline std::vector<std::vector<double>> convertUblasToStd(
    const std::vector<boost::numeric::ublas::vector<double>>& ublasVecs) 
{
    std::vector<std::vector<double>> stdVecs;
    stdVecs.reserve(ublasVecs.size());

    for (const auto& ublasVec : ublasVecs) {
        std::vector<double> tempVec(ublasVec.size());
        std::copy(ublasVec.begin(), ublasVec.end(), tempVec.begin());
        stdVecs.push_back(std::move(tempVec));
    }

    return stdVecs;
}
/**
 * @brief Type conversion std::vector<std::vector<double>> to std::vector<boost::numeric::ublas::vector<double>> 

*/
inline std::vector<boost::numeric::ublas::vector<double>> convertStdToUblas(
    const std::vector<std::vector<double>>& stdVecs) 
{
    std::vector<boost::numeric::ublas::vector<double>> ublasVecs;
    ublasVecs.reserve(stdVecs.size());

    for (const auto& stdVec : stdVecs) {
        boost::numeric::ublas::vector<double> tempVec(stdVec.size());
        std::copy(stdVec.begin(), stdVec.end(), tempVec.begin());
        ublasVecs.push_back(std::move(tempVec));
    }

    return ublasVecs;
}
    
inline std::vector<double> linSpace(double min, double max,int n)
{
    std::vector<double> a;
    if(n<1){n=1;}
    a.resize(n);
    for(int i=0;i<n;++i){a[i]=min+(max-min)*i/(n-1);}
    return a;
};

inline double stdDev(std::vector<double> v)
{
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double m =  sum / v.size();

    double accum = 0.0;
    std::for_each (v.begin(), v.end(), [&](const double d) {
        accum += (d - m) * (d - m);
    });

    return sqrt(accum / (v.size()));

};
std::vector<std::pair<double, double>> readCoordinatesFromCSV(const std::string& filePath);

void loadDataFP(std::string filename, std::vector<std::vector<double>> &data);
#endif