#ifndef __HELPER_HPP_
#define __HELPER_HPP_
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <string>
#include <filesystem>
#include<fstream>
std::vector<std::vector<double>> convertUblasToStd(
    const std::vector<boost::numeric::ublas::vector<double>>& ublasVecs) ;

std::vector<boost::numeric::ublas::vector<double>> convertStdToUblas(
    const std::vector<std::vector<double>>& stdVecs);
std::string generateFilePath(const char* dir, const char* file);
bool createOrOverwriteDirectory(const std::string& path);
std::string createShapeDirectory(const std::string& parentDir, const std::string& shapeBaseName);
std::string generateExperimentFilePath(std::string experiment_dir, const char* dir, const char* file);
bool copyFileToDirectory(const std::string& sourceFilePath, const std::string& destinationDirectory, const std::string& destinationFilename);
bool saveMetricsToFile(double specificity, double sensitivity, const std::string& filePath) ;
#endif