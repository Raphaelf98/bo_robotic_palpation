#include "helper.hpp"
/**
 * @file
 * 
 */



/**
* @brief Returns absolute file paths.
*/
std::string generateFilePath(const char* dir, const char* file)
{   
    std::filesystem::path p = std::filesystem::current_path();
    p = p.parent_path();
    std::string working_dir_path = p;
    std::string rel_data_path = dir;
    std::string scattered_data_file = file;
    std::string path = working_dir_path+rel_data_path+scattered_data_file;
    std::cout<<"PATH  "<<path<<std::endl;
    return path;
}
/**
* @brief Returns absolute file paths of experiment directories.
*/
std::string generateExperimentFilePath(std::string experiment_dir, const char* dir, const char* file)
{   

    std::string rel_data_path = dir;
    std::string scattered_data_file = file;
    std::string path = experiment_dir+rel_data_path+scattered_data_file;
    
    return path;
}
/**
* @brief Create or Overwirte directory
*/
bool createOrOverwriteDirectory(const std::string& path) {
    try {
        std::filesystem::path dirPath(path);

        // Check if directory exists
        if (std::filesystem::exists(dirPath) && std::filesystem::is_directory(dirPath)) {
            // Remove the existing directory and all its contents
            std::filesystem::remove_all(dirPath);
        }

        // Create a new directory
        if (std::filesystem::create_directory(dirPath)) {
            std::cout << "Directory created: " << path << std::endl;
            return true;
        } else {
            std::cout << "Failed to create directory: " << path << std::endl;
            return false;
        }
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return false;
    }
}
/**
* @brief Creates directory according to shape name.
*/
std::string createShapeDirectory(const std::string& parentDir, const std::string& shapeBaseName) {
    int dirNumber = 1;
    std::string dirPath;

    // Make sure the parent directory exists
    if (!std::filesystem::exists(parentDir) || !std::filesystem::is_directory(parentDir)) {
        std::cerr << "Parent directory does not exist or is not a directory: " << parentDir << std::endl;
        return "";
    }

    // Find the next available directory number
    do {
        dirPath = std::filesystem::path(parentDir) / (shapeBaseName + std::to_string(dirNumber));
        dirNumber++;
    } while (std::filesystem::exists(dirPath));

    // Create the new directory
    std::filesystem::create_directory(dirPath);

    std::cout << "Created directory: " << dirPath << std::endl;
    std::string path =  dirPath;
    return path;
}
/**
* @brief copyFileToDirectory
*/
bool copyFileToDirectory(const std::string& sourceFilePath, const std::string& destinationDirectory, const std::string& destinationFilename) {
    // Construct the full source and destination file paths
    std::filesystem::path srcFilePath = std::filesystem::path(sourceFilePath);
    std::filesystem::path destFilePath = std::filesystem::path(destinationDirectory) / destinationFilename;

    // Ensure the source file exists
    if (!std::filesystem::exists(srcFilePath) || !std::filesystem::is_regular_file(srcFilePath)) {
        std::cerr << "Source file does not exist or is not a regular file: " << srcFilePath << std::endl;
        return false;
    }

    // Ensure the destination directory exists
    if (!std::filesystem::exists(destinationDirectory)) {
        std::cerr << "Destination directory does not exist: " << destinationDirectory << std::endl;
        return false;
    }

    // Copy the file
    try {
        std::filesystem::copy_file(srcFilePath, destFilePath, std::filesystem::copy_options::overwrite_existing);
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error copying file: " << e.what() << std::endl;
        return false;
    }

    std::cout << "File copied successfully to " << destFilePath << std::endl;
    return true;
}
/**
* @brief saveMetricsToFile
*/
bool saveMetricsToFile(int pairNumber, double specificity, double sensitivity, const std::string& filePath) {
    std::ofstream file(filePath, std::ios::app); // Open in append mode

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        return false;
    }

    file << "Contour #" << pairNumber << std::endl;
    file << "Specificity: " << specificity << std::endl;
    file << "Sensitivity: " << sensitivity << std::endl << std::endl;

    file.close();
    return true;
}
/**
* @brief saveFileToCSV
*/
bool saveFileToCSV(const std::string& filename, const std::vector<std::vector<double>>& data) {
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
/**
* @brief readCoordinatesFromCSV
*/
std::vector<std::pair<double, double>> readCoordinatesFromCSV(const std::string& filePath) {
    std::vector<std::pair<double, double>> points;
    std::ifstream file(filePath);
    std::string line;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string x_str, y_str;
            double x, y;

            std::getline(ss, x_str, ',');
            std::getline(ss, y_str);

            x = std::stod(x_str);
            y = std::stod(y_str);

            points.push_back(std::make_pair(x, y));
        }
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filePath << std::endl;
    }

    return points;
}