#include<iostream>
#include<bayesopt/bayesopt.hpp>
#include<bayesopt/parameters.hpp>
#include <stdlib.h>

#include"param_loader.hpp"
#include "testfunctions.hpp"
#include "helper.hpp"
#include "contour.hpp"
using namespace alglib;

enum ShapeType {
    SHAPE_CIRCLE,
    SHAPE_TRIANGLE,
    SHAPE_RECTANGLE,
    SHAPE_TWOCIRCLES,
    SHAPE_UNKNOWN // for unrecognized strings
};
ShapeType getShapeType(const std::string& shape) {
    if (shape == "Circle") return SHAPE_CIRCLE;
    if (shape == "Triangle") return SHAPE_TRIANGLE;
    if (shape == "Rectangle") return SHAPE_RECTANGLE;
    if (shape == "TwoCircles") return SHAPE_TWOCIRCLES;
    return SHAPE_UNKNOWN;
}

int main(int argc, char *argv[])
{
   
  bayesopt::Parameters par;
  TumorModelParameters model_parameters;
  ContourParameters contour_parameters;
  //Either load bayesian optimization parameters from file or set them manually
  //TODO adjust file parsing
  std::string config_path = generateFilePath(CONFIG_PATH,"");
  std::string bo_parameters = FILE_BO_PARAMETERS;
  bo_parameters = config_path+bo_parameters;
  if (bayesopt::utils::ParamLoader::load(bo_parameters, par))
  {
      std::cout << "Found "+ bo_parameters << std::endl;
  }
  else
  {
    par = initialize_parameters_to_default();
    
    par.n_iterations = 60;
    par.n_init_samples = 10;
    par.crit_name = "cEI";
    par.epsilon = 3;
    par.random_seed = 10;
    par.init_method = 3;
    par.mean.name = "mZero";
    par.force_jump = 0;
    par.kernel.name = "kSEARD";
    par.kernel.hp_mean[0] = 0.08;
    par.kernel.hp_std[0] = 1.0;
    par.n_inner_iterations = 500; 
    par.verbose_level = 1;
     std::cout << "Using default paramters" << std::endl;
  }
  std::string contour_params_path = FILE_CONTOUR_PARAMETERS;
  contour_params_path = config_path+contour_params_path;
  if (contour_parameters.loadContourParameters(contour_params_path, contour_parameters))
  {
      std::cout << "Found contour_parameters.txt" << std::endl;
      contour_parameters.PrintParameters();
      
  }
  else
  {
    std::cout<<"Could not load contour_parameters.txt"<< std::endl;
  }
  std::string tumor_params_path = FILE_TUMOR_MODEL_PARAMETERS;
  tumor_params_path = config_path+tumor_params_path;
  if (model_parameters.loadModelParameters(tumor_params_path, model_parameters))
  {
      std::cout << "Found tumor_model_parameters.txt" << std::endl;
       model_parameters.printParameters();
  }
  else
  {
    std::cout<<"Could not load tumor_model_parameters"<< std::endl;
  }

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <Shape>\n";
        return 1;
    }

    std::string arg = argv[1];
    std::string experiment_path = generateFilePath(DATA_PATH,"");
    std::string dirName = createShapeDirectory(experiment_path,arg);
    std::string log_dir = LOG_PATH;
    log_dir = dirName + log_dir;
    createOrOverwriteDirectory(log_dir);
    std::string results_dir = RESULTS_PATH;
    results_dir = dirName + results_dir;
    createOrOverwriteDirectory(results_dir);
    //std::string arg = "TwoCircles";


    std::string contour_params_file = FILE_CONTOUR_PARAMETERS;
    
    std::string model_params_file = FILE_TUMOR_MODEL_PARAMETERS;
    std::string params_dir = PARAMS_PATH;
    params_dir = dirName + params_dir;
    createOrOverwriteDirectory(params_dir);
    std::cout<<"params_dir: "<<params_dir<<std::endl;
    //model_params_file = log_dir + model_params_file;
    std::string bo_paramters = FILE_BO_PARAMETERS;
    bo_parameters = params_dir+bo_paramters;
    bayesopt::utils::ParamLoader::save(bo_parameters,par);
    copyFileToDirectory(contour_params_path, params_dir, contour_params_file);
    copyFileToDirectory(tumor_params_path, params_dir, model_params_file);

    ShapeType type = getShapeType(arg);
   
    
    std::unique_ptr<Shape> shape;
    switch (type) {

        case SHAPE_CIRCLE:
            std::cout << "Running Experiment on Circle" << std::endl;
             
            shape = std::make_unique<Circle>(par,model_parameters); 
           
            break;

        case SHAPE_TRIANGLE:
            std::cout << "Running Experiment on Triangle" << std::endl;
            
            shape = std::make_unique<Triangle>(par,model_parameters);                                        
            
            break;

        case SHAPE_RECTANGLE:
            std::cout << "Running Experiment on Rectangle" << std::endl;
            
            shape = std::make_unique<Rectangle>(par,model_parameters);
            break;

        case SHAPE_TWOCIRCLES:
            std::cout << "Running Experiment on Two Circles" << std::endl;
            
            shape = std::make_unique<TwoCircles>(par,model_parameters);
            break;
            
        default:
            std::cout << "Unknown Shape: " << arg << std::endl;
    }
    
    shape->saveGroundTruth(contour_parameters.c_points,log_dir);
    Contour contour(shape.get(),contour_parameters, dirName);
    //run bayesian optimization
    contour.runGaussianProcess();
    contour.computeCluster();
    //contour.exploreContour();
    //contour.approximateContour();
}