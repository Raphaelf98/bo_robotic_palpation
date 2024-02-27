#include<iostream>
#include<bayesopt/bayesopt.hpp>
#include<bayesopt/parameters.hpp>
#include <stdlib.h>

#include"param_loader.hpp"
#include "display2dgp.hpp"
#include "testfunctions.hpp"
#include "helper.hpp"
#include "contour.hpp"
using namespace bayesopt;


// Unfortunately OpenGL functions require no parameters, so the object
// has to be global.
DisplayHeatMap2D GLOBAL_MATPLOT;

void display( void ){ GLOBAL_MATPLOT.display(); }
void reshape( int w,int h ){ GLOBAL_MATPLOT.reshape(w,h); }
void idle( void ) { glutPostRedisplay(); } 

void mouse(int button, int state, int x, int y ){ GLOBAL_MATPLOT.mouse(button,state,x,y); }
void motion(int x, int y ){ GLOBAL_MATPLOT.motion(x,y); }
void passive(int x, int y ){ GLOBAL_MATPLOT.passivemotion(x,y); }

void keyboard(unsigned char key, int x, int y)
{
    GLOBAL_MATPLOT.keyboard(key, x, y); 
    if(key=='r')   //Toogle run/stop
    { 
	    GLOBAL_MATPLOT.toogleRUN();
    }
    if(key=='s')   //Activate one step
    { 
	    GLOBAL_MATPLOT.setSTEP();
    }
 
}
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
int main(int argc, char* argv[])
{

  bayesopt::Parameters par;
  TumorModelParameters model_parameters;
  ContourParamters contour_parameters;
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
    
    if (argc != 2) 
    {
        std::cerr << "Usage: " << argv[0] << " <Shape>\n";
        return 1;
    }

    std::string arg = argv[1];
    //std::string arg = "TwoCircles";
    
    std::string experiment_path = generateFilePath(DATA_PATH,"");
    std::string dirName = createShapeDirectory(experiment_path,arg);
    std::string log_dir = LOG_PATH;
    log_dir = dirName + log_dir;
    
    createOrOverwriteDirectory(log_dir);
    std::string results_dir = RESULTS_PATH;
    results_dir = dirName + results_dir;
    createOrOverwriteDirectory(results_dir);

    


    std::string contour_params_file = FILE_CONTOUR_PARAMETERS;
    std::string model_params_file = FILE_TUMOR_MODEL_PARAMETERS;
    std::string params_dir = PARAMS_PATH;
    params_dir = dirName + params_dir;
    createOrOverwriteDirectory(params_dir);
    
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
            shape = std::make_unique<Circle>(par,model_parameters.circle_low,model_parameters.circle_high,model_parameters.circle_radius,
                                                    model_parameters.circle_x_trans,model_parameters.circle_y_trans,model_parameters.circle_epsilon,model_parameters.circle_noise); 
           
            break;

        case SHAPE_TRIANGLE:
            std::cout << "Running Experiment on Triangle" << std::endl;
            shape = std::make_unique<Triangle>(par,model_parameters.triangle_low,model_parameters.triangle_high,model_parameters.triangle_radius,
                                                    model_parameters.triangle_x_trans,model_parameters.triangle_y_trans,model_parameters.triangle_epsilon, model_parameters.triangle_noise); 
            
            break;

        case SHAPE_RECTANGLE:
            std::cout << "Running Experiment on Rectangle" << std::endl;
            shape = std::make_unique<Rectangle>(par,model_parameters.rectangle_low,model_parameters.rectangle_high,model_parameters.rectangle_radius,
                                                    model_parameters.rectangle_x_trans,model_parameters.rectangle_y_trans,model_parameters.rectangle_epsilon, model_parameters.rectangle_noise); 
            
            break;

        case SHAPE_TWOCIRCLES:
            std::cout << "Running Experiment on Two Circles" << std::endl;
            shape = std::make_unique<TwoCircles>(par,model_parameters.two_circles_low,model_parameters.two_circles_high,model_parameters.two_circles_radius_1,model_parameters.two_circles_radius_2,
                                                    model_parameters.two_circles_x_trans_1,model_parameters.two_circles_x_trans_2, model_parameters.two_circles_y_trans_1, model_parameters.two_circles_y_trans_2,
                                                            model_parameters.two_circles_epsilon, model_parameters.two_circles_noise);
            
            break;
            
        default:
            std::cout << "Unknown Shape: " << arg << std::endl;
    }
    
    shape->saveGroundTruth(contour_parameters.c_points,log_dir);

  Contour contour(shape.get(),contour_parameters, dirName);
   //Contour contour(shape.get(),10);
    GLOBAL_MATPLOT.init(&contour,2);
      
  glutInit(&argc, argv);
  glutCreateWindow(50,50,900,650);
  glutDisplayFunc( display );
  glutReshapeFunc( reshape );
  glutIdleFunc( idle );
  glutMotionFunc( motion );
  glutMouseFunc( mouse );
  glutPassiveMotionFunc(passive);    
  glutKeyboardFunc( keyboard );     
  
  glutMainLoop();    
  

  return 0;
}
