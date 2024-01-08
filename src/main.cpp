#include<iostream>
#include<bayesopt/bayesopt.hpp>
#include<bayesopt/parameters.hpp>
#include <stdlib.h>
#include "param_loader.hpp"
#include "display2dgp.hpp"
#include "testfunctions.hpp"
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

void loadModel(bayesopt::utils::FileParser &fp, TumorModelParameters &cp)
{
    fp.readOrWrite("triangle_low", cp.triangle_low);
    fp.readOrWrite("triangle_high", cp.triangle_high);
    fp.readOrWrite("triangle_radius", cp.triangle_radius);
    fp.readOrWrite("triangle_x_trans", cp.triangle_x_trans);
    fp.readOrWrite("triangle_y_trans", cp.triangle_y_trans);
    fp.readOrWrite("triangle_epsilon", cp.triangle_epsilon);

    fp.readOrWrite("rectangle_low", cp.rectangle_low);
    fp.readOrWrite("rectangle_high", cp.rectangle_high);
    fp.readOrWrite("rectangle_radius", cp.rectangle_radius);
    fp.readOrWrite("rectangle_x_trans", cp.rectangle_x_trans);
    fp.readOrWrite("rectangle_y_trans", cp.rectangle_y_trans);
    fp.readOrWrite("rectangle_epsilon", cp.rectangle_epsilon);

    fp.readOrWrite("circle_low", cp.circle_low);
    fp.readOrWrite("circle_high", cp.circle_high);
    fp.readOrWrite("circle_radius", cp.circle_radius);
    fp.readOrWrite("circle_x_trans", cp.circle_x_trans);
    fp.readOrWrite("circle_y_trans", cp.circle_y_trans);
    fp.readOrWrite("circle_epsilon", cp.circle_epsilon);

    fp.readOrWrite("two_circles_low", cp.two_circles_low);
    fp.readOrWrite("two_circles_high", cp.two_circles_high);
    fp.readOrWrite("two_circles_radius_1", cp.two_circles_radius_1);
    fp.readOrWrite("two_circles_radius_2", cp.two_circles_radius_2);
    fp.readOrWrite("two_circles_x_trans_1", cp.two_circles_x_trans_1);
    fp.readOrWrite("two_circles_x_trans_2", cp.two_circles_x_trans_2);
    fp.readOrWrite("two_circles_y_trans_1", cp.two_circles_y_trans_1);
    fp.readOrWrite("two_circles_y_trans_2", cp.two_circles_y_trans_2);
    fp.readOrWrite("two_circles_epsilon", cp.two_circles_epsilon);

    

}
void loadContour(bayesopt::utils::FileParser &fp, ContourParamters &cp)
{
    fp.readOrWrite("n_exploration_directions", cp.n_exploration_directions);
    fp.readOrWrite("c_points", cp.c_points);
    fp.readOrWrite("means_shift_bandwidth", cp.means_shift_bandwidth);
    fp.readOrWrite("lim_steps", cp.lim_steps);
    fp.readOrWrite("threshold_multiplier", cp.threshold_multiplier);
}

bool loadModelParameters(std::string filename, TumorModelParameters &cp)
{
    bayesopt::utils::FileParser fp(filename);
    if(!fp.fileExists())
    {
        return false;
    }
    
    fp.openInput();
    loadModel(fp,cp);
    return true;

}

bool loadContourParameters(std::string filename, ContourParamters &cp)
{
    bayesopt::utils::FileParser fp(filename);
    if(!fp.fileExists())
    {
        return false;
    }
    
    fp.openInput();
    loadContour(fp,cp);
    return true;
}

int main(int argc, char* argv[])
{
  bayesopt::Parameters par;
  TumorModelParameters model_parameters;
  ContourParamters contour_parameters;
  //Either load bayesian optimization parameters from file or set them manually
  //TODO adjust file parsing
  if (bayesopt::utils::ParamLoader::load("/home/raphael/robolab/displaygp/config/bo_parameters.txt", par))
  {
      std::cout << "Found bo_parameters.txt" << std::endl;
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
    
  }
  bayesopt::utils::ParamLoader::save("triangle_display_params.txt", par);

  if (loadContourParameters("/home/raphael/robolab/displaygp/config/contour_parameters.txt", contour_parameters))
  {
      std::cout << "Found contour_parameters.txt" << std::endl;
  }
  else
  {
    std::cout<<"Could not load contour_parameters.txt"<< std::endl;
  }

  if (loadModelParameters("/home/raphael/robolab/displaygp/config/tumor_model_parameters.txt", model_parameters))
  {
      std::cout << "Found tumor_model_parameters.txt" << std::endl;
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
  // std::string arg = "Triangle";
    ShapeType type = getShapeType(arg);
   
   
    std::unique_ptr<Shape> shape;
    
    switch (type) {

        case SHAPE_CIRCLE:
            std::cout << "Running Experiment on Circle" << std::endl;
            shape = std::make_unique<SmoothCircle>(par,model_parameters.circle_low,model_parameters.circle_high,model_parameters.circle_radius,
                                                    model_parameters.circle_x_trans,model_parameters.circle_y_trans,model_parameters.circle_epsilon); 
            break;

        case SHAPE_TRIANGLE:
            std::cout << "Running Experiment on Triangle" << std::endl;
            shape = std::make_unique<Triangle>(par,model_parameters.triangle_low,model_parameters.triangle_high,model_parameters.triangle_radius,
                                                    model_parameters.triangle_x_trans,model_parameters.triangle_y_trans,model_parameters.triangle_epsilon); 
            break;

        case SHAPE_RECTANGLE:
            std::cout << "Running Experiment on Rectangle" << std::endl;
            shape = std::make_unique<Rectangle>(par,model_parameters.rectangle_low,model_parameters.rectangle_high,model_parameters.rectangle_radius,
                                                    model_parameters.rectangle_x_trans,model_parameters.rectangle_y_trans,model_parameters.rectangle_epsilon); 
            break;

        case SHAPE_TWOCIRCLES:
            std::cout << "Running Experiment on Two Circles" << std::endl;
            shape = std::make_unique<TwoCircles>(par,model_parameters.two_circles_low,model_parameters.two_circles_high,model_parameters.two_circles_radius_1,model_parameters.two_circles_radius_2,
                                                    model_parameters.two_circles_x_trans_1,model_parameters.two_circles_x_trans_2,  model_parameters.two_circles_y_trans_1, model_parameters.two_circles_y_trans_2,
                                                            model_parameters.rectangle_epsilon);
            break;
            
        default:
            std::cout << "Unknown Shape: " << arg << std::endl;
    }
    
    
  Contour contour(shape.get(),contour_parameters);
  // Contour contour(shape.get(),10);
    GLOBAL_MATPLOT.init(&contour,2);
      
  glutInit(&argc, argv);
  glutCreateWindow(50,50,800,650);
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
