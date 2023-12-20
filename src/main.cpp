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



int main(int argc, char* argv[])
{
  bayesopt::Parameters par;
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


  if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <Shape>\n";
        return 1;
    }

    std::string arg = argv[1];
    ShapeType type = getShapeType(arg);
   
   
    std::unique_ptr<Shape> shape;
    
    switch (type) {

        case SHAPE_CIRCLE:
            std::cout << "Running Experiment on Circle" << std::endl;
            shape = std::make_unique<SmoothCircle>(par,1,2, 0.1,0.5,0.5,0.1); 
            break;

        case SHAPE_TRIANGLE:
            std::cout << "Running Experiment on Triangle" << std::endl;
            shape = std::make_unique<Triangle>(par,1,2, 0.1,0.5,0.5,0.2); 
            break;

        case SHAPE_RECTANGLE:
            std::cout << "Running Experiment on Rectangle" << std::endl;
            shape = std::make_unique<Rectangle>(par,1,2, 0.1,0.5,0.5,0.1); 
            break;

        case SHAPE_TWOCIRCLES:
            std::cout << "Running Experiment on Two Circles" << std::endl;
            shape = std::make_unique<TwoCircles>(par,1,2, 0.05,0.1,0.2,0.7,0.8,0.3,0.1); 
            break;
            
        default:
            std::cout << "Unknown Shape: " << arg << std::endl;
    }
    
    
  Contour contour(shape.get(),10);
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
