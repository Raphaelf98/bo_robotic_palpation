#include<iostream>
#include<bayesopt/bayesopt.hpp>
#include<bayesopt/parameters.hpp>
#include <stdlib.h>
#include "param_loader.hpp"
#include "display2dgp.hpp"
#include "testfunctions.hpp"

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


int main(int nargs, char *args[])
{
  bayesopt::Parameters par;
  if(nargs > 1)
  {
    if(!bayesopt::utils::ParamLoader::load(args[1], par)){
        std::cout << "ERROR: provided file \"" << args[1] << "\" does not exist" << std::endl;
        return -1;
    }
    
  }
  else{
    par = initialize_parameters_to_default();
    
    par.n_iterations = 60;
    par.n_init_samples = 10;
    par.crit_name = "cEI";
    par.epsilon = 3;
    
    
    par.random_seed = 10;
    par.init_method = 3;
    
    //par.crit_name = "cLCB";
    //par.epsilon = 10;
    par.mean.name = "mZero";
    par.force_jump = 0;
    par.kernel.name = "kSEARD";
    par.kernel.hp_mean[0] = 0.08;
    par.kernel.hp_std[0] = 1.0;
    par.n_inner_iterations = 500;
    //set_surrogate(&par,"sStudentTProcessNIG");

    //par.l_type = L_MCMC;
    //par.sc_type = SC_MAP;

    
    par.verbose_level = 1;
    
  }
  bayesopt::utils::ParamLoader::save("triangle_display_params.txt", par);
  boost::scoped_ptr<Triangle> triangle(new Triangle(par));
  GLOBAL_MATPLOT.init(triangle.get(),2);
  //boost::scoped_ptr<BraninNormalized> Branin(new BraninNormalized(par));
  //GLOBAL_MATPLOT.init(Branin.get(),2);
  //boost::scoped_ptr<Circle> circle(new Circle(par));
  //GLOBAL_MATPLOT.init(circle.get(),2);

  vectord sv(2);  
  sv(0) = 0.1239; sv(1) = 0.8183;
  GLOBAL_MATPLOT.setSolution(sv);
  
  sv(0) = 0.5428; sv(1) = 0.1517;
  GLOBAL_MATPLOT.setSolution(sv);

  sv(0) = 0.9617; sv(1) = 0.1650;
  GLOBAL_MATPLOT.setSolution(sv);
  
  glutInit(&nargs, args);
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
