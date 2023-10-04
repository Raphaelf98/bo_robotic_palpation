#include<iostream>
#include<bayesopt/bayesopt.hpp>
#include<bayesopt/parameters.hpp>
#include <stdlib.h>
#include "param_loader.hpp"
#include "displaygp.hpp"
#include "testfunctions.hpp"
using namespace bayesopt;

class Example1D: public ContinuousModel
{
 public:
  
  Example1D(bayesopt::Parameters param):ContinuousModel(1,param) 
  {
     // My constructor 
  }
  double evaluateSample( const boost::numeric::ublas::vector<double> &query ) 
  {
     double x = query(0);
     
     return (x-5.0)*(x-5.0) + 20;
  }
  bool checkReachability( const boost::numeric::ublas::vector<double> &query )
  { 
     // My restrictions here 
     return true;
  }
};

class ExampleDiscrete: public ContinuousModel
{
 public:
  
  ExampleDiscrete(bayesopt::Parameters param):ContinuousModel(1,param) 
  {
     // My constructor 
  }
  double evaluateSample( const boost::numeric::ublas::vector<double> &query ) 
  {   
      
     double x = query(0);
     std::cout<<"Sample probe at: "<< x ;
     if(x <= 1)
     {
         return -1.0;
     }
     if(x > 1 && x <= 5)
     {
         return -5.0;
     }
     else {
      return -2.0;
     }
     
  }
  bool checkReachability( const boost::numeric::ublas::vector<double> &query )
  { 
     // My restrictions here 
     return true;
  }
};

// Unfortunately OpenGL functions require no parameters, so the object
// has to be global.
bayesopt::utils::DisplayProblem2D GLOBAL_MATPLOT;

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
  if(nargs > 1){
    if(!bayesopt::utils::ParamLoader::load(args[1], par)){
        std::cout << "ERROR: provided file \"" << args[1] << "\" does not exist" << std::endl;
        return -1;
    }
  }
  else{
    par = initialize_parameters_to_default();
    par.n_iterations = 100;
    par.n_init_samples = 2;
    par.n_iter_relearn = 1;
    par.random_seed = 10;
    //set_surrogate(&par,"sStudentTProcessNIG");

    par.l_type = L_MCMC;
    par.sc_type = SC_MAP;
    par.verbose_level = 1;
    //bayesopt::utils::ParamLoader::save("bo_branin_display.txt", par);
  }
  
  boost::scoped_ptr<BraninNormalized> branin(new BraninNormalized(par));
  GLOBAL_MATPLOT.init(branin.get(),2);

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
