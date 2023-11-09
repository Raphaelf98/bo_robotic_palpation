#include<iostream>
#include<bayesopt/bayesopt.hpp>
#include<bayesopt/parameters.hpp>
#include <stdlib.h>
#include "param_loader.hpp"
#include "displaygp.hpp"
#include "testfunctions.hpp"
#include<thread>
#include<mutex>
#include <condition_variable>
#include <chrono>
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
     
     if(x <= 0.3)
     {
         return -1.0;
     }
     if (x > 0.3 && x <= 0.4)
     {
        return evaluateSmoothStep(0.3,0.4,-1.0,-3.0,x);
     }
     if(x > 0.4 && x <= 0.5)
     {
         return -3.0;
     }
     if (x > 0.5 && x <= 0.6)
     {
        return evaluateSmoothStep(0.5,0.6,-3.0,-2.0,x);
     }
     else {
      return -2.0;
     }
     
  }
  double evaluateSmoothStep(const double x_l, const double x_r, const double y_l, const double y_r, const double x){
    double a = (-2*y_l + 2*y_r)/(x_l*x_l*x_l - 3*x_l*x_l*x_r + 3*x_l*x_r*x_r - x_r*x_r*x_r);
    double b = (3*x_l*y_l - 3*x_l*y_r + 3*x_r*y_l - 3*x_r*y_r)/(x_l*x_l*x_l - 3*x_l*x_l*x_r + 3*x_l*x_r*x_r - x_r*x_r*x_r);
    double c = (-6*x_l*x_r*y_l + 6*x_l*x_r*y_r)/(x_l*x_l*x_l - 3*x_l*x_l*x_r + 3*x_l*x_r*x_r - x_r*x_r*x_r);
    double d = (x_l*x_l*x_l*y_r - 3*x_l*x_l*x_r*y_r + 3*x_l*x_r*x_r*y_l - x_r*x_r*x_r*y_l)/(x_l*x_l*x_l - 3*x_l*x_l*x_r + 3*x_l*x_r*x_r - x_r*x_r*x_r);
    return a*x*x*x + b*x*x + c*x + d;
  }
  bool checkReachability( const boost::numeric::ublas::vector<double> &query )
  { 
     // My restrictions here 
     return true;
  }
};

// Unfortunately OpenGL functions require no parameters, so the object has to be global.
bayesopt::utils::DisplayProblem1D GLOBAL_MATPLOT;

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

int menu()
{
  std::string input;
  int option = 0;
  while ((option < 1) || (option > 5))
    {
      std::cout << "Please select an option for the parameters:\n\n"
		<< "  1- Default parameters.\n"
		<< "  2- Student t process model.\n"
		<< "  3- Combined kernel.\n"
		<< "  4- Lower Confidence Bound.\n"
		<< "  5- A-optimality criteria.\n\n"
		<< "Select [1-5]>";

      std::cin >> input;
      std::istringstream is(input);
      is >> option;
    }
  return option;
}
void dostuff(int nargs, char *args[])
{
  bopt_params parameters = initialize_parameters_to_default();
  parameters.n_init_samples = 7;
  parameters.n_iterations = 100;
  parameters.verbose_level = 2;

  switch( menu() )
    {
    case 1: break;
    case 2: 
      {
	set_surrogate(&parameters,"sStudentTProcessNIG"); 
	parameters.n_iter_relearn = 5;
	break;
      }
    case 3:   
      { 
	set_kernel(&parameters,"kSum(kPoly3,kRQISO)");
	double mean[128] = {1, 1, 1, 1};
	double std[128] = {5, 5, 5, 5};
	size_t nhp = 4;
	memcpy(parameters.kernel.hp_mean, mean, nhp * sizeof(double));
	memcpy(parameters.kernel.hp_std,std, nhp * sizeof(double));
	parameters.kernel.n_hp = nhp;
	break;
      }
    case 4:
      set_criteria(&parameters,"cLCB");
      parameters.crit_params[0] = 5;
      parameters.n_crit_params = 1;
      break;      
    case 5:
      set_criteria(&parameters,"cAopt");
      parameters.n_crit_params = 0;
      break;
    default:
      break;
    };

  bayesopt::Parameters parameters_class(parameters);
  boost::scoped_ptr<ExampleDiscrete> opt(new ExampleDiscrete(parameters_class));
  GLOBAL_MATPLOT.init(opt.get(),1);

  glutInit(&nargs, args);
  glutCreateWindow(50,50,800,650);
  //Call is redirected to DISPLAY() function call
  glutDisplayFunc( display );
  glutReshapeFunc( reshape );
  glutIdleFunc( idle );
  glutMotionFunc( motion );
  glutMouseFunc( mouse );
  glutPassiveMotionFunc(passive);    
  glutKeyboardFunc( keyboard );        
  glutMainLoop();    



}
class opt{

  std::mutex mtx;
  int shared_data;
  std::condition_variable cv;
  bool received;
  bool sent;
  
void sender(){
  
  int i = 0;
  while(i != 10)
  {
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    shared_data++;
    {
      std::unique_lock<std::mutex> lock1(mtx);
      
      cv.wait(lock1, [&]{ return this->received; });
      std::cout<<"Set data to " << shared_data<<std::endl;

      sent = true;
      received = false;
    }
    
    cv.notify_one();
   
    i++;

  }
    
}
  bool waitsend(){
    return sent;
}
void receiver()
{

  int i = 0;
  std::cout<<"Wait for data..."<<std::endl;
  while(true)
  {
    {
    std::unique_lock<std::mutex> lock(mtx);
    cv.wait(lock, [&]{ return this->sent; });
    std::cout<<"Received data "<< shared_data<<std::endl;
    i++;
    received = true;
    sent = false;
     if (shared_data >= 10) break; 
    lock.unlock();
    }
    cv.notify_one();
  }
    
}
void dostuff()
{
  bopt_params parameters = initialize_parameters_to_default();
  parameters.n_init_samples = 7;
  parameters.n_iterations = 100;
  parameters.verbose_level = 2;

  bayesopt::Parameters parameters_class(parameters);
  boost::scoped_ptr<ExampleDiscrete> opt(new ExampleDiscrete(parameters_class));
  GLOBAL_MATPLOT.init(opt.get(),1);
  int nargs;
  char *args[1];
  glutInit(&nargs, args);
  glutCreateWindow(50,50,800,650);
  //Call is redirected to DISPLAY() function call
  glutDisplayFunc( display );
  glutReshapeFunc( reshape );
  glutIdleFunc( idle );
  glutMotionFunc( motion );
  glutMouseFunc( mouse );
  glutPassiveMotionFunc(passive);    
  glutKeyboardFunc( keyboard );        
  glutMainLoop();    



}
public:
opt()
  {
    shared_data = 0;
    received = true;
    sent = false;
  }
void establish()
{
  std::thread t1(&opt::sender, this);
  std::thread t2(&opt::receiver, this);
  std::thread t3(&opt::dostuff, this);
  t1.join();
  t2.join();
  t3.join();



}
};

int main(int nargs, char *args[])
{
  
  opt opt_;
  opt_.establish();
  
  return 0;
}
