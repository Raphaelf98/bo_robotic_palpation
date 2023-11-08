#include "contour.hpp"
#include "param_loader.hpp"
#include<bayesopt/parameters.hpp>
using namespace alglib;

int main(int argc, char **argv)
{   
    bayesopt::Parameters par;
 
  
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
    //boost::scoped_ptr<TwoCircles> twoCircles(new TwoCircles(par));
    //Contour contour(twoCircles.get(),100);
    //boost::scoped_ptr<Circle> circle(new Circle(par));
    //Contour contour(circle.get(),100);
    boost::scoped_ptr<SmoothCircle> smoothCircle(new SmoothCircle(par));
    Contour contour(smoothCircle.get(),100);
    //run bayesian optimization
    contour.runGaussianProcess();
    contour.computeCluster();
    //dummy data
    contour.labelData();
    //contour.exploreContour();
    //contour.approximateContour();
    return 0;
}

