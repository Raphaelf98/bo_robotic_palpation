#include "tumorModel.hpp"

 Triangle::Triangle(bayesopt::Parameters par):
    ContinuousModel(2,par) {}



  double Triangle::evaluateSample( const vectord& xin)
  {
     if (xin.size() != 2)
      {
	std::cout << "WARNING: This only works for 2D inputs." << std::endl
		  << "WARNING: Using only first two components." << std::endl;
      }
    double x = xin(0);
    double y = xin(1);
    
    
      if (y >= 0.5 && y <= 1.5*x && y <= -1.5*x +1.5)
    {
      return 0.0;
    }
    if (y >= 0.4 && y <= 1.5*x+0.2 && y <= -1.5*x +1.7)
    {
      return 5.0;
    }
    else{
      return 10.0;
    }
  }

  bool Triangle::checkReachability(const vectord &query)
  {
    return true;
    };
    
  Circle::Circle(bayesopt::Parameters par):
    ContinuousModel(2,par) {}

  double Circle::evaluateSample( const vectord& xin)
  {
     if (xin.size() != 2)
      {
	std::cout << "WARNING: This only works for 2D inputs." << std::endl
		  << "WARNING: Using only first two components." << std::endl;
      }
    double x = xin(0);
    double y = xin(1);

    if ((x-0.5)*(x-0.5)+ (y-0.5)*(y-0.5) <= 0.01)

    {
      return 1.0;
    }
    
    else
    {
      return 2.0;
    }
  
   //return (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5);
    
  }

  bool Circle::checkReachability(const vectord &query)
  {return true;};


TwoCircles::TwoCircles(bayesopt::Parameters par):
    ContinuousModel(2,par) {}

  
  double TwoCircles::evaluateSample( const vectord& xin){
    if (xin.size() != 2)
      {
	std::cout << "WARNING: This only works for 2D inputs." << std::endl
		  << "WARNING: Using only first two components." << std::endl;
      }
    double x = xin(0);
    double y = xin(1);

    if ((x-0.2)*(x-0.2)+ (y-0.5)*(y-0.5) <= 0.01)

    {
      return 1.0;
    }
    if ((x-0.8)*(x-0.8)+ (y-0.5)*(y-0.5) <= 0.03)

    {
      return 1.0;
    }
    
    else
    {
      return 2.0;
    }
  
  }
  

    
  bool TwoCircles::checkReachability(const vectord &query)
  {
    return true;
  };




