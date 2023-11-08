#include "tumorModel.hpp"
#include <random>

GaussianNoise::GaussianNoise(double mean, double std):dist_(mean,std){}
double GaussianNoise::noise()
{
  return dist_(generator_);
}


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




SmoothCircle::SmoothCircle(bayesopt::Parameters par):
    ContinuousModel(2,par), gaussianNoise_(0.0,0.005){}

  double SmoothCircle::smoothstep(double edge0, double edge1, double x) 
  {
    // Scale, bias and saturate x to 0..1 range
    x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0); 
    // Evaluate polynomial
    return x * x * (3 - 2 * x);
  } 
  double SmoothCircle::clamp(double x, double lowerlimit, double upperlimit) {
    if (x < lowerlimit)
        x = lowerlimit;
    if (x > upperlimit)
        x = upperlimit;
    return x;
}

  double SmoothCircle::evaluateSample( const vectord& xin)
  {
    if (xin.size() != 2)
    {
	    std::cout << "WARNING: This only works for 2D inputs." << std::endl;
    }
    double x = xin(0);
    double y = xin(1);
    double y_0 = 0.5;
    double x_0 = 0.5;
    double r = sqrt((x - x_0) * (x - x_0) + (y - y_0) * (y - y_0));
    double radius = 0.1;
    double epsilon = 0.1;

    if (r <= radius) 
    {
      return 0.0 + gaussianNoise_.noise();
    } else if (r >= radius + epsilon) 
    {
      return 1.0 + gaussianNoise_.noise();
    } else 
    {
      // smoothstep to interpolate between 1 and 0 for the edge of the circle
      return smoothstep(radius, radius + epsilon, r);
    }

    
  }
  
  
  bool SmoothCircle::checkReachability(const vectord &query)
  {return true;};

