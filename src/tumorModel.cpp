#include "tumorModel.hpp"
#include <random>

GaussianNoise::GaussianNoise(double mean, double std):dist_(mean,std){}
double GaussianNoise::noise()
{
  return dist_(generator_);
}

PolarPolygon::PolarPolygon(size_t num_verices, double radius, double x_translation, double y_translation):
                num_vertices_(num_verices),radius_(radius), x_trans_(x_translation), y_trans_(y_translation)
{
  if(num_verices>2){
      computeVertices();
  }
  
}
void PolarPolygon::computeVertices()
{
  for(size_t i= 0; i<num_vertices_; i++)
  {
    vertices_x_.push_back(std::make_pair(cos((M_PI*2/num_vertices_)*i), cos( (M_PI*2/num_vertices_)*(i+1))));
    vertices_y_.push_back(std::make_pair( sin((M_PI*2/num_vertices_)*i), sin( (M_PI*2/num_vertices_)*(i+1))));
    radians_lookup_.push_back((M_PI*2/num_vertices_)*i);
  }
  radians_lookup_.push_back(M_PI*2);
}
double PolarPolygon::rTheta_(const std::pair<double,double> &x,const std::pair<double,double> &y,const double &theta )
{
  return ((x.second*y.first- x.first*y.second)/(x.second-x.first)) * 1/( sin(theta) - ((y.second-y.first)/(x.second-x.first))*cos(theta) );
}
double PolarPolygon::polar_x_(const double &theta, const double &rad )
{
  return  rad * cos(theta);
}
double PolarPolygon::polar_y_(const double &theta, const double &rad )
{
  return rad * sin(theta);
}
double PolarPolygon::r_polar_(const double &x,const double &y)
{
  return sqrt(y*y + x*x);
}
double PolarPolygon::theta_polar_(const double &x,const double &y)
{
    double angle = atan2(y, x);
    if (angle < 0 )
    {
      return 2*M_PI + angle;
    }
    else 
    {
      return angle;
    }
    


}
double PolarPolygon::insidePolygon(const double &x, const double & y){
 double r = r_polar_(x-x_trans_,y-y_trans_);
 double theta = theta_polar_(x-x_trans_,y-y_trans_);
 if(isnan(theta))
 {
  return 0.0;
 }
 //TODO Compare radius to line segment radius
 double rad = abs(evaluate(theta,radius_));
 return rad;
}
double PolarPolygon::outsidePolygon(const double &x, const double & y, double epsilon){
double r =  r_polar_(x-x_trans_,y-y_trans_);
 double theta =   theta_polar_(x-x_trans_,y-y_trans_);
 if(isnan(theta))
 {
  return 0.0;
 }
 //TODO Compare radius to line segment radius
 double rad = abs(evaluate(theta, radius_ + epsilon));
 return rad;
}
double PolarPolygon::insideCircle(const double &x, const double & y){
 double r = r_polar_(x-x_trans_,y-y_trans_);
 double theta = theta_polar_(x-x_trans_,y-y_trans_);
 if(isnan(theta))
 {
  return 0.0;
 }
 //TODO Compare radius to line segment radius
 double rad = abs(evaluateCircle(theta,radius_));
 return rad;
}
double PolarPolygon::outsideCircle(const double &x, const double & y, double epsilon){
double r =  r_polar_(x-x_trans_,y-y_trans_);
 double theta =   theta_polar_(x-x_trans_,y-y_trans_);
 if(isnan(theta))
 {
  return 0.0;
 }
 //TODO Compare radius to line segment radius
 double rad = abs(evaluateCircle(theta, radius_ + epsilon));
 return rad;
}

double PolarPolygon::evaluateCircle(const double & theta, double rad)
{

  return 1;

}
double PolarPolygon::evaluate(const double & theta, double rad)
{
  int idx = 0;
  for (size_t i = 0; i < radians_lookup_.size(); i++)
  {
      
      if (radians_lookup_[i] > theta)
      {
        idx = i; 
        break;
      }

  }
 double rt = rTheta_(std::make_pair( polar_x_(radians_lookup_[idx-1], rad) ,polar_x_(radians_lookup_[idx],rad)), std::make_pair( polar_y_(radians_lookup_[idx-1],rad), polar_y_(radians_lookup_[idx],rad) ), theta);
  return rt;

}
 Triangle::Triangle(bayesopt::Parameters par):
    ContinuousModel(2,par), triangle_(PolarPolygon(3, 0.25, 0.5,0.5)),  gaussianNoise_(0.0,0.005) {}

  double Triangle::smoothstep(double edge0, double edge1, double x) 
  {
    // Scale, bias and saturate x to 0..1 range
    x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0); 
    // Evaluate polynomial
    return x * x * (3 - 2 * x);
  } 
  double Triangle::clamp(double x, double lowerlimit, double upperlimit) {
    if (x < lowerlimit)
        x = lowerlimit;
    if (x > upperlimit)
        x = upperlimit;
    return x;
}

  double Triangle::evaluateSample( const vectord& xin)
  {
     if (xin.size() != 2)
      {
	std::cout << "WARNING: This only works for 2D inputs." << std::endl
		  << "WARNING: Using only first two components." << std::endl;
      }
    double x = xin(0);
    double y = xin(1);
    double epsilon = 0.3;
    double y_0 = 0.5;
    double x_0 = 0.5;
    double r = sqrt((x - x_0)*(x - x_0) + (y - y_0)*(y - y_0));

    double radius_inner = triangle_.insidePolygon(x,y);
    double radius_outer = triangle_.outsidePolygon(x,y, epsilon);
    if(r < radius_inner){
      return 0+ gaussianNoise_.noise();
    }
    else if (r >= radius_outer) 
    {
      return 1.0+ gaussianNoise_.noise();
    } else 
    {
      // smoothstep to interpolate between 1 and 0 for the edge of the circle
      return smoothstep(radius_inner,radius_outer, r);
    }

   
  }

  bool Triangle::checkReachability(const vectord &query)
  {
    return true;
    };


 Rectangle::Rectangle(bayesopt::Parameters par):
    ContinuousModel(2,par), rectangle_(PolarPolygon(4, 0.25, 0.5,0.5)),  gaussianNoise_(0.0,0.005) {}

  double Rectangle::smoothstep(double edge0, double edge1, double x) 
  {
    // Scale, bias and saturate x to 0..1 range
    x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0); 
    // Evaluate polynomial
    return x * x * (3 - 2 * x);
  } 
  double Rectangle::clamp(double x, double lowerlimit, double upperlimit) {
    if (x < lowerlimit)
        x = lowerlimit;
    if (x > upperlimit)
        x = upperlimit;
    return x;
}

  double Rectangle::evaluateSample( const vectord& xin)
  {
     if (xin.size() != 2)
      {
	std::cout << "WARNING: This only works for 2D inputs." << std::endl
		  << "WARNING: Using only first two components." << std::endl;
      }
    double x = xin(0);
    double y = xin(1);
    double epsilon = 0.3;
    double y_0 = 0.5;
    double x_0 = 0.5;
    double r = sqrt((x - x_0)*(x - x_0) + (y - y_0)*(y - y_0));

    double radius_inner = rectangle_.insidePolygon(x,y);
    double radius_outer = rectangle_.outsidePolygon(x,y, epsilon);
    if(r < radius_inner){
      return 0+ gaussianNoise_.noise();
    }
    else if (r >= radius_outer) 
    {
      return 1.0+ gaussianNoise_.noise();
    } else 
    {
      // smoothstep to interpolate between 1 and 0 for the edge of the circle
      return smoothstep(radius_inner,radius_outer, r);
    }

   
  }

  bool Rectangle::checkReachability(const vectord &query)
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
   
    if ((x - 0.5) * (x - 0.5)+ (y - 0.5) * (y - 0.5) <= 0.01)
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
    ContinuousModel(2,par), gaussianNoise_(0.0,0.005), circle_(PolarPolygon(1, 0.01,0.2,0.5)) {}


  double TwoCircles::smoothstep(double edge0, double edge1, double x) 
  {
    // Scale, bias and saturate x to 0..1 range
    x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0); 
    // Evaluate polynomial
    return x * x * (3 - 2 * x);
  } 
  double TwoCircles::clamp(double x, double lowerlimit, double upperlimit) {
    if (x < lowerlimit)
        x = lowerlimit;
    if (x > upperlimit)
        x = upperlimit;
    return x;
}  
  double TwoCircles::evaluateSample( const vectord& xin){
    if (xin.size() != 2)
      {
	std::cout << "WARNING: This only works for 2D inputs." << std::endl
		  << "WARNING: Using only first two components." << std::endl;
      }
    double x = xin(0);
    double y = xin(1);

    double y_0 = 0.5;
    double x_0 = 0.2;
    double y_1 = 0.5;
    double x_1 = 0.8;

    double r0 = sqrt((x - x_0)*(x - x_0) + (y - y_0)*(y - y_0));
    double r1 = sqrt((x - x_1)*(x - x_1) + (y - y_1)*(y - y_1));
    double radius_0 = 0.01;
    double radius_1 = 0.03;

    double epsilon = 0.1;

    double radius_inner = circle_.insidePolygon(x,y);
    double radius_outer = circle_.outsidePolygon(x,y, epsilon);

    if(r0 < radius_inner){
      return 0+ gaussianNoise_.noise();
    }
    else if (r0 >= radius_outer) 
    {
      return 1.0+ gaussianNoise_.noise();
    } else 
    {
      // smoothstep to interpolate between 1 and 0 for the edge of the circle
      return smoothstep(radius_inner,radius_outer, r0);
    }

    /*
    if ((x - 0.2)*(x - 0.2)+ (y-0.5)*(y-0.5) <= 0.01)

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
    */
  
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
    double r = sqrt((x - x_0)*(x - x_0) + (y - y_0)*(y - y_0));
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

