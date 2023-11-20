#include "tumorModel.hpp"
#include <random>

GaussianNoise::GaussianNoise(double mean, double std):dist_(mean,std){}
double GaussianNoise::noise()
{
  return dist_(generator_);
}

PolarPolygon::PolarPolygon(size_t num_verices, double radius, double x_translation, double y_translation, double epsilon):
                num_vertices_(num_verices),radius_(radius), x_trans_(x_translation), y_trans_(y_translation), epsilon_(epsilon)
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

 double rad = abs(evaluate(theta, radius_ + epsilon));
 return rad;
}
double PolarPolygon::insideCircle(const double &x, const double & y){

 
 return radius_;
}
double PolarPolygon::outsideCircle(const double &x, const double & y, double epsilon){

 
 return radius_ + epsilon;
}


double PolarPolygon::evaluate(const double & theta, double rad)
{
  int idx = 0;
  //deals with undefined behaviour for theta = 0
  if(theta == 0){
    return rad;
  }

  for (size_t i = 0; i < radians_lookup_.size(); i++)
  {
      
      if (radians_lookup_[i] >= theta)
      {
        idx = i; 
        break;
      }

  }
  
 double rt = rTheta_(std::make_pair( polar_x_(radians_lookup_[idx-1], rad) ,polar_x_(radians_lookup_[idx],rad)), std::make_pair( polar_y_(radians_lookup_[idx-1],rad), polar_y_(radians_lookup_[idx],rad) ), theta);
  return rt;

}
//s e [0,1] (parametric funtion variable)
double  PolarPolygon::f_x_(const double& s)
{
 return x_trans_ + polar_x_(s*2*M_PI, evaluate(s*2*M_PI, radius_+0.5*epsilon_) );
}
double  PolarPolygon::f_y_(const double& s)
{
return y_trans_ + polar_y_(s*2*M_PI, evaluate(s*2*M_PI, radius_+0.5*epsilon_) );
}
std::function<double(const double&)> PolarPolygon::fParametric_x()
{
  return std::bind(&PolarPolygon::f_x_,this, std::placeholders::_1);
}

std::function<double(const double&)> PolarPolygon::fParametric_y()
{
  return std::bind(&PolarPolygon::f_y_,this, std::placeholders::_1);
}

Shape::Shape(bayesopt::Parameters par):ContinuousModel(2,par){

}

double Shape::smoothstep(const double x_l, const double x_r, const double x, const double y_l, const double y_r )
  {
    double a = (-2*y_l + 2*y_r)/(x_l*x_l*x_l - 3*x_l*x_l*x_r + 3*x_l*x_r*x_r - x_r*x_r*x_r);
    double b = (3*x_l*y_l - 3*x_l*y_r + 3*x_r*y_l - 3*x_r*y_r)/(x_l*x_l*x_l - 3*x_l*x_l*x_r + 3*x_l*x_r*x_r - x_r*x_r*x_r);
    double c = (-6*x_l*x_r*y_l + 6*x_l*x_r*y_r)/(x_l*x_l*x_l - 3*x_l*x_l*x_r + 3*x_l*x_r*x_r - x_r*x_r*x_r);
    double d = (x_l*x_l*x_l*y_r - 3*x_l*x_l*x_r*y_r + 3*x_l*x_r*x_r*y_l - x_r*x_r*x_r*y_l)/(x_l*x_l*x_l - 3*x_l*x_l*x_r + 3*x_l*x_r*x_r - x_r*x_r*x_r);
    return a*x*x*x + b*x*x + c*x + d;
  }
  double Shape::clamp(double x, double lowerlimit, double upperlimit) {
    if (x < lowerlimit)
        x = lowerlimit;
    if (x > upperlimit)
        x = upperlimit;
    return x;
}



 Triangle::Triangle(bayesopt::Parameters par):
   Shape(par), triangle_(PolarPolygon(3, 0.2, 0.5,0.5,0.1)),  gaussianNoise_(0.0,0.005) {}

  std::function<double (const double&)> Triangle::f_x()
  {
    return triangle_.fParametric_x();
  }
  std::function<double (const double&)> Triangle::f_y()
  {
    return triangle_.fParametric_y();
  }


  double Triangle::evaluateSample( const vectord& xin)
  {
     if (xin.size() != 2)
      {
	    std::cout << "WARNING: This only works for 2D inputs." << std::endl
		  << "WARNING: Using only first two components." << std::endl;
      }
      double low= 0.0;
      double high = 1.0;
      double x = xin(0);
      double y = xin(1);
      double epsilon = 0.1;
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
        return smoothstep(radius_inner,radius_outer, r,low,high);
      }

   
  }

  bool Triangle::checkReachability(const vectord &query)
  {
    return true;
    };


 Rectangle::Rectangle(bayesopt::Parameters par):
   Shape(par), rectangle_(PolarPolygon(4, 0.15, 0.5,0.5,0.1)),  gaussianNoise_(0.0,0.00) {}

  

  double Rectangle::evaluateSample( const vectord& xin)
  {
     if (xin.size() != 2)
      {
	std::cout << "WARNING: This only works for 2D inputs." << std::endl
		  << "WARNING: Using only first two components." << std::endl;
      }
      double low= 1.0;
      double high = 2.0;
    double x = xin(0);
    double y = xin(1);
    double epsilon = 0.1;
    double y_0 = 0.5;
    double x_0 = 0.5;
    double r = sqrt((x - x_0)*(x - x_0) + (y - y_0)*(y - y_0));

    double radius_inner = rectangle_.insidePolygon(x,y);
    double radius_outer = rectangle_.outsidePolygon(x,y, epsilon);
    if(r < radius_inner){
      return 1 + gaussianNoise_.noise();
    }
    else if (r >= radius_outer) 
    {
      return 2 + gaussianNoise_.noise();
    } else 
    {
      // smoothstep to interpolate between 1 and 0 for the edge of the circle
      return smoothstep(radius_inner,radius_outer, r,low,high);
    }

   
  }
 std::function<double (const double&)> Rectangle::f_x()
  {
    return rectangle_.fParametric_x();
  }
  std::function<double (const double&)> Rectangle::f_y()
  {
    return rectangle_.fParametric_y();
  }


  bool Rectangle::checkReachability(const vectord &query)
  {
    return true;
    };
    

    /*
  Circle::Circle(bayesopt::Parameters par):
   Shape(par) {}

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
 std::function<double (const double&)> Circle::f_x()
  {
    return triangle_.fParametric_x();
  }
  std::function<double (const double&)> Circle::f_y()
  {
    return triangle_.fParametric_y();
  }

  bool Circle::checkReachability(const vectord &query)
  {return true;};

*/
TwoCircles::TwoCircles(bayesopt::Parameters par,double r_1, double r_2,double x_t_1,double x_t_2,
                        double y_t_1,double y_t_2,double epsilon):
                        Shape(par),
                       gaussianNoise_(0.0,0.005), 
                      r_1_(r_1), r_2_(r_2), x_t_1_(x_t_1), x_t_2_(x_t_2), y_t_1_(y_t_1), y_t_2_(y_t_2),epsilon_(epsilon),
                      circle_1_(PolarPolygon(1, r_1_,x_t_1_,y_t_1_,0.1)),circle_2_(PolarPolygon(1, r_2_,x_t_2_,y_t_2_, 0.1)) {}


  
  double TwoCircles::evaluateSample( const vectord& xin){
    if (xin.size() != 2)
      {
	std::cout << "WARNING: This only works for 2D inputs." << std::endl
		  << "WARNING: Using only first two components." << std::endl;
      }
    double x = xin(0);
    double y = xin(1);
     double low= 0.0;
    double high = 1.0;
    double r1 = sqrt((x - x_t_1_)*(x - x_t_1_) + (y - y_t_1_)*(y - y_t_1_));
    double r2 = sqrt((x - x_t_2_)*(x - x_t_2_) + (y - y_t_2_)*(y - y_t_2_));

    double radius_inner_1 = circle_1_.insideCircle(x,y);
    double radius_outer_1 = circle_1_.outsideCircle(x,y, epsilon_);
    double radius_inner_2 = circle_2_.insideCircle(x,y);
    double radius_outer_2 = circle_2_.outsideCircle(x,y, epsilon_);

    if(r1 < r2)
    {
    if(r1 < radius_inner_1 )
    {
      return 0 + gaussianNoise_.noise();
    }
    if (r1 >= radius_outer_1 ){
      return 1.0 + gaussianNoise_.noise();
    } 
    if (r1 > radius_inner_1  && r1 < radius_outer_1)
    {
      return smoothstep(radius_inner_1,radius_outer_1, r1,low,high);
    }
   
    }
    else
    {
      if(r2 < radius_inner_2)
    {
      return 0 + gaussianNoise_.noise();
    }
    if (r2 >= radius_outer_2 )
    {
      return 1.0 + gaussianNoise_.noise();
    } 
    if (r2 > radius_inner_2  && r2 < radius_outer_2)
    {
        return smoothstep(radius_inner_2,radius_outer_2, r2,low,high);
    }
    }
  }
  std::function<double (const double&)> TwoCircles::f_x()
  {
    return circle_1_.fParametric_x();
  }
  std::function<double (const double&)> TwoCircles::f_y()
  {
    return circle_1_.fParametric_y();
  }
    
  bool TwoCircles::checkReachability(const vectord &query)
  {
    return true;
  };

SmoothCircle::SmoothCircle(bayesopt::Parameters par):
    Shape(par), gaussianNoise_(0.0,0.005){}

  

  double SmoothCircle::evaluateSample( const vectord& xin)
  {
    if (xin.size() != 2)
    {
	    std::cout << "WARNING: This only works for 2D inputs." << std::endl;
    }
    double low= 0.0;
    double high = 1.0;
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
      return smoothstep(radius, radius + epsilon, r, low, high);
    }

    
  }
std::function<double (const double&)> SmoothCircle::f_x()
  {
    return circle_.fParametric_x();
  }
  std::function<double (const double&)> SmoothCircle::f_y()
  {
    return circle_.fParametric_y();
  }
  bool SmoothCircle::checkReachability(const vectord &query)
  {return true;};

