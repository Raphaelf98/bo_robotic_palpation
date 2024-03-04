#include <random>

#include "tumorModel.hpp"
#include"parameters.hpp"
#include"fileparser.hpp"

/**
 * @file
 * 
 */
GaussianNoise::GaussianNoise(double mean, double std):dist_(mean,std){}
double GaussianNoise::noise()
{
  return dist_(generator_);
}

PolarPolygon::PolarPolygon(size_t num_vertices, double radius, double x_translation, double y_translation):
                num_vertices_(num_vertices),radius_(radius), x_trans_(x_translation), y_trans_(y_translation)
{
  if(num_vertices>2)
  {
      computeVertices_();
      circle_ = false;
  }
  else
  {
    std::cout<<"Shape is circle"<<std::endl;
    circle_ = true;
  }
  
}
void PolarPolygon::computeVertices_()
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
double PolarPolygon::polygonRadius(const double &x, const double & y)
{
  if(!circle_)
  {
    double theta = theta_polar_(x-x_trans_,y-y_trans_);
    if(isnan(theta))
    {
     return 0.0;
    }

    return abs(evaluate_(theta,radius_));
  }
  else return radius_;
}
double PolarPolygon::evaluate_(const double & theta, const double & rad)
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
//s e [0,1] (parametric function variable)
double  PolarPolygon::f_x_(const double& s)
{ 
  //make sure to return different values for polygon and circle
  if(circle_)
  {
    return x_trans_ + polar_x_(s*2*M_PI, radius_);
  }
  else
  {
    return x_trans_ + polar_x_(s*2*M_PI, evaluate_(s*2*M_PI, radius_) );
  }
 
}
double  PolarPolygon::f_y_(const double& s)
{
  if(circle_)
  {
    return y_trans_ + polar_y_(s*2*M_PI, radius_);
  }
  else
  {
    return y_trans_ + polar_y_(s*2*M_PI, evaluate_(s*2*M_PI, radius_) );
  }
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


 Triangle::Triangle(bayesopt::Parameters par, double low,double high, double radius, double x_trans, double y_trans, double epsilon, double noise):
   Shape(par),  gaussianNoise_(0.0,noise) 
   {
      low_stiffness_ = low;
      high_stiffness_= high;
      radius_ = radius;
      x_trans_= x_trans;
      y_trans_ = y_trans;
      epsilon_ = epsilon;
      inner_triangle_ = PolarPolygon(3,radius_,x_trans_,y_trans_);
      outer_triangle_ = PolarPolygon(3,radius_+epsilon,x_trans_,y_trans_);
      groundTruth_triangle_ = PolarPolygon(3,radius_+epsilon*0.5,x_trans_,y_trans_);
   }
   Triangle::Triangle(bayesopt::Parameters par, TumorModelParameters &model_parameters):
   Shape(par),  gaussianNoise_(0.0,model_parameters.triangle_noise) 
   {
      low_stiffness_ = model_parameters.triangle_low;
      high_stiffness_= model_parameters.triangle_high;
      radius_ = model_parameters.triangle_radius;
      x_trans_= model_parameters.triangle_x_trans;
      y_trans_ = model_parameters.triangle_y_trans;
      epsilon_ = model_parameters.triangle_epsilon;

      inner_triangle_ = PolarPolygon(3,radius_,x_trans_,y_trans_);
      outer_triangle_ = PolarPolygon(3,radius_+epsilon_,x_trans_,y_trans_);
      groundTruth_triangle_ = PolarPolygon(3,radius_+epsilon_*0.5,x_trans_,y_trans_);
   }

  std::function<double (const double&)> Triangle::f_x()
  {
    return groundTruth_triangle_.fParametric_x();
  }
  std::function<double (const double&)> Triangle::f_y()
  {
    return groundTruth_triangle_.fParametric_y();
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
      double r = sqrt((x - x_trans_)*(x - x_trans_) + (y - y_trans_)*(y - y_trans_));

      double radius_inner = inner_triangle_.polygonRadius(x,y);
      double radius_outer = outer_triangle_.polygonRadius(x,y);

      if(r < radius_inner){
        return -high_stiffness_+ gaussianNoise_.noise();
      }
      else if (r >= radius_outer) 
      {
        return -low_stiffness_+ gaussianNoise_.noise();
      } else 
      {
        // smoothstep to interpolate between 1 and 0 for the edge of the circle
        return smoothstep(radius_inner, radius_outer, r, -high_stiffness_,-low_stiffness_);
      }
  }

  bool Triangle::checkReachability(const vectord &query)
  {
    return true;
    };

void Triangle::saveGroundTruth(const size_t c_points, std::string file_path)
{
      std::vector<double> cX=linSpace(0,1,c_points);
      std::vector<double> cY=linSpace(0,1,c_points);
      std::vector<std::vector<double>> cZ(c_points, std::vector<double>(c_points));
      vectord q_(2);
      q_(0) = cX[0]; 
      q_(1) = cY[0];
      

      for(int i=0;i<c_points;++i)
	    {  
	        for(int j=0;j<c_points;++j)
	        {
	            vectord q(2);
	            q(0) = cX[j]; 
              q(1) = cY[i];
	            cZ[i][j]= evaluateSample(q);
	        }
          
	    }
      std::cout<<"Write gound truth heatmap to file: " <<file_path<<std::endl;
      std::string gt = FILE_GROUND_TRUTH_HEATMAP;
      file_path = file_path+gt;
      saveFileToCSV(file_path,cZ );
} 
 Rectangle::Rectangle(bayesopt::Parameters par,double low,double high, double radius, double x_trans, double y_trans, double epsilon, double noise):
   Shape(par), gaussianNoise_(0.0,noise) 
   {
      low_stiffness_ = low;
      high_stiffness_= high;
      radius_ = radius;
      x_trans_= x_trans;
      y_trans_ = y_trans;
      epsilon_ = epsilon;
      inner_rectangle_ = PolarPolygon(4,radius_,x_trans_,y_trans_);
      outer_rectangle_ = PolarPolygon(4,radius_+epsilon,x_trans_,y_trans_);
      groundTruth_rectangle_ =PolarPolygon(4,radius_+epsilon*0.5,x_trans_,y_trans_);

   }
  Rectangle::Rectangle(bayesopt::Parameters par, TumorModelParameters &model_parameters):
   Shape(par),  gaussianNoise_(0.0,model_parameters.rectangle_noise) 
   {
      low_stiffness_ = model_parameters.rectangle_low;
      high_stiffness_= model_parameters.rectangle_high;
      radius_ = model_parameters.rectangle_radius;
      x_trans_= model_parameters.rectangle_x_trans;
      y_trans_ = model_parameters.rectangle_y_trans;
      epsilon_ = model_parameters.rectangle_epsilon;

      inner_rectangle_ = PolarPolygon(4,radius_,x_trans_,y_trans_);
      outer_rectangle_ = PolarPolygon(4,radius_+epsilon_,x_trans_,y_trans_);
      groundTruth_rectangle_ = PolarPolygon(4,radius_+epsilon_*0.5,x_trans_,y_trans_);
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
    
    double r = sqrt((x - x_trans_)*(x - x_trans_) + (y - y_trans_)*(y - y_trans_));

    double radius_inner = inner_rectangle_.polygonRadius(x,y);
    double radius_outer = outer_rectangle_.polygonRadius(x,y);
    if(r < radius_inner)
    {
      return -high_stiffness_ + gaussianNoise_.noise();
    }
    else if (r >= radius_outer) 
    {
      return -low_stiffness_ + gaussianNoise_.noise();
    } 
    else 
    {
      // smoothstep to interpolate between 1 and 0 for the edge of the circle
      return smoothstep(radius_inner,radius_outer, r,-high_stiffness_,-low_stiffness_) + gaussianNoise_.noise();
    }

   
  }
 std::function<double (const double&)> Rectangle::f_x()
  {
    return groundTruth_rectangle_.fParametric_x();
  }
  std::function<double (const double&)> Rectangle::f_y()
  {
    return groundTruth_rectangle_.fParametric_y();
  }


  bool Rectangle::checkReachability(const vectord &query)
  {
    return true;
    };
void Rectangle::saveGroundTruth(const size_t c_points, std::string file_path)
{
      std::vector<double> cX=linSpace(0,1,c_points);
      std::vector<double> cY=linSpace(0,1,c_points);
      std::vector<std::vector<double>> cZ(c_points, std::vector<double>(c_points));
      vectord q_(2);
      q_(0) = cX[0]; 
      q_(1) = cY[0];
      

      for(int i=0;i<c_points;++i)
	    {  
	        for(int j=0;j<c_points;++j)
	        {
	            vectord q(2);
	            q(0) = cX[j]; 
              q(1) = cY[i];
	            cZ[i][j]= evaluateSample(q);
	        }
          
	    }
      std::cout<<"Write gound truth heatmap to file: " <<file_path<<std::endl;
      std::string gt = FILE_GROUND_TRUTH_HEATMAP;
      file_path = file_path+gt;
      saveFileToCSV(file_path,cZ );
}     

   
TwoCircles::TwoCircles(bayesopt::Parameters par,double low_stiffness, double high_stiffness ,double r_1, double r_2,double x_t_1,double x_t_2,
                        double y_t_1,double y_t_2,double epsilon,double noise):
                        Shape(par), circle_count_(1),
                       gaussianNoise_(0.0,noise), 
                      r_1_(r_1), r_2_(r_2), x_t_1_(x_t_1), x_t_2_(x_t_2), y_t_1_(y_t_1), y_t_2_(y_t_2),epsilon_(epsilon),
                      inner_circle_1_(PolarPolygon(1, r_1_,x_t_1_,y_t_1_)),inner_circle_2_(PolarPolygon(1, r_2_,x_t_2_,y_t_2_)),
                      outer_circle_1_(PolarPolygon(1, r_1_+epsilon,x_t_1_,y_t_1_)),outer_circle_2_(PolarPolygon(1, r_2_+epsilon,x_t_2_,y_t_2_)),
                      groundTruth_circle_1_(PolarPolygon(1, r_1_+0.5*epsilon,x_t_1_,y_t_1_)),groundTruth_circle_2_(PolarPolygon(1, r_2_+0.5*epsilon,x_t_2_,y_t_2_))
                      {

                          low_stiffness_ = low_stiffness;
                          high_stiffness_ = high_stiffness;
                      }

TwoCircles::TwoCircles(bayesopt::Parameters par, TumorModelParameters &model_parameters):
   Shape(par),  
   gaussianNoise_(0.0,model_parameters.two_circles_noise) , 
   circle_count_(1)
   {
      low_stiffness_ = model_parameters.two_circles_low;
      high_stiffness_= model_parameters.two_circles_high;
      epsilon_ = model_parameters.two_circles_epsilon;

      r_1_ = model_parameters.two_circles_radius_1;
      x_t_1_= model_parameters.two_circles_x_trans_1;
      y_t_1_ = model_parameters.two_circles_y_trans_1;
      
      r_2_ = model_parameters.two_circles_radius_2;
      x_t_2_= model_parameters.two_circles_x_trans_2;
      y_t_2_ = model_parameters.two_circles_y_trans_2;

      inner_circle_1_ = PolarPolygon(1,r_1_,x_t_1_,y_t_1_);
      outer_circle_1_ = PolarPolygon(1,r_1_+epsilon_,x_t_1_,y_t_1_);
      groundTruth_circle_1_ = PolarPolygon(1,r_1_+epsilon_*0.5,x_t_1_,y_t_1_);

      inner_circle_2_ = PolarPolygon(1,r_2_,x_t_2_,y_t_2_);
      outer_circle_2_ = PolarPolygon(1,r_2_+epsilon_,x_t_2_,y_t_2_);
      groundTruth_circle_2_ = PolarPolygon(1,r_2_+epsilon_*0.5,x_t_2_,y_t_2_);
   }
  
  
  double TwoCircles::evaluateSample( const vectord& xin)
  {
    if (xin.size() != 2)
      {
	std::cout << "WARNING: This only works for 2D inputs." << std::endl
		  << "WARNING: Using only first two components." << std::endl;
      }
    double x = xin(0);
    double y = xin(1);
     
    double r1 = sqrt((x - x_t_1_)*(x - x_t_1_) + (y - y_t_1_)*(y - y_t_1_));
    double r2 = sqrt((x - x_t_2_)*(x - x_t_2_) + (y - y_t_2_)*(y - y_t_2_));

    double radius_inner_1 = inner_circle_1_.polygonRadius(x,y);
    double radius_outer_1 = outer_circle_1_.polygonRadius(x,y);

    double radius_inner_2 = inner_circle_2_.polygonRadius(x,y);
    double radius_outer_2 = outer_circle_2_.polygonRadius(x,y);

    if(r1 < r2)
    {
    if(r1 < radius_inner_1 )
    {
      return -high_stiffness_ + gaussianNoise_.noise();
    }
    if (r1 >= radius_outer_1 ){
      return -low_stiffness_ + gaussianNoise_.noise();
    } 
    if (r1 > radius_inner_1  && r1 < radius_outer_1)
    {
      return smoothstep(radius_inner_1,radius_outer_1, r1,-high_stiffness_,-low_stiffness_)+ gaussianNoise_.noise();
    }
    else 0.0;
   
    }
    else
    {
      if(r2 < radius_inner_2)
    {
      return -high_stiffness_ + gaussianNoise_.noise();
    }
    if (r2 >= radius_outer_2 )
    {
      return -low_stiffness_ + gaussianNoise_.noise();
    } 
    if (r2 > radius_inner_2  && r2 < radius_outer_2)
    {
      return smoothstep(radius_inner_2,radius_outer_2, r2,-high_stiffness_,-low_stiffness_)+ gaussianNoise_.noise();
    }
    else 0.0;
    }
  }
/*
Returns function pointer 
*/
std::function<double (const double&)> TwoCircles::f_x() 
  {
    if(circle_count_ == 1)
    {
      circle_count_++;
      return groundTruth_circle_1_.fParametric_x();
      

    } 
    if(circle_count_ == 2)
    {
      return groundTruth_circle_2_.fParametric_x();
    } 
    else nullptr;
  }
std::function<double (const double&)> TwoCircles::f_y() 
  {
    if(circle_count_ == 1)
    {
      
      return groundTruth_circle_1_.fParametric_y();
      
    } 
    if(circle_count_ == 2)
    {
      return groundTruth_circle_2_.fParametric_y();
    } 
    else nullptr;
  }
    
  bool TwoCircles::checkReachability(const vectord &query)
  {
    return true;
  };
void TwoCircles::saveGroundTruth(const size_t c_points, std::string file_path)
{
      std::vector<double> cX=linSpace(0,1,c_points);
      std::vector<double> cY=linSpace(0,1,c_points);
      std::vector<std::vector<double>> cZ(c_points, std::vector<double>(c_points));
      vectord q_(2);
      q_(0) = cX[0]; 
      q_(1) = cY[0];
      

      for(int i=0;i<c_points;++i)
	    {  
	        for(int j=0;j<c_points;++j)
	        {
	            vectord q(2);
	            q(0) = cX[j]; 
              q(1) = cY[i];
	            cZ[i][j]= evaluateSample(q);
	        }
          
	    }
      std::cout<<"Write gound truth heatmap to file: " <<file_path<<std::endl;
      std::string gt = FILE_GROUND_TRUTH_HEATMAP;
      file_path = file_path+gt;
      saveFileToCSV(file_path,cZ );
}     

Circle::Circle(bayesopt::Parameters par,double low,double high, double radius, double x_trans, double y_trans, double epsilon,double noise):
    Shape(par),
    gaussianNoise_(0.0,noise)
    {
      low_stiffness_ = low;
      high_stiffness_= high;
      radius_ = radius;
      x_trans_= x_trans;
      y_trans_ = y_trans;
      epsilon_ = epsilon;
      
      groundTruth_circle_ = PolarPolygon(1,radius_+epsilon*0.5,x_trans_,y_trans_);


    }
Circle::Circle(bayesopt::Parameters par, TumorModelParameters &model_parameters):
   Shape(par),  gaussianNoise_(0.0,model_parameters.circle_noise) 
   {
      low_stiffness_ = model_parameters.circle_low;
      high_stiffness_= model_parameters.circle_high;
      radius_ = model_parameters.circle_radius;
      x_trans_= model_parameters.circle_x_trans;
      y_trans_ = model_parameters.circle_y_trans;
      epsilon_ = model_parameters.circle_epsilon;

      
      groundTruth_circle_ = PolarPolygon(1,radius_+epsilon_*0.5,x_trans_,y_trans_);
   }
  
  double Circle::evaluateSample( const vectord& xin)
  {
    if (xin.size() != 2)
    {
	    std::cout << "WARNING: This only works for 2D inputs." << std::endl;
    }
    
    double x = xin(0);
    double y = xin(1);
   
    double r = sqrt((x - x_trans_)*(x - x_trans_) + (y - y_trans_)*(y - y_trans_));
    
    if (r <= radius_) 
    {
      return -(high_stiffness_ + gaussianNoise_.noise());

    } 
    else if (r >= radius_ + epsilon_) 
    {
      return -(low_stiffness_ + gaussianNoise_.noise());
    } 
    else 
    {
      // smoothstep to interpolate between 1 and 0 for the edge of the circle
      return smoothstep(radius_, radius_ + epsilon_, r, -high_stiffness_, -low_stiffness_)+ gaussianNoise_.noise();
    }
    
  }
std::function<double (const double&)> Circle::f_x()
  {
    return groundTruth_circle_.fParametric_x();
  }
  std::function<double (const double&)> Circle::f_y()
  {
    return groundTruth_circle_.fParametric_y();
  }
  bool Circle::checkReachability(const vectord &query)
  {return true;};

void Circle::saveGroundTruth(const size_t c_points, std::string file_path)
{
      std::vector<double> cX=linSpace(0,1,c_points);
      std::vector<double> cY=linSpace(0,1,c_points);
      std::vector<std::vector<double>> cZ(c_points, std::vector<double>(c_points));
      vectord q_(2);
      q_(0) = cX[0]; 
      q_(1) = cY[0];
      

      for(int i=0;i<c_points;++i)
	    {  
	        for(int j=0;j<c_points;++j)
	        {
	            vectord q(2);
	            q(0) = cX[j]; 
              q(1) = cY[i];
	            cZ[i][j]= evaluateSample(q);
	        }
          
	    }
      std::cout<<"Write gound truth heatmap to file: " <<file_path<<std::endl;
      std::string gt = FILE_GROUND_TRUTH_HEATMAP;
      file_path = file_path+gt;
      saveFileToCSV(file_path,cZ );
}     
