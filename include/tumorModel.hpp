#ifndef _TUMORMODEL_HPP_
#define _TUMORMODEL_HPP_
#include "bayesopt/bayesopt.hpp"
#include <random>
class GaussianNoise{
  public:
    GaussianNoise(double mean, double std);
    double noise();
  private:
    std::default_random_engine generator_;
    std::normal_distribution<double> dist_;

};
class PolarPolygon{
    public:
    PolarPolygon(){};
    PolarPolygon(size_t num_verices, double radius, double x_translation, double y_translation);
    double evaluate(const double & theta, double rad);
    double insidePolygon(const double &x, const double & y);
    double outsidePolygon(const double &x, const double & y, double epsilon);
    private:
    double radius_, x_trans_, y_trans_;
    size_t num_vertices_;
    std::vector<std::pair<double,double>> vertices_x_;
    std::vector<std::pair<double,double>> vertices_y_;
    std::vector<double> radians_lookup_;
    double polar_x_(const double &theta, const double &rad);
    double polar_y_(const double &theta, const double &rad);
    double r_polar_(const double &x,const double &y);
    double theta_polar_(const double &x,const double &y);
    void computeVertices();
    
    double rTheta_(const std::pair<double,double> &x,const std::pair<double,double> &y,const double &thetha );
    
};
class Triangle: public bayesopt::ContinuousModel
  {
    GaussianNoise gaussianNoise_;
    PolarPolygon triangle_;
  public:

    Triangle(bayesopt::Parameters par);
  
    double evaluateSample( const vectord& xin);
    double smoothstep(double edge0, double edge1, double x);
    double clamp(double x, double lowerlimit, double upperlimit);
    
    
  bool checkReachability(const vectord &query);

};
class Rectangle: public bayesopt::ContinuousModel
  {
    GaussianNoise gaussianNoise_;
    PolarPolygon rectangle_;
  public:
    Rectangle(bayesopt::Parameters par);
  
    double evaluateSample( const vectord& xin);
    double smoothstep(double edge0, double edge1, double x);
    double clamp(double x, double lowerlimit, double upperlimit);

    
  bool checkReachability(const vectord &query);

};
class Circle: public bayesopt::ContinuousModel
  {
  public:
    Circle(bayesopt::Parameters par);
  
    double evaluateSample( const vectord& xin);
  

    
  bool checkReachability(const vectord &query);

};
class TwoCircles: public bayesopt::ContinuousModel
  {
  public:
    TwoCircles(bayesopt::Parameters par);
  
    double evaluateSample( const vectord& xin);
    double smoothstep(double edge0, double edge1, double x);
    double clamp(double x, double lowerlimit, double upperlimit);
    

    
  bool checkReachability(const vectord &query);
  private:
    GaussianNoise gaussianNoise_;
};
class SmoothCircle: public bayesopt::ContinuousModel
  {
  public:
    SmoothCircle(bayesopt::Parameters par);
  
    double evaluateSample( const vectord& xin);
    
    double smoothstep(double edge0, double edge1, double x);
    double clamp(double x, double lowerlimit, double upperlimit);
    
  bool checkReachability(const vectord &query);
  private:
    GaussianNoise gaussianNoise_;
};
#endif
