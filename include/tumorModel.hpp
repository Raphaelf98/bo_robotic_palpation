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
class Triangle: public bayesopt::ContinuousModel
  {
  public:
    Triangle(bayesopt::Parameters par);
  
    double evaluateSample( const vectord& xin);
  

    
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
  

    
  bool checkReachability(const vectord &query);

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
