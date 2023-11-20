#ifndef _TUMORMODEL_HPP_
#define _TUMORMODEL_HPP_
#include "bayesopt/bayesopt.hpp"
#include <random>
#include <functional>
/*Real world measurements are subject to noise. Therefore noise is modeled with zero mean Gauassian noise*/

class GaussianNoise{
  public:
    GaussianNoise(double mean, double std);
    double noise();
  private:
    std::default_random_engine generator_;
    std::normal_distribution<double> dist_;

};
/**/
class PolarPolygon{
    public:
    PolarPolygon(){};
    PolarPolygon(size_t num_verices, double radius, double x_translation, double y_translation, double epsilon);
    double evaluate(const double & theta, double rad);
    double evaluateCircle(const double & theta, double rad);
    double insideCircle(const double &x, const double & y);
    double outsideCircle(const double &x, const double & y, double epsilon);
    double insidePolygon(const double &x, const double & y);
    double outsidePolygon(const double &x, const double & y, double epsilon);
    std::function<double(const double&)> fParametric_x();
    std::function<double(const double&)> fParametric_y();
    private:
    double radius_, x_trans_, y_trans_;
    size_t num_vertices_;
    double epsilon_;
    std::vector<std::pair<double,double>> vertices_x_;
    std::vector<std::pair<double,double>> vertices_y_;
    std::vector<double> radians_lookup_;
    double f_x_(const double &s);
    double f_y_(const double &s);
    double polar_x_(const double &theta, const double &rad);
    double polar_y_(const double &theta, const double &rad);
    double r_polar_(const double &x,const double &y);
    double theta_polar_(const double &x,const double &y);
    void computeVertices();
   
    double rTheta_(const std::pair<double,double> &x,const std::pair<double,double> &y,const double &thetha );
    
};

class Shape : public bayesopt::ContinuousModel{
  public:
  Shape(bayesopt::Parameters par);
  
  virtual std::function<double (const double&)> f_x()=0;
  virtual std::function<double (const double&)> f_y()=0;
  virtual double evaluateSample( const vectord& xin){return 0;}
  double smoothstep(double edge0, double edge1, double x, double low, double high);
  double clamp(double x, double lowerlimit, double upperlimit);
  
};
class Triangle: public Shape
  {
    GaussianNoise gaussianNoise_;
    PolarPolygon triangle_;
  public:

    Triangle(bayesopt::Parameters par);
  
    double evaluateSample( const vectord& xin);
    virtual std::function<double (const double&)> f_x();
    virtual std::function<double (const double&)> f_y();
    
  bool checkReachability(const vectord &query);

};
class Rectangle: public Shape
  {
    GaussianNoise gaussianNoise_;
    PolarPolygon rectangle_;
  public:
    Rectangle(bayesopt::Parameters par);
  
    double evaluateSample( const vectord& xin);
    virtual std::function<double (const double&)> f_x();
    virtual std::function<double (const double&)> f_y();

    
  bool checkReachability(const vectord &query);

};
/*
class Circle: public Shape
  {
  public:
    Circle(bayesopt::Parameters par);
  
    double evaluateSample( const vectord& xin);
  virtual std::function<double (const double&)> f_x();
    virtual std::function<double (const double&)> f_y();

    
  bool checkReachability(const vectord &query);

};*/
class TwoCircles: public Shape
  {
  public:
    TwoCircles(bayesopt::Parameters par,double r_1, double r_2,double x_t_1,double x_t_2,
                        double y_t_1,double y_t_2,double epsilon);
  
    double evaluateSample( const vectord& xin);
    virtual std::function<double (const double&)> f_x();
    virtual std::function<double (const double&)> f_y();
    
  bool checkReachability(const vectord &query);
  private:
  double r_1_, r_2_, x_t_1_, x_t_2_, y_t_1_, y_t_2_,epsilon_ ;
  PolarPolygon circle_1_ ,circle_2_;
    GaussianNoise gaussianNoise_;
};
class SmoothCircle: public Shape
  {
  public:
    SmoothCircle(bayesopt::Parameters par);
  
    double evaluateSample( const vectord& xin);
   virtual std::function<double (const double&)> f_x();
    virtual std::function<double (const double&)> f_y();
    
  bool checkReachability(const vectord &query);
  private:
    GaussianNoise gaussianNoise_;
    PolarPolygon circle_;
};
#endif
