#ifndef _TUMORMODEL_HPP_
#define _TUMORMODEL_HPP_
#include "bayesopt/bayesopt.hpp"

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
#endif
