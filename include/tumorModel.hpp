#ifndef _TUMORMODEL_HPP_
#define _TUMORMODEL_HPP_
#include "bayesopt/bayesopt.hpp"
#include <random>
#include <functional>
#include"helper.hpp"

/**
 * @brief *Real world measurements are subject to noise. Therefore noise is modeled with zero mean Gauassian noise.

 */
class GaussianNoise{
  public:
    /**
 * @brief Constructor for Gaussian Noise class.
 * @param[in] mean Mean value for Gaussian distribution.
 * @param[in] std Standard deviation of Gaussian distribution.
 */
    GaussianNoise(double mean, double std);
    /**
 * @brief Get noise.
 * @retval Returns floating point noise.
 */
    double noise();
  private:
    std::default_random_engine generator_;
    std::normal_distribution<double> dist_;

};
/**
 * @brief PolarPolygon class.
 * 
 */
class PolarPolygon{
    public:
       /**
 * @brief Default constructor.
 */
    PolarPolygon(){};
       /**
 * @brief Constructor for PolarPolygon class. Constructs either polygon for num_vertices greater 2 or circle for num_vertices smaller or equal to 2.
 * @param[in] num_verices Number of vertices. If number <=2, circle will be constructed, number >2 polygon.
 * @param[in] radius distance from center to vertices or radius of circle.
 * @param[in] x_translation X coordinate of center.
 * @param[in] y_translation Y coordinate of center.
 * @param[in] epsilon TO BE REMOVED
 */
    PolarPolygon(size_t num_verices, double radius, double x_translation, double y_translation, double epsilon);
      /**
 * @brief Constructor for Gaussian Noise class.
 * @param[in] mean Mean value for Gaussian distribution.
 * @param[in] std Standard deviation of Gaussian distribution.
 * @retval Returns floating point noise.
 */
    double evaluate(const double & theta, double rad);
       /**
 * @brief Constructor for Gaussian Noise class.
 * @param[in] mean Mean value for Gaussian distribution.
 * @param[in] std Standard deviation of Gaussian distribution.
 * @retval Returns floating point noise.
 */
    double evaluateCircle(const double & theta, double rad);
       /**
 * @brief Constructor for Gaussian Noise class.
 * @param[in] mean Mean value for Gaussian distribution.
 * @param[in] std Standard deviation of Gaussian distribution.
 * @retval Returns floating point noise.
 */
    double insideCircle(const double &x, const double & y);
       /**
 * @brief Constructor for Gaussian Noise class.
 * @param[in] mean Mean value for Gaussian distribution.
 * @param[in] std Standard deviation of Gaussian distribution.
 * @retval Returns floating point noise.
 */
    double outsideCircle(const double &x, const double & y, double epsilon);
    /**
 * @brief Constructor for Gaussian Noise class.
 * @param[in] mean Mean value for Gaussian distribution.
 * @param[in] std Standard deviation of Gaussian distribution.
 * @retval Returns floating point noise.
 */
    double insidePolygon(const double &x, const double & y);
    /**
 * @brief Constructor for Gaussian Noise class.
 * @param[in] mean Mean value for Gaussian distribution.
 * @param[in] std Standard deviation of Gaussian distribution.
 * @retval Returns floating point noise.
 */
    double outsidePolygon(const double &x, const double & y, double epsilon);
    /**
 * @brief Constructor for Gaussian Noise class.
 
 * @retval Returns floating point noise.
 */
    std::function<double(const double&)> fParametric_x();
        /**
 * @brief Constructor for Gaussian Noise class.
 
 * @retval Returns floating point noise.
 */
    std::function<double(const double&)> fParametric_y();

    private:
    bool circle_;
    double radius_, x_trans_, y_trans_;
    size_t num_vertices_;
    double epsilon_;
    std::vector<std::pair<double,double>> vertices_x_;
    std::vector<std::pair<double,double>> vertices_y_;
    std::vector<double> radians_lookup_;
        /**
 * @brief Constructor for Gaussian Noise class.
 
 * @retval Returns floating point noise.
 */
    double f_x_(const double &s);
        /**
 * @brief Constructor for Gaussian Noise class.
 
 * @retval Returns floating point noise.
 */
    double f_y_(const double &s);
        /**
 * @brief Constructor for Gaussian Noise class.
 
 * @retval Returns floating point noise.
 */
    double polar_x_(const double &theta, const double &rad);
        /**
 * @brief Constructor for Gaussian Noise class.
 
 * @retval Returns floating point noise.
 */
    double polar_y_(const double &theta, const double &rad);
        /**
 * @brief Constructor for Gaussian Noise class.
 
 * @retval Returns floating point noise.
 */
    double r_polar_(const double &x,const double &y);
        /**
 * @brief Constructor for Gaussian Noise class.
 
 * @retval Returns floating point noise.
 */
    double theta_polar_(const double &x,const double &y);
        /**
 * @brief Constructor for Gaussian Noise class.
 
 * @retval Returns floating point noise.
 */
    void computeVertices();
       /**
 * @brief Constructor for Gaussian Noise class.
 
 * @retval Returns floating point noise.
 */
    double rTheta_(const std::pair<double,double> &x,const std::pair<double,double> &y,const double &thetha );
    
};
/**
 * @brief Shape base class.
 * 
 */
class Shape : public bayesopt::ContinuousModel{
  public:
  Shape(bayesopt::Parameters par);
  
  virtual std::function<double (const double&)> f_x()=0;
  virtual std::function<double (const double&)> f_y()=0;
  
  virtual double evaluateSample( const vectord& xin){return 0;}
  double smoothstep(double edge0, double edge1, double x, double low, double high);
  double clamp(double x, double lowerlimit, double upperlimit);
  virtual void saveGroundTruth(const size_t c_points, std::string file_path){}
  protected:
  double low_stiffness_, high_stiffness_;
  double x_trans_, y_trans_, radius_;
  double epsilon_;
};
/**
 * @brief Triangle child class.
 * 
 */
class Triangle: public Shape
  {
    GaussianNoise gaussianNoise_;
    PolarPolygon triangle_;
  public:

    Triangle(bayesopt::Parameters par, double low, double high, double radius, double x_trans, double y_trans, double epsilon,double noise);
  
    double evaluateSample( const vectord& xin);
    virtual std::function<double (const double&)> f_x();
    virtual std::function<double (const double&)> f_y();
    virtual void saveGroundTruth(const size_t c_points, std::string file_path);
    bool checkReachability(const vectord &query);

};
/**
 * @brief Rectangle child class.
 * 
 */
class Rectangle: public Shape
  {
    GaussianNoise gaussianNoise_;
    PolarPolygon rectangle_;
  public:
    Rectangle(bayesopt::Parameters par,double low,double high, double radius, double x_trans, double y_trans, double epsilon, double noise);
  
    double evaluateSample( const vectord& xin);
    virtual std::function<double (const double&)> f_x();
    virtual std::function<double (const double&)> f_y();
    virtual void saveGroundTruth(const size_t c_points, std::string file_path);
    
  bool checkReachability(const vectord &query);

};
/**
 * @brief TwoCircles child class.
 * 
 */
class TwoCircles: public Shape
  {
  public:
    TwoCircles(bayesopt::Parameters par,double low_stiffness, double high_stiffness , double r_1, double r_2,double x_t_1,double x_t_2,
                        double y_t_1,double y_t_2,double epsilon,double noise);
  
    double evaluateSample( const vectord& xin);
    virtual std::function<double (const double&)> f_x();
    virtual std::function<double (const double&)> f_y();
    virtual void saveGroundTruth(const size_t c_points, std::string file_path);
    bool checkReachability(const vectord &query);
  private:
    size_t circle_count_;

    double r_1_, r_2_, x_t_1_, x_t_2_, y_t_1_, y_t_2_,epsilon_ ;
    PolarPolygon circle_1_ ,circle_2_;
    GaussianNoise gaussianNoise_;
};
/**
 * @brief SmoothCircle child class.
 * 
 */
class SmoothCircle: public Shape
  {
  public:
    SmoothCircle(bayesopt::Parameters par,double low,double high, double radius, double x_trans, double y_trans, double epsilon, double noise);

    double evaluateSample( const vectord& xin);
    virtual std::function<double (const double&)> f_x();
    virtual std::function<double (const double&)> f_y();
    virtual void saveGroundTruth(const size_t c_points, std::string file_path);
    bool checkReachability(const vectord &query);
  private:
   
    GaussianNoise gaussianNoise_;
    PolarPolygon circle_;
};
#endif
