#ifndef _TUMORMODEL_HPP_
#define _TUMORMODEL_HPP_
#include "bayesopt/bayesopt.hpp"
#include <random>
#include <functional>
#include"helper.hpp"
/**
 * @file
 * 
 */
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
 * @brief Constructor of PolarPolygon class. Constructs either polygon for num_vertices greater 2 or circle for num_vertices smaller or equal to 2.
 * @param[in] num_vertices Number of vertices. If number <=2, circle will be constructed, number >2 polygon.
 * @param[in] radius distance from center to vertices or radius of circle.
 * @param[in] x_translation X coordinate of center.
 * @param[in] y_translation Y coordinate of center.
 
 */
    PolarPolygon(size_t num_vertices, double radius, double x_translation, double y_translation);

   /**
   * @brief Used for check whether point lies inside polygon or not. 
   * Takes a set of coordinates x and y and computes polar angle from it. 
   * Angle is used to query radius of polygon outline.
   * @param[in] x X coordinate in global reference frame
   * @param[in] y Y coordinate in global reference frame
   * @retval Returns radius(theta) of polygon for given points in cartesian coordinates.
   */
   double polygonRadius(const double &x, const double & y);
  
    /**
   * @brief Returns function pointer to parametric representation of polygon in x.Values range from 0 to 1.
   * 
   * @retval Function pointer of parametric polygon in x.
   */
    std::function<double(const double&)> fParametric_x();
     /**
   * @brief Returns function pointer to parametric representation of polygon in y.Values range from 0 to 1.
   * 
   * @retval Function pointer of parametric polygon in y.
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
   * @param[in] theta Polar angle theta
   * @param[in] rad Polar radius
   * @retval Returns radius r(theta) of polygon.
   */
   double evaluate_(const double & theta, double rad);
    /**
   * @brief Function of which pointer is returned by fParametric_x.
   * @param[in] s Parametric function variable s e [0,1].
   * @retval Returns x coordinate of polygon in absolute coordinates.
   */
    double f_x_(const double &s);
    /**
   * @brief Function of which pointer is returned by fParametric_y.
   * @param[in] s Parametric function variable s e [0,1].
   * @retval Returns y coordinate of polygon in absolute coordinates.
   */
    double f_y_(const double &s);
    /**
   * @brief Computes x coordinate given radius and theta in translated reference frame (center of polygon).
   * @param[in] theta Polar angle theta
   * @param[in] rad Polar radius 
   * @retval Returns x coordinate of polygon in absolute coordinates.
   */
    double polar_x_(const double &theta, const double &rad);
     /**
   * @brief Computes y coordinate given radius and theta in translated reference frame (center of polygon).
   * @param[in] theta Polar angle theta
   * @param[in] rad Polar radius 
   * @retval Returns y coordinate of polygon in absolute coordinates.
   */
    double polar_y_(const double &theta, const double &rad);
   /**
   * @brief Computes polar radius given x and y coordinates in translated reference frame (center of polygon).
   * @param[in] x X coordinate in shape reference frame 
   * @param[in] y Y coordinate in shape reference frame 
   * @retval Returns polar radius
   */
    double r_polar_(const double &x,const double &y);
    /**
   * @brief Computes polar angle given x and y coordinates in translated reference frame (center of polygon).
   * @param[in] x X coordinate in shape reference frame 
   * @param[in] y Y coordinate in shape reference frame 
   * @retval Returns polar angle
   */
    double theta_polar_(const double &x,const double &y);
    /**
   * @brief Function used to compute lookup table with radius segments that connect vertices.
   */
    void computeVertices_();
    /**
   * @brief Takes a pair of point coordinates and outputs r(theta) for straight line connecting the two points.
   * @param[in] x Pair of coordinates x1 and x2 
   * @param[in] y Pair of coordinates y1 and y2  
   * @retval Returns polar r(theta) for straight line connecting point 1 and point 2.
   */
    double rTheta_(const std::pair<double,double> &x, const std::pair<double,double> &y, const double &thetha );
    
};
/**
 * @brief Shape base class. Holds virtual functions that are unique to derived child class and need to be implemented by them.
 * Inherits from bayesopt ContinousModel class.
 * 
*/
class Shape : public bayesopt::ContinuousModel
{
  public:
   /**
   * @brief Constructor of shape base class.
   * @param[in] par Takes bayesopt parameters  
   
   */
  Shape(bayesopt::Parameters par);
  /**
   * @brief Virtual function that is implemented by derived class. 
   * Returns parametric function of tumor model in x as function pointer.
   * @retval Function pointer of parametric function in x.
   */
  virtual std::function<double (const double&)> f_x()=0;
   /**
   * @brief Virtual function that is implemented by derived class. 
   * Returns parametric function of tumor model in y as function pointer.
   * @retval Function pointer of parametric function in y.
   */
  virtual std::function<double (const double&)> f_y()=0;
  /**
   * @brief Function called by bayesian optimization in bayes opt to query tumor model. Virtual Function that is implemented by child classes.
   * When deploying algorithm to hardware, stiffness measurements need to be generated in this function.
   * @param[in] x Takes bayesopt parameters
   * @retval stiffness value/cost of global function that is optimized.
   */
  virtual double evaluateSample(const vectord& xin){return 0;}
  /**
   * @brief Checks reachability of query. Virtual function that is implemented by derived class. 
   * Note: function is not implemented in child classes but might be of value for implementation with hardware.
   * @param[in] query Vector with x and y coordinates.
   * @retval True if query is reachable.
   */
  virtual bool checkReachability(const vectord &query){return 0;}
 /**
   * @brief This function implements a smooth step as a third degree polynomial. The polynomial is constrained through three conditions
   *  1. f(x_1) = y_1 and f(x_2) = y_2
   *  2. f'(x_1) = 0 and f'(x_2) = 0
   *   
   * @param[in] x_l Left x coordinate
   * @param[in] x_r  Right x coordinate
   *  @param[in] x X value a which smooth step is evaluated.
   * @param[in] y_l  Left y coordinate
   *  @param[in] y_r Right y coordinate
   * @retval f(x) function value of smoothstep.
   */
  double smoothstep(const double x_l, const double x_r, const double x, const double y_l, const double y_r );
   /**
   * @brief Virtual function implemented by child class. Saves ground truth to file with resolution of c_points.
   * @param[in] c_points Takes bayesopt parameters
   * @param[in] file_path absolute path to file. Filename is generated in method.
   */
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
    PolarPolygon inner_triangle_, outer_triangle_, groundTruth_triangle_;
   public:
   /**
   * @brief Triangle constructor. Instantiates PolarPolygon with three vertices.
   * @param[in] par Parameters needed to run bayesian optimization.
   * @param[in] low Low stiffness value in Model
   * @param[in] high High stiffness value in Model
   * @param[in] radius Size of Triangle. Radius of circle connecting all vertices.
   * @param[in] x_trans Center of Triangle x coordinate x_trans e [0,1]
   * @param[in] y_trans Center of Triangle y coordinate y_trans e [0,1]
   * @param[in] epsilon Defines smoothness transition from low to high stiffness values. 
   * @param[in] noise Standard deviation of Gaussian noise 
   */
    Triangle(bayesopt::Parameters par, double low, double high, double radius, double x_trans, double y_trans, double epsilon,double noise);
   /**
   * @brief Function called by bayesian optimization in bayes opt to query tumor model.
   * 
   * @param[in] xin Takes bayesopt parameters
   * @retval stiffness value/cost of global function that is optimized.
   */
    double evaluateSample( const vectord& xin);
     /**
   * @brief Returns parametric function of tumor model in x as function pointer.
   * @retval Function pointer of parametric function in x.
   */
    virtual std::function<double (const double&)> f_x();
     /**
   * @brief Returns parametric function of tumor model in y as function pointer.
   * @retval Function pointer of parametric function in y.
   */
    virtual std::function<double (const double&)> f_y();
    /**
   * @brief Saves ground truth to file with resolution of c_points.
   * @param[in] c_points Takes bayesopt parameters
   * @param[in] file_path absolute path to file. Filename is generated in method.
   */
    virtual void saveGroundTruth(const size_t c_points, std::string file_path);
   /**
   * @brief Checks reachability of query. Virtual function that is implemented by derived class. 
   * Note: function is not implemented in child classes but might be of value for implementation with hardware.
   * @param[in] query Vector with x and y coordinates.
   * @retval True if query is reachable.
   */
    virtual bool checkReachability(const vectord &query);

   };
/**
 * @brief Rectangle child class.
 * 
 */
class Rectangle: public Shape
  {
    GaussianNoise gaussianNoise_;
    PolarPolygon inner_rectangle_, outer_rectangle_, groundTruth_rectangle_;
  public:
  /**
   * @brief Rectangle constructor. Instantiates PolarPolygon with four vertices.
   * @param[in] par Parameters needed to run bayesian optimization.
   * @param[in] low Low stiffness value in Model
   * @param[in] high High stiffness value in Model
   * @param[in] radius Size of Rectangle. Radius of circle connecting all vertices.
   * @param[in] x_trans Center of Rectangle x coordinate x_trans e [0,1]
   * @param[in] y_trans Center of Rectangle y coordinate y_trans e [0,1]
   * @param[in] epsilon Defines smoothness transition from low to high stiffness values. 
   * @param[in] noise Standard deviation of Gaussian noise 
   */
    Rectangle(bayesopt::Parameters par,double low,double high, double radius, double x_trans, double y_trans, double epsilon, double noise);
  /**
   * @brief Function called by bayesian optimization in bayes opt to query tumor model.
   * 
   * @param[in] xin Takes bayesopt parameters
   * @retval stiffness value/cost of global function that is optimized.
   */
    double evaluateSample( const vectord& xin);
    /**
   * @brief Returns parametric function of tumor model in x as function pointer.
   * @retval Function pointer of parametric function in x.
   */
    virtual std::function<double (const double&)> f_x();
     /**
   * @brief Returns parametric function of tumor model in y as function pointer.
   * @retval Function pointer of parametric function in y.
   */
    virtual std::function<double (const double&)> f_y();
     /**
   * @brief Saves ground truth to file with resolution of c_points.
   * @param[in] c_points Takes bayesopt parameters
   * @param[in] file_path absolute path to file. Filename is generated in method.
   */
    virtual void saveGroundTruth(const size_t c_points, std::string file_path);
     /**
   * @brief Checks reachability of query. Virtual function that is implemented by derived class. 
   * Note: function is not implemented in child classes but might be of value for implementation with hardware.
   * @param[in] query Vector with x and y coordinates.
   * @retval True if query is reachable.
   */
  bool checkReachability(const vectord &query);

};
/**
 * @brief TwoCircles child class.
 * 
 */
class TwoCircles: public Shape
  {
  public:
   /**
   * @brief Rectangle constructor. Instantiates PolarPolygon with four vertices.
   * @param[in] par Parameters needed to run bayesian optimization.
   * @param[in] low Low stiffness value in Model
   * @param[in] high High stiffness value in Model
   * @param[in] r_1 Size of first circle. 
   * @param[in] r_2 Size of second circle. 
   * @param[in] x_t_1 X coordinate of center of first circle [0,1]
   * @param[in] x_t_2 X coordinate of center of second circle [0,1]
   * @param[in] y_t_1 Y coordinate of center of first circle [0,1]
   * @param[in] y_t_2 Y coordinate of center of second circle [0,1]
   * @param[in] epsilon Defines smoothness transition from low to high stiffness values for both circles. 
   * @param[in] noise Standard deviation of Gaussian noise 
   */
    TwoCircles(bayesopt::Parameters par,double low_stiffness, double high_stiffness , double r_1, double r_2,double x_t_1,double x_t_2,
                        double y_t_1,double y_t_2,double epsilon,double noise);
   /**
   * @brief Function called by bayesian optimization in bayes opt to query tumor model.
   * 
   * @param[in] xin Takes bayesopt parameters
   * @retval stiffness value/cost of global function that is optimized.
   */
    double evaluateSample( const vectord& xin);
    /**
   * @brief Returns parametric function of tumor model in x as function pointer.
   * @retval Function pointer of parametric function in x.
   */
    virtual std::function<double (const double&)> f_x();
   /**
   * @brief Returns parametric function of tumor model in y as function pointer.
   * @retval Function pointer of parametric function in y.
   */
    virtual std::function<double (const double&)> f_y();
   /**
   * @brief Saves ground truth to file with resolution of c_points.
   * @param[in] c_points Takes bayesopt parameters
   * @param[in] file_path absolute path to file. Filename is generated in method.
   */
    virtual void saveGroundTruth(const size_t c_points, std::string file_path);
   /**
   * @brief Checks reachability of query. Virtual function that is implemented by derived class. 
   * Note: function is not implemented in child classes but might be of value for implementation with hardware.
   * @param[in] query Vector with x and y coordinates.
   * @retval True if query is reachable.
   */
    bool checkReachability(const vectord &query);
  private:
    size_t circle_count_;

    double r_1_, r_2_, x_t_1_, x_t_2_, y_t_1_, y_t_2_,epsilon_ ;
    PolarPolygon inner_circle_1_ ,inner_circle_2_;
    PolarPolygon outer_circle_1_ ,outer_circle_2_;
    PolarPolygon groundTruth_circle_1_ ,groundTruth_circle_2_;
    GaussianNoise gaussianNoise_;
};
/**
 * @brief Circle child class.
 * 
 */
class Circle: public Shape
  {
  public:
  /**
   * @brief Circle constructor. Instantiates PolarPolygon with less than 3 vertices (circle).
   * @param[in] par Parameters needed to run bayesian optimization.
   * @param[in] low Low stiffness value in Model
   * @param[in] high High stiffness value in Model
   * @param[in] radius Size of Rectangle. Radius of circle connecting all vertices.
   * @param[in] x_trans Center of Rectangle x coordinate x_trans e [0,1]
   * @param[in] y_trans Center of Rectangle y coordinate y_trans e [0,1]
   * @param[in] epsilon Defines smoothness transition from low to high stiffness values. 
   * @param[in] noise Standard deviation of Gaussian noise 
   */
    Circle(bayesopt::Parameters par,double low,double high, double radius, double x_trans, double y_trans, double epsilon, double noise);
   /**
   * @brief Function called by bayesian optimization in bayes opt to query tumor model.
   * 
   * @param[in] xin Takes bayesopt parameters
   * @retval stiffness value/cost of global function that is optimized.
   */
    double evaluateSample( const vectord& xin);
     /**
   * @brief Returns parametric function of tumor model in x as function pointer.
   * @retval Function pointer of parametric function in x.
   */
    virtual std::function<double (const double&)> f_x();
    /**
   * @brief Returns parametric function of tumor model in y as function pointer.
   * @retval Function pointer of parametric function in y.
   */
    virtual std::function<double (const double&)> f_y();
    /**
   * @brief Saves ground truth to file with resolution of c_points.
   * @param[in] c_points Takes bayesopt parameters
   * @param[in] file_path absolute path to file. Filename is generated in method.
   */
    virtual void saveGroundTruth(const size_t c_points, std::string file_path);
    /**
   * @brief Checks reachability of query. Virtual function that is implemented by derived class. 
   * Note: function is not implemented in child classes but might be of value for implementation with hardware.
   * @param[in] query Vector with x and y coordinates.
   * @retval True if query is reachable.
   */
    bool checkReachability(const vectord &query);
  private:
   
    GaussianNoise gaussianNoise_;
    PolarPolygon inner_circle_, outer_circle_,groundTruth_circle_;
};
#endif
