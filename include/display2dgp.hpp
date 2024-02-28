
#ifndef _DISPLAY2DGP_HPP_
#define _DISPLAY2DGP_HPP_


#include "matplotpp.h"
#include "tumorModel.hpp"
#include "contour.hpp"
/**
 * @file
 * 
 */
  /**
  * @brief Enum holding states for state machine used by DisplayHeatMap2D class to visualize Algorithm.
  *
  * 
  */
    enum RunningStatus
      {
	    RUN, STEP, STOP, NOT_READY, PLOT_CENTROIDS, PLOT_CONTOUR_POINTS, PLOT_CONTOUR_APPX
      };
   
   /**
  * @brief This class implements a state machine to perform the visualization of the following:
  *  1. Approximation through a gaussian process(GP) of an underlying tumor model
  *  2. Computation of maximum stiffness regions on posterior of GP in order to find tumor centroids
  *  3  Computation of points attributed to tumor contour by exploration of Centroids 
  *  4  Spline approximation of tumor contour
  *
  * 
  */
    class DisplayHeatMap2D :public MatPlot
    { 
    private:

      RunningStatus status; //enum to set current state
      size_t state_ii; // current iteration
      Contour *contour_; // pointer to contour class object

      std::vector<double> lx,ly; //last point
      std::vector<double> cx, cy; //temp variable to store last point
      std::vector<double> cX,cY; 
      std::vector<std::vector<double>> z_, c_, std_; //containers for posterior, eval criterion, standard deviation
      size_t c_points; //c_points x c_points defines size of domain computation is performed on 
      std::vector<std::vector<double> > cZ; //container to store ground truth
      double gTHigh, gTLow, PHigh, PLow, StdHigh, StdLow, CVHigh, CVLow;
      void setHighLow(double &low, double &high, double val);
      //Cluster plot variables
      std::vector<double> clusterx_, clsutery_;
      bool computeClusters; 
      //Contour points plot variables
      bool computeContourPoints;
      std::vector<double> cpointsx_, cpointsy_;
      //Contour approximation plot variables
      bool computeSplinePoints;
      std::vector<std::vector<double>> splinex_, spliney_;
    public:
    /**
    * @brief DisplayHeatMap2D class constructor.
    *
    * @param[in] bopt_model* BayesOptBase model pointer to run optimization on.
    * @param[in] cp Contour parameters used to define functionality of Contour class.
    * @param[in] experiment_path Current working directory path as string.
    * 
    * 
    */
      DisplayHeatMap2D();
     /**
    * @brief Query ground truth at c_points*c_points positions and store result to class member.
    * 
    * 
    */
      void loadGroundTruth();
      /**
    * @brief Initialize Bayesian optimization and instantiate containers for data visualization.
    * @param[in] *contour Contour class object used perform Bayesian optimization, centroid computation, centroid exploration/contour point search and contour approximation.
    * @param[in] dim dimensionality of problem.
    * 
    */
      void init(Contour *contour,size_t dim);

      /**
    * @brief Switch state to STEP optimization. 
    */
      void setSTEP();

        /**
    * @brief Swtich state to RUN. Performs continuous optimization until all samples are queried.
    */
      void toogleRUN();

        /**
    * @brief  State-machine. Perform update on gaussian process until N iterations have passed. Switches states automatically and runs centroid computation, contour point search and contour approximation successively.
    * Data for visualization is stored in class members and later displayed in plotData
    */
      void updateData();

        /**
    * @brief Holds Matplotpp plotting functionality and updates subplots used to display functionality of algorithm.
    */
      void plotData();
        /**
    * @brief Implements Matplotpp's virtual DISPLAY function. Calls updateData and plotData re-iteratively.
    */
      void DISPLAY();
     
    };
  
#endif
