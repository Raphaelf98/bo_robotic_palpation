
#ifndef _DISPLAY2DGP_HPP_
#define _DISPLAY2DGP_HPP_


#include "matplotpp.h"
#include "tumorModel.hpp"
#include "contour.hpp"
    /*
    States for state-machine 
    */
    enum RunningStatus
      {
	    RUN, STEP, STOP, NOT_READY, PLOT_CENTROIDS, PLOT_CONTOUR_POINTS, PLOT_CONTOUR_APPX
      };
    /*
    This class implements a state machine to perform the visualization of the following:
      1. Approximation through a gaussian process(GP) of an underlying tumor model
      2. Computation of maximum stiffness regions on posterior of GP in order to find tumor centroids
      3  Computation of points attributed to tumor contour by exploration of Centroids 
      4  Spline approximation of tumor contour

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
      size_t c_points; //c_points x c_points defines size of domain computation is perfromed on 
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
      DisplayHeatMap2D();
     
      void loadGroundTruth();
      void init(Contour *contour,size_t dim);
      void setSTEP();
      void toogleRUN();
      void updateData();
      void plotData();
      void DISPLAY();
     
    };
  
#endif
