
#ifndef _DISPLAY2DGP_HPP_
#define _DISPLAY2DGP_HPP_


#include "matplotpp.h"
#include "fileParser.hpp"
#include "tumorModel.hpp"
#include "contour.hpp"

    enum RunningStatus
      {
	    RUN, STEP, STOP, NOT_READY, PLOT_CENTROIDS, PLOT_CONTOUR_POINTS, PLOT_CONTOUR_APPX
      };

    class DisplayHeatMap2D :public MatPlot
    { 
    private:
      RunningStatus status;
      size_t state_ii;
      Contour *contour_;
      std::vector<double> lx,ly;
      std::vector<double> cx, cy;
      std::vector<double> solx, soly;
      size_t c_points;
      std::vector<double> cX,cY;
      std::vector<std::vector<double> > cZ;
      double gTHigh, gTLow, PHigh, PLow, StdHigh, StdLow, CVHigh, CVLow;
      void setHighLow(double &low, double &high, double val);
      //Cluster plot
      std::vector<double> clusterx_, clsutery_;
      bool computeClusters;
      //Contour points plot
      bool computeContourPoints;
      std::vector<double> cpointsx_, cpointsy_;
      //Contour approximation plot
      bool computeSplinePoints;
      std::vector<std::vector<double>> splinex_, spliney_;
    public:
      DisplayHeatMap2D();
      void setSolution(vectord sol);
      void loadGroundTruth();
      void init(Contour *contour,size_t dim);
      void setSTEP();
      void toogleRUN();
      
      void DISPLAY();
      //Methods for Mean Shift algorithm
     
    };
  
#endif
