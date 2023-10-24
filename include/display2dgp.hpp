
#ifndef _DISPLAY2DGP_HPP_
#define _DISPLAY2DGP_HPP_


#include "matplotpp.h"
#include "fileParser.hpp"
#include "tumorModel.hpp"
    
    enum RunningStatus
      {
	    RUN, STEP, STOP, NOT_READY, RUN_MS, 
      };

    class DisplayHeatMap2D :public MatPlot
    { 
    private:
      RunningStatus status;
      size_t state_ii;
      bayesopt::BayesOptBase* bopt_model;
      std::vector<double> lx,ly;
      std::vector<double> cx, cy;
      std::vector<double> solx, soly;
      size_t c_points;
      std::vector<double> cX,cY;
      std::vector<std::vector<double> > cZ;
      double gTHigh, gTLow, PHigh, PLow, StdHigh, StdLow, CVHigh, CVLow;
      void setHighLow(double &low, double &high, double val);
    public:
      DisplayHeatMap2D();
      void setSolution(vectord sol);
      void prepareContourPlot();
      void init(bayesopt::BayesOptBase* bopt, size_t dim);
      void setSTEP();
      void toogleRUN();
      void DISPLAY();
      //Methods for Mean Shift algorithm
     
    };
  
#endif
