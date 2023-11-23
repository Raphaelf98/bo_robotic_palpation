#include"display2dgp.hpp"
#include "dataset.hpp"
#include "prob_distribution.hpp"
#include <stdlib.h>
#include "meanShift.hpp"


 




    DisplayHeatMap2D::DisplayHeatMap2D(): 
      MatPlot(), cx(1), cy(1), c_points(100), cX(c_points), 
      cY(c_points), cZ(c_points,std::vector<double>(c_points)), computeClusters(1), computeContourPoints(1), computeSplinePoints(1)
    {
      status = NOT_READY;
      
    }

    void DisplayHeatMap2D::setSolution(vectord sol)
    {
      solx.push_back(sol(0));
      soly.push_back(sol(1));
    }
    void DisplayHeatMap2D::setHighLow(double &low, double &high, double val){
      if (val > high)
      {
          high = val;
      }
      if (val < low)
      {
          low = val;
      }
    }
    void DisplayHeatMap2D::loadGroundTruth()
    {
      cX=linspace(0,1,c_points);
      cY=linspace(0,1,c_points);
      vectord q_(2);
      q_(0) = cX[0]; 
      q_(1) = cY[0];
      gTLow = contour_->evaluateGaussianProcess(q_);
      gTHigh = contour_->evaluateGaussianProcess(q_);

      for(int i=0;i<c_points;++i)
	    {  
	        for(int j=0;j<c_points;++j)
	        {
	            vectord q(2);
	            q(0) = cX[j]; q(1) = cY[i];
	            cZ[i][j]= contour_->evaluateGaussianProcess(q );
              setHighLow(gTLow,gTHigh,cZ[i][j]);
              
	        }
          
	    }
    }

    void DisplayHeatMap2D::init(Contour *contour, size_t dim)
    {
      contour_ = contour;
      if (dim != 2) 
	    { 
	      throw std::invalid_argument("This display only works "
				      "for 2D problems"); 
	    } 
      contour_->prepareGaussianProcess();
     
      loadGroundTruth();
       contour_->getInitialSamples(lx,ly);
  
      
      state_ii = 0;    
      status = STOP;
    };

    void DisplayHeatMap2D::setSTEP()
    {
      if (status != NOT_READY)
	    {
	        status = STEP;
	    }
    };

    void DisplayHeatMap2D::toogleRUN()
    {
      if (status != NOT_READY)
	  {
	    if(status != RUN)
	    {
	      status = RUN;
	    }
	    else
	    {
	      status = STOP;
	    }
	  }
    }
    
   
    void DisplayHeatMap2D::DISPLAY()
    {

      if (status != NOT_READY)
	  {
      std::vector<std::vector<double>> z{c_points,std::vector<double>(c_points)};
	    
	    title("Press r to run and stop, s to run a step and q to quit.");
      subplot(2,2,2);
      title("Ground Truth");
	    //contour(cX,cY,cZ,50);                         // Contour plot (50 lines)
	    jet();
        
      // To generate pseudo color plot:r
      pcolor(cX,cY,cZ);
      //pcolor(cX,cY,z);

      // To delete edge lines:
      set("EdgeColor","none");
    
      // To add color bar
      colorbar(gTLow,gTHigh);
      subplot(2,2,1);
      title("Prediction and Samples");
      plot(cx,cy);set("g");set("o");set(4);         // Data points as black star
	    //plot(solx,soly);set("r"); set("o");set(4); 
    
        
	    if ((status != STOP) && (state_ii < contour_->number_of_step_runs))
	    {
	      // We are moving. Next iteration
	      ++state_ii;
          
	      contour_->stepRunGaussianProcess(); 
        
	      const vectord last = contour_->getLastSample();
	      //GP subplot
	      cx[0] = last(0);
	      cy[0] = last(1);
	      
	      if (!lx.empty())
		    {	
		      plot(lx,ly);set("k");set("o");set(4);         // Data points as black star
		    }
	      
	      lx.push_back(last(0));
	      ly.push_back(last(1));
	      
	      if (status == STEP) 
        { 
          status = STOP; 
        }
	    }	    
	      else
	      {
          plot(lx,ly);set("k");set("o");set(4);         // Data points as black star
         
	      }
        std::vector<double> x,y;
        int n=c_points;
        x = linspace(0,1,n);
        y = linspace(0,1,n);
        vectord q(2);
        bayesopt::ProbabilityDistribution* pd = contour_->getPredictionGaussianProcess(q);
        
        PHigh = pd->getMean();
        PLow = pd->getMean();
        std::vector<std::vector<double>> c{c_points,std::vector<double>(c_points)};
        //check 
        StdHigh = 2*pd->getStd();
        StdLow = 2*pd->getStd();
        std::vector<std::vector<double>> std{c_points,std::vector<double>(c_points)};
        CVHigh = -contour_->evaluateCriteriaGaussianProcess(q);
        CVLow = -contour_->evaluateCriteriaGaussianProcess(q);
        for(size_t i=0; i<n; ++i)
	      {
        
          for (size_t j=0; j<n; ++j)
          {
            q(1) = y[i]; q(0) = x[j];
            bayesopt::ProbabilityDistribution* pd = contour_->getPredictionGaussianProcess(q);
            c[i][j] = -contour_->evaluateCriteriaGaussianProcess(q);     //Criteria value
            z[i][j] = pd->getMean(); //Expected value
            std[i][j] = pd->getStd();
            setHighLow(PLow,PHigh,z[i][j]);
            setHighLow(StdLow,StdHigh,std[i][j]);
            setHighLow(CVLow,CVHigh,c[i][j]);
          }

	      }

        if ((status != STOP) && (state_ii == contour_->number_of_step_runs))
        {
          std::cout<<"write file"<<std::endl;
            FileParser fp("/home/raphael/robolab/displaygp/build/posterior.txt");
          fp.open(0);
          fp.write_stdvecOfvec("posterior", z);
          ++state_ii;
          status = PLOT_CENTROIDS;
        }

        if(status == PLOT_CENTROIDS | status == PLOT_CONTOUR_POINTS | status==PLOT_CONTOUR_APPX)
        {
          if(computeClusters)
          {
            contour_->computeCluster();
            computeClusters=false;
          }
          std::vector<Point> clusters = contour_->getClusters();
          for (auto &c : clusters)
          {
            clusterx_.push_back(c.x);
            clsutery_.push_back(c.y);

          }
          plot(clusterx_,clsutery_);set("r");set("o");set(4);         // Data points as black star
          status = PLOT_CONTOUR_POINTS;
        }
        if(status == PLOT_CONTOUR_POINTS | status==PLOT_CONTOUR_APPX)
        {
          if(computeContourPoints)
          {
            contour_->exploreContour();
            computeContourPoints = false;
          }
            std::vector<Point> cpoints = contour_->getContourPoints();
          for (auto &c : cpoints)
          {
            cpointsx_.push_back(c.x);
            cpointsy_.push_back(c.y);

          }
          plot(cpointsx_,cpointsy_);set("r");set("o");set(4);
          status = PLOT_CONTOUR_APPX;

        }
        if (status == PLOT_CONTOUR_APPX)
        {
          if(computeSplinePoints)
          {
            contour_->approximateContour();
            computeSplinePoints = false;

            SplineInterpolant_ptr_pair_vec spline_pairs = contour_->getSplineInterpolant();
            
            size_t num_vertices = 100;
            splinex_ = std::vector<std::vector<double>>(spline_pairs.size(), std::vector<double>(num_vertices,0));
            spliney_ = std::vector<std::vector<double>>(spline_pairs.size(), std::vector<double>(num_vertices,0));
            for (size_t i = 0; i < spline_pairs.size(); i++)
            {
              std::vector<double> tmp;
               std::shared_ptr<alglib::spline1dinterpolant> s1 = spline_pairs[i].first;
               std::shared_ptr<alglib::spline1dinterpolant> s2= spline_pairs[i].second;
               for (size_t j = 0; j < num_vertices; j++)
               {   
                 //counter-clockwise 
                 double x = spline1dcalc(*s1,  1.0-j * 1.0/((double)num_vertices));
                 double y = spline1dcalc(*s2,  1.0-j * 1.0/((double)num_vertices));
                 splinex_[i][j]=x;
                 spliney_[i][j]=y;
                 
             }
            }
            
          }
          for (size_t i = 0; i < splinex_.size(); i++)
            {
              line(splinex_[i],spliney_[i]);set("k");
            }
        }
       
        jet();
        
        // To generate pseudo color plot:r
        //pcolor(cX,cY,cZ);
        pcolor(cX,cY,z);

        // To delete edge lines:
        set("EdgeColor","none");
    
        // To add color bar
        colorbar(PLow,PHigh);

      subplot(2,2,3);
      
      title("Criteria Value");
       jet();
        
        // To generate pseudo color plot:r
        //pcolor(cX,cY,cZ);
        pcolor(cX,cY,c);

        // To delete edge lines:
        set("EdgeColor","none");
    
        // To add color bar
        colorbar(CVLow,CVHigh);
        subplot(2,2,4);
      title("Standard Deviation");
       jet();
        
        // To generate pseudo color plot:r
        //pcolor(cX,cY,cZ);
        pcolor(cX,cY,std);

        // To delete edge lines:
        set("EdgeColor","none");
    
        // To add color bar
        colorbar(StdLow,StdHigh);
        

	}

    };



