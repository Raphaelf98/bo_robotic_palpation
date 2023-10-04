#include"display2dgp.hpp"
#include "dataset.hpp"
#include "prob_distribution.hpp"
#include <stdlib.h>
   


  Triangle::Triangle(bayesopt::Parameters par):
    ContinuousModel(2,par) {}

  double Triangle::evaluateSample( const vectord& xin)
  {
     if (xin.size() != 2)
      {
	std::cout << "WARNING: This only works for 2D inputs." << std::endl
		  << "WARNING: Using only first two components." << std::endl;
      }
    double x = xin(0);
    double y = xin(1);
    
    if (y >= 0.5 && y <= 1.5*x && y <= -1.5*x +1.5)
    {
      return 0.0;
    }
    if (y >= 0.4 && y <= 1.5*x+0.2 && y <= -1.5*x +1.7)
    {
      return 5.0;
    }
    else{
      return 10.0;
    }
    
    
   //return (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5);
    
  }

  bool Triangle::checkReachability(const vectord &query)
  {return true;};










    DisplayHeatMap2D::DisplayHeatMap2D(): 
      MatPlot(), cx(1), cy(1), c_points(100), cX(c_points),
      cY(c_points), cZ(c_points,std::vector<double>(c_points))
    {
      status = NOT_READY;
    }

    void DisplayHeatMap2D::setSolution(vectord sol)
    {
      solx.push_back(sol(0));
      soly.push_back(sol(1));
    }

    void DisplayHeatMap2D::prepareContourPlot()
    {
      cX=linspace(0,1,c_points);
      cY=linspace(0,1,c_points);
      
      for(int i=0;i<c_points;++i)
	    {
	        for(int j=0;j<c_points;++j)
	        {
	            vectord q(2);
	            q(0) = cX[j]; q(1) = cY[i];
	            cZ[i][j]= bopt_model->evaluateSample(q);
	        }
	    }
    }

    void DisplayHeatMap2D::init(bayesopt::BayesOptBase* bopt, size_t dim)
    {
      if (dim != 2) 
	{ 
	  throw std::invalid_argument("This display only works "
				      "for 2D problems"); 
	}
      
      bopt_model = bopt;
      prepareContourPlot();
      
      bopt->initializeOptimization();
      size_t n_points = bopt->getData()->getNSamples();
      for (size_t i = 0; i<n_points;++i)
	    {
	      const vectord last = bopt->getData()->getSampleX(i);
	      lx.push_back(last(0));
	      ly.push_back(last(1));
	    }
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
	    size_t nruns = bopt_model->getParameters()->n_iterations;
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
      colorbar();
      subplot(2,2,1);
      title("Prediction and Samples");
      plot(cx,cy);set("g");set("o");set(4);         // Data points as black star
	    plot(solx,soly);set("r"); set("o");set(4); 
    
        
	    if ((status != STOP) && (state_ii < nruns))
	    {
	      // We are moving. Next iteration
	      ++state_ii;
          
	      bopt_model->stepOptimization(); 
	      const vectord last = bopt_model->getData()->getLastSampleX();
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
        std::vector<std::vector<double>> z{c_points,std::vector<double>(c_points)};
        std::vector<std::vector<double>> c{c_points,std::vector<double>(c_points)};
    
        for(size_t i=0; i<n; ++i)
	      {
        
          for (size_t j=0; j<n; ++j)
          {
            q(1) = y[i]; q(0) = x[j];
            bayesopt::ProbabilityDistribution* pd = bopt_model->getPrediction(q);
            c[i][j] = -bopt_model->evaluateCriteria(q);     //Criteria value
            z[i][j] = pd->getMean(); //Expected value
          }

	      }
      
       jet();
        
        // To generate pseudo color plot:r
        //pcolor(cX,cY,cZ);
        pcolor(cX,cY,z);

        // To delete edge lines:
        set("EdgeColor","none");
    
        // To add color bar
        colorbar();

      subplot(2,2,3);
      title("Criteria Value");
       jet();
        
        // To generate pseudo color plot:r
        //pcolor(cX,cY,cZ);
        pcolor(cX,cY,c);

        // To delete edge lines:
        set("EdgeColor","none");
    
        // To add color bar
        colorbar();
        

	}
    };



