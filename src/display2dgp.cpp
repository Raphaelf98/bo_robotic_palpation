#include"display2dgp.hpp"
#include "dataset.hpp"
#include "prob_distribution.hpp"
#include <stdlib.h>
#include "meanShift.hpp"


 




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
    void DisplayHeatMap2D::prepareContourPlot()
    {
      cX=linspace(0,1,c_points);
      cY=linspace(0,1,c_points);
      vectord q_(2);
      q_(0) = cX[0]; 
      q_(1) = cY[0];
      gTLow = bopt_model->evaluateSample(q_);
      gTHigh = bopt_model->evaluateSample(q_);

      for(int i=0;i<c_points;++i)
	    {   std::cout<<"[";
	        for(int j=0;j<c_points;++j)
	        {
	            vectord q(2);
	            q(0) = cX[j]; q(1) = cY[i];
	            cZ[i][j]= bopt_model->evaluateSample(q);
              setHighLow(gTLow,gTHigh,cZ[i][j]);
              std::cout<<cZ[i][j]<< ", ";
	        }
          std::cout<<"]"<<std::endl;
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
      std::vector<std::vector<double>> z{c_points,std::vector<double>(c_points)};
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
      colorbar(gTLow,gTHigh);
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
        bayesopt::ProbabilityDistribution* pd = bopt_model->getPrediction(q);
        
        PHigh = pd->getMean();
        PLow = pd->getMean();
        std::vector<std::vector<double>> c{c_points,std::vector<double>(c_points)};
        StdHigh = 2*pd->getStd();
        StdLow = 2*pd->getStd();
        std::vector<std::vector<double>> std{c_points,std::vector<double>(c_points)};
        CVHigh = -bopt_model->evaluateCriteria(q);
        CVLow = -bopt_model->evaluateCriteria(q);
        for(size_t i=0; i<n; ++i)
	      {
        
          for (size_t j=0; j<n; ++j)
          {
            q(1) = y[i]; q(0) = x[j];
            bayesopt::ProbabilityDistribution* pd = bopt_model->getPrediction(q);
            c[i][j] = -bopt_model->evaluateCriteria(q);     //Criteria value
            z[i][j] = pd->getMean(); //Expected value
            std[i][j] = pd->getStd();
            setHighLow(PLow,PHigh,z[i][j]);
            setHighLow(StdLow,StdHigh,std[i][j]);
            setHighLow(CVLow,CVHigh,c[i][j]);
          }

	      }
        if ((status != STOP) && (state_ii == nruns))
        {
          std::cout<<"write file"<<std::endl;
            FileParser fp("/home/raphael/robolab/displaygp/build/posterior.txt");
          fp.open(0);
          fp.write_stdvecOfvec("posterior", z);
          ++state_ii;
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



