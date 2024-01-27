#include <stdlib.h>

#include "meanShift.hpp"
#include "helper.hpp"
#include"display2dgp.hpp"
#include "dataset.hpp"
#include "prob_distribution.hpp"
 //Constructor
DisplayHeatMap2D::DisplayHeatMap2D(): 
  MatPlot(), cx(1), cy(1), 
  computeClusters(1), 
  computeContourPoints(1), 
  computeSplinePoints(1)
{
  status = NOT_READY;
}

/*
reassignes bounds low and high if val is smaller than low and larger than high.
*/    
void DisplayHeatMap2D::setHighLow(double &low, double &high, double val)
{
  if (val > high)
  {
      high = val;
      
  }
  if (val < low)
  {
      low = val;
     
  }
}
/*
Queries tumor model and stores its values
*/
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
          cZ[i][j]= contour_->evaluateGaussianProcess(q);
          setHighLow(gTLow,gTHigh,cZ[i][j]);
          
      }
      
  }
}
/*
Initialize size of containers, assign contour pointer and call initialization functions
*/
void DisplayHeatMap2D::init(Contour *contour, size_t dim)
{
  contour_ = contour;
  c_points = contour_->getCPoints();
  cX = std::vector<double>(c_points); 
  cY = std::vector<double>(c_points); 
  cZ = std::vector<std::vector<double>>(c_points, std::vector<double>(c_points));
  z_ = std::vector<std::vector<double>>(c_points, std::vector<double>(c_points));
  c_ = std::vector<std::vector<double>>(c_points, std::vector<double>(c_points));
  std_ = std::vector<std::vector<double>>(c_points, std::vector<double>(c_points));
  if (dim != 2) 
  { 
    throw std::invalid_argument("This display only works "
		      "for 2D problems"); 
  } 
  contour_->prepareGaussianProcess(); //initialize gaussian process
 
  loadGroundTruth(); 
   contour_->getInitialSamples(lx,ly);//get initial samples from GP
  
  state_ii = 0;    
  status = STOP;
};
/*
Switch state to step
*/
void DisplayHeatMap2D::setSTEP()
{
  if (status != NOT_READY)
  {
    status = STEP;
  }
}
/*
Sets state to run or stop.
*/
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


void DisplayHeatMap2D::updateData()
{
      
    //1. get new Sample 
    if ((status != STOP) && (state_ii < contour_->number_of_step_runs))
	  {
	      // We are moving. Next iteration
        ++state_ii;
	      contour_->stepRunGaussianProcess(); 
	      const vectord last = contour_->getLastSample();
	      //Update current last sample
	      cx[0] = last(0);
	      cy[0] = last(1);
	      //Update list of all samples collected
	      lx.push_back(last(0));
	      ly.push_back(last(1));
	      if (status == STEP) 
        { 
          status = STOP; 
        }
	  }	    
    /*
    Read posterior, criterion and standard deviation for each iteration. 
    Readjust high and low values for each aspect to adapt colorbar in heatmap later.
    */
    if(state_ii <= contour_->number_of_step_runs)
    {
        std::vector<double> x,y;
        int n=c_points;
        x = linspace(0,1,n);
        y = linspace(0,1,n);
        vectord q(2);
        bayesopt::ProbabilityDistribution* pd = contour_->getPredictionGaussianProcess(q);
        PHigh = pd->getMean();
        PLow = pd->getMean();
        //check 
        StdHigh = 2*pd->getStd();
        StdLow = 2*pd->getStd();
        CVHigh = -contour_->evaluateCriteriaGaussianProcess(q);
        CVLow = -contour_->evaluateCriteriaGaussianProcess(q);
        for(size_t i=0; i<n; ++i)
	      {
        
          for (size_t j=0; j<n; ++j)
          {
            q(1) = y[i]; q(0) = x[j];
            bayesopt::ProbabilityDistribution* pd = contour_->getPredictionGaussianProcess(q);
            c_[i][j] = -contour_->evaluateCriteriaGaussianProcess(q);     //Criteria value
            z_[i][j] = pd->getMean(); //Expected value
            std_[i][j] = pd->getStd();
            setHighLow(PLow,PHigh,z_[i][j]);
            setHighLow(StdLow,StdHigh,std_[i][j]);
            setHighLow(CVLow,CVHigh,c_[i][j]);
          }
	      }
    }
     
    //Update Centroids Data
    if(status == PLOT_CENTROIDS )
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
    // Data points as black star
      status = PLOT_CONTOUR_POINTS;
    }
    //Update Contour Data
    if(status == PLOT_CONTOUR_POINTS)
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
      status = PLOT_CONTOUR_APPX;
    }
    //Update Spline Approximation Data
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
    }
    /*
    Write posterior to file after number of iterations is reached
    */
    if ((status != STOP) && (state_ii == contour_->number_of_step_runs))
    {
      std::cout<<"write file"<<std::endl;
      std::string results_path = RESULTS_PATH;
      std::string posterior_path = FILE_POSTERIOR;
      std::string experiment_dir = contour_->getResultsPath();
      posterior_path = experiment_dir+results_path+posterior_path;
      bayesopt::utils::FileParser fp(posterior_path);
      fp.open(0);
      std::vector<boost::numeric::ublas::vector<double>> ublasVecs = convertStdToUblas(z_);
      fp.write_vecOfvec("posterior", ublasVecs);
      ++state_ii;
      status = PLOT_CENTROIDS;
    }
        

	   
}
/*
Calls plotting funcitons to plot
1. tumor model
2. posterior
3. standard deviation
4. criterion
5. samples
6. centroids
7. contour points
8. contour approximation
*/
void DisplayHeatMap2D::plotData()
{
    title("Press r to run and stop, s to run a step and q to quit.");
    //Plot Data in Ground Truth Plot
    subplot(2,2,2);
      title("Ground Truth");

	    jet();  
      pcolor(cX,cY,cZ);
      // To delete edge lines:
      set("EdgeColor","none");
      // To add color bar
      colorbar(gTLow,gTHigh);
    subplot(2,2,3);
    title("Criteria Value");
     jet();
      
      // To generate pseudo color plot:r
      //pcolor(cX,cY,cZ);
      pcolor(cX,cY,c_);
      // To delete edge lines:
      set("EdgeColor","none");
  
      // To add color bar
      colorbar(CVLow,CVHigh);
      subplot(2,2,4);
    title("Standard Deviation");
     jet();
      
      // To generate pseudo color plot:r
      //pcolor(cX,cY,cZ);
      pcolor(cX,cY,std_);
      // To delete edge lines:
      set("EdgeColor","none");
  
      // To add color bar
      colorbar(StdLow,StdHigh);
      //Plot DATA in Prediction and Samples Plot
    subplot(2,2,1);
    title("Prediction and Samples");
    jet();
    
    plot(cx,cy);set("g");set("o");set(8);        
    //Plot Samples
    if (!lx.empty())
	  {	
	    plot(lx,ly);set("k");set("o");set(4);         // Data points as black star
	  }
    
    //Plot Cluster Centers
    if (!clusterx_.empty() && !clsutery_.empty() )
	  {	
      plot(clusterx_,clsutery_);set("r");set("o");set(4);
    }
    //Plot Contour
    if(!cpointsx_.empty() && !cpointsy_.empty() )
    {
      plot(cpointsx_,cpointsy_);set("r");set("o");set(4);
    }
    //Plot Spline Approximation
    if(!splinex_.empty() && !spliney_.empty() )
    {
      for (size_t i = 0; i < splinex_.size(); i++)
      {
        line(splinex_[i],spliney_[i]);set("k");
      }
    }
      //Plot Prediction
    pcolor(cX,cY,z_);
    set("EdgeColor","none");
    
    // To delete edge lines:
    
    // To add color bar
    colorbar(PLow,PHigh);
        
}
// Implements Matplotpp's virtual DISPLAY function
void DisplayHeatMap2D::DISPLAY()
{
  if (status != NOT_READY)
	{
    updateData();
    plotData();
	}
};



