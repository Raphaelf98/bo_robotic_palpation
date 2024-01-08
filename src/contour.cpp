#include "contour.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <fstream>
#include "dataset.hpp"
#include "prob_distribution.hpp"
#include "fileParser.hpp"
#include "param_loader.hpp"
#include<bayesopt/parameters.hpp>
/*
TODO
-manage parameter parsing for:
    -n_exploration_directions_
    -c_points_
    -bandwidth_
    -lim_steps_
    -multiplier_
*/
Contour::Contour(bayesopt::BayesOptBase* bopt_model, size_t n_exploration_directions):
        bopt_model_(bopt_model),
        n_exploration_directions_(n_exploration_directions),
        c_points_(100),
        bandwidth_(0.05),
        samples_(c_points_),
        n_directions_(n_exploration_directions),
        //stiffness_threshold_(0.01),
        lim_steps_(1000),
        n_samples_(n_directions_+1),
        c_(c_points_,std::vector<double>(c_points_)),
        multiplier_(3)

{
    number_of_step_runs = bopt_model_->getParameters()->n_iterations ;
    total_number_of_iterations_ = bopt_model_->getParameters()->n_init_samples + number_of_step_runs;
    y_values_.reserve(total_number_of_iterations_);
     bayesopt::utils::ParamLoader::save("parameters_stored.txt", *(bopt_model->getParameters()));
}


Contour::~Contour()
{
}
void Contour::runGaussianProcess(){
    size_t state_ii = 0;
    
    bopt_model_->initializeOptimization();
    size_t nruns = bopt_model_->getParameters()->n_iterations;
    
    while(state_ii < nruns)
    {
        ++state_ii;
        bopt_model_->stepOptimization();

    }

}
void Contour::computeCluster()
{   

    std::vector<double> x,y;
    vectord q(2);
    x = linSpace(0,1,c_points_);
    y = linSpace(0,1,c_points_);
    for(size_t i=0; i<c_points_; ++i)
	{
          for (size_t j=0; j<c_points_; ++j)
          {
            q(1) = y[i]; q(0) = x[j];
            bayesopt::ProbabilityDistribution* pd = bopt_model_->getPrediction(q);
            c_[i][j] = pd->getMean();   
          }
    }
    mean_shift_ = MeanShift(c_, bandwidth_, samples_);
    mean_shift_.meanshift_mlpack();

    std::vector<std::vector<double>> vec = mean_shift_.getCentroids();
    clusters_ = std::vector<Point>(vec[0].size());
    
    //TODO: make more elegant
    for(size_t j = 0; j<vec[0].size(); j++ )
    {
        
        clusters_[j].x = vec[0][j];
        clusters_[j].y = vec[1][j];
        std::cout<<"Cluster #"<<j<< " X: "<<clusters_[j].x<< " Y: "<<clusters_[j].y<<std::endl;
    }
    
}
std::vector<Point> Contour::getClusters()
{   
    return clusters_;

}
void Contour::labelData_()
{
    //collect all sampled stiffness values from gaussian process
    bayesopt::Dataset* data = const_cast<bayesopt::Dataset*>(bopt_model_->getData());
    for(size_t i = 0; i<total_number_of_iterations_; i++)
    {
        y_values_.push_back(data->getSampleY(i));
    }
    //cluster them to filter stiffness values associated with tumor
    K_means kmeans(y_values_);
    kmeans.cluster();
    //identify label for centroid with lowest stiffness values
    std::vector<double> c = kmeans.getCentroids();
    //Assumes that centroid index corresponds to label
    std::vector<double>::iterator min_idx_it = std::min_element(c.begin(), c.end());
    size_t min_idx = std::distance(c.begin(), min_idx_it);
    std::cout<<"Label/Index for max stiffness: "<<min_idx<<std::endl;
    for (auto &i : c)
    {
        std::cout<<"Centroid # "<<i<< std::endl;
        
    }
    //retrieve cluster center of tumor stiffness
    tumor_stiffness_mean_ = c[min_idx];
    std::cout<<"Tumor Stiffness mean is  "<<tumor_stiffness_mean_<< std::endl;
    std::vector<std::pair<double,size_t>> labels = kmeans.getAssignments();
    for(auto &a : labels){
        
        if(a.second == min_idx){

            tumor_stiffness_vec_.push_back(a.first);
            std::cout<<"Stiffness in Tumor cluster "<< a.first<<std::endl;
        }
        else
        {
            std::cout<<"Surrounding stiffness"<< a.first<<std::endl;
        }
        
    }
    
}
void Contour::computeStiffnessThreshold_()
{
     std::cout<<"STD: "<<stdDev(tumor_stiffness_vec_)<< std::endl;
    threshold_ = multiplier_*stdDev(tumor_stiffness_vec_)+tumor_stiffness_mean_;
    std::cout<<"THRESHOLD: "<< threshold_<< std::endl;
}
bool Contour::contourPoint_(double &stiffness)
{
    return stiffness >= threshold_;
}
void Contour::exploreContour(){
    /*
    TODO: 
    */
    labelData_();
    computeStiffnessThreshold_();

    for(auto &p1 : clusters_)
    {
        
        double stepsize = 0.003;
        
        double delta_theta = 360 /static_cast<double>(n_directions_);
        //explore
        vectord q1(2);
	    q1(0) = p1.x; q1(1) = p1.y;

	    std::cout<<"Compute contour points for center at: X: " <<p1.x <<" Y: "<<p1.y<<std::endl;
    
        double theta = 0.0;
        std::vector<Point> contour_vec{n_directions_+1};
        std::vector<std::vector<double>> stiffness_vec{n_directions_+1};
        for (size_t i = 0; i < n_directions_+1; i++)
        {   
            theta = delta_theta* i;
            double delta_x = stepsize * sin(theta*M_PI/180.0);
            double delta_y = stepsize * cos(theta*M_PI/180.0);
            int j = 0;
            Point p = p1;
            double old_stiffness;
            int contour_position_idx;

            while(j < lim_steps_)
            {   
                p.x +=  delta_x;
                p.y +=  delta_y;

                vectord q(2);
                q(0) = p.x; 
                q(1) = p.y;
                double new_stiffness = bopt_model_->evaluateSample(q);
                stiffness_vec[i].push_back(new_stiffness);
                if (j>0)
                {
                    if (contourPoint_(new_stiffness))
                    //if(new_stiffness > 1.5)
                    {   
                        contour_vec[i] = p;
                        std::cout<<"Contour point at "<<theta<< " at: X: " <<p.x <<" Y: "<<p.y<< " with step: dx / dy "<< delta_x<< " / " << delta_y<<std::endl;

                        contour_position_idx = j;
                        j = lim_steps_ +1;
                    }
                }
                j++;
                if(j == lim_steps_)
                {
                    std::cout<<"No contour point for theta at "<< theta<<std::endl;
                }
                old_stiffness = new_stiffness;
            }



        }
        contours_.push_back(contour_vec);
    }
}
std::vector<Point> Contour::getContourPoints()
{
    std::vector<Point> contour_points;

    for(auto &c : contours_)
    { 
        contour_points.insert(contour_points.end(), c.begin(),c.end());
    }
    return contour_points;
}
void Contour::approximateContour()
{
    int file_num = 1;
    //Iterate over all 
    //TODO implement file write in extra function
    for(std::vector<Point> &contour : contours_)
    {
        double t[n_samples_];
        std::cout<<"sample size: "<<sizeof(t)/sizeof(double)<<std::endl;
        double f_1[n_samples_];
        double f_2[n_samples_];
        std::string file_contour_points = "/home/raphael/robolab/displaygp/config/contour_points_" + std::to_string(file_num) + ".csv";
        std::ofstream file(file_contour_points);
        if (!file.is_open()) 
        {
            std::cerr << "Failed to open the file." << std::endl;
        }

        file << "X,Y\n";

        for (size_t i = 0; i < n_samples_; i++)
        {
            f_1[i] =contour[i].x;
            f_2[i] =contour[i].y;
            t[i] = i * 1.0/((double)n_samples_-1);
            std::cout<<"sample # "<<t[i] << " at X: " <<f_1[i]<< " Y: "<<f_2[i] <<std::endl;
            file <<f_1[i]  << "," << f_2[i]  << "\n";
        }
        file.close();
        alglib::real_1d_array x ;
        alglib::real_1d_array y ;
        alglib::real_1d_array theta;
        theta.setcontent(n_samples_, t);
        x.setcontent(n_samples_, f_1);
        y.setcontent(n_samples_, f_2);
        alglib::real_1d_array x2;
        alglib::real_1d_array y2;
        alglib::spline1dinterpolant s1;
        alglib::spline1dinterpolant s2;

        // Build B-spline
        alglib::spline1dbuildcubic(theta, x, s1);
        alglib::spline1dbuildcubic(theta, y, s2);
        //Save spline contour to contours_vector 

        // Evaluate the B-spline at x=2.5
        double v1 = alglib::spline1dcalc(s1, 0.5);
        double v2 = alglib::spline1dcalc(s2, 0.5);
        std::cout << "Value at theta=0.5: "<<" X: " << v1 << " Y: "  << v2<< std::endl;
        size_t n_plotting_samples = 1000;
        std::string file_ = "/home/raphael/robolab/displaygp/config/contour_" + std::to_string(file_num) + ".csv";
        std::ofstream file_c(file_);

        if (!file_c.is_open()) 
        {
            std::cerr << "Failed to open the file." << std::endl;
        }

        file_c << "X,Y\n";

        for (size_t i = 0; i < n_plotting_samples; i++) 
        {
            double increment = static_cast<double>(i)/static_cast<double>(n_plotting_samples-1);
            double v1 = alglib::spline1dcalc(s1,increment);
            double v2 = alglib::spline1dcalc(s2, increment);
            file_c << v1  << "," << v2 << "\n";
        }

        file_c.close();
        file_num++;

        auto ptr1 = std::make_shared<alglib::spline1dinterpolant>(s1);
        auto ptr2 = std::make_shared<alglib::spline1dinterpolant>(s2);

        // Move them into the vector
        spline_contours_.push_back(std::make_pair(ptr1, ptr2));
    }
}
SplineInterpolant_ptr_pair_vec Contour::getSplineInterpolant()
{
    
    return spline_contours_;
}

bayesopt::vectord Contour::getLastSample()
{
    return bopt_model_->getData()->getLastSampleX();
}

void Contour::getInitialSamples( std::vector<double> &samples_x, std::vector<double> &samples_y )
{
     size_t n_points = bopt_model_->getData()->getNSamples();
      for (size_t i = 0; i<n_points;++i)
	    {
	      const vectord last = bopt_model_->getData()->getSampleX(i);
	      samples_x.push_back(last(0));
	      samples_y.push_back(last(1));
	    }

}

bayesopt::ProbabilityDistribution* Contour::getPredictionGaussianProcess(const bayesopt::vectord &q)
{
    return bopt_model_->getPrediction(q);
}
double Contour::evaluateGaussianProcess(const bayesopt::vectord &q)
{
    return bopt_model_->evaluateSample(q);

}

void Contour::prepareGaussianProcess()
{
    bopt_model_->initializeOptimization();
    size_t nruns = bopt_model_->getParameters()->n_iterations;
  
}
void  Contour::stepRunGaussianProcess()
{
    bopt_model_->stepOptimization();
}
double Contour::evaluateCriteriaGaussianProcess(const bayesopt::vectord &q)
{
    return bopt_model_->evaluateCriteria(q);

}
