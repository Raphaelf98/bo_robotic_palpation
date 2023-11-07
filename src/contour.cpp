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

Contour::Contour(bayesopt::BayesOptBase* bopt_model, size_t n_exploration_directions):
        bopt_model_(bopt_model),
        n_exploration_directions_(n_exploration_directions),
        c_points_(100),
        bandwidth_(0.05),
        samples_(c_points_),
        n_directions_(8),
        stiffness_threshold_(0.01),
        lim_steps_(1000),
        n_samples_(n_directions_+1),
        c_(c_points_,std::vector<double>(c_points_))
{

}

Contour::~Contour()
{
}
void Contour::runGaussianProcess(){
    size_t state_ii = 0;
    
    bopt_model_->initializeOptimization();
    size_t nruns = bopt_model_->getParameters()->n_iterations;
    std::vector<double> x,y;
    vectord q(2);
    x = linspace(0,1,c_points_);
    y = linspace(0,1,c_points_);
    while(state_ii < nruns)
    {
        ++state_ii;
        bopt_model_->stepOptimization();

    }
    for(size_t i=0; i<c_points_; ++i)
	{
        
          for (size_t j=0; j<c_points_; ++j)
          {
            q(1) = y[i]; q(0) = x[j];
            bayesopt::ProbabilityDistribution* pd = bopt_model_->getPrediction(q);
            c_[i][j] = pd->getMean();   
            
          }

	      
    }

}
void Contour::computeCluster()
{   
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

void Contour::exploreContour(){
    /*
    TODO: 
    */
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
                    if (abs(new_stiffness-old_stiffness) >= stiffness_threshold_)
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
void Contour::approximateContour()
{
    int file_num = 1;
    for(std::vector<Point> &contour : contours_)
    {
    double t[n_samples_];
    std::cout<<"sample size: "<<sizeof(t)/sizeof(double)<<std::endl;
    double f_1[n_samples_];
    double f_2[n_samples_];
    for (size_t i = 0; i < n_samples_; i++)
    {
        f_1[i] =contour[i].x;
        f_2[i] =contour[i].y;
        t[i] = i * 1.0/((double)n_samples_-1);
        std::cout<<"sample # "<<t[i] << " at X: " <<f_1[i]<< " Y: "<<f_2[i] <<std::endl;
    }

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
    std::ofstream file(file_);

    if (!file.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
    }

    file << "X,Y\n";

    for (size_t i = 0; i < n_plotting_samples; i++) {
        double increment = static_cast<double>(i)/static_cast<double>(n_plotting_samples-1);
        double v1 = alglib::spline1dcalc(s1,increment);
        double v2 = alglib::spline1dcalc(s2, increment);
        file << v1  << "," << v2 << "\n";
    }

    file.close();
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

