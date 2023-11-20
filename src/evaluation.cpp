#include "evaluation.hpp"
#include <iostream>
#include <cmath>
#include <limits>
#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interpolation.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <memory>
#include "dataset.hpp"
#include "prob_distribution.hpp"
#include "fileParser.hpp"
#include "param_loader.hpp"
#include<bayesopt/parameters.hpp>

using namespace alglib;
namespace ublas = boost::numeric::ublas;

ContourPairAnalyzer::ContourPairAnalyzer(SplineInterpolant_ptr_pair &f_parametric_A, SplineInterpolant_ptr_pair &f_parametric_B, double domain_area):domain_area_(domain_area)
{   
    f_parametric_A_.first = f_parametric_A.first;
    f_parametric_A_.second = f_parametric_A.second;
    f_parametric_B_.first = f_parametric_B.first;
    f_parametric_B_.second = f_parametric_B.second;
}   
ContourPairAnalyzer::~ContourPairAnalyzer()
{
}
void ContourPairAnalyzer::analyzeContours(){
    Polygon_2 A;
    Polygon_2 B;
    polygonizeSpline(f_parametric_A_, A,1000);
  
    std::cout << "AREA OF GROUND TRUTH:  = "<< A.area()<<std::endl;
    polygonizeSpline(f_parametric_B_, B,1000);
   
    std::cout << "AREA OF APPROXIMATION:  = "<< B.area()<<std::endl;

    computeFalseNegative(A,B);
    computeFalsePositive(A,B);
    computeTrueNegative(domain_area_, A,B);
    computeTruePositive(A,B);
    computeSensitiviy(A,B);
    computeSpecificity(domain_area_,A,B);
}
bool ContourPairAnalyzer::polygonizeSpline( SplineInterpolant_ptr_pair &spline_pair, Polygon_2 &P, size_t num_vertices = 1000)
{   
    
    std::cout<<"Number of vertices: "<<num_vertices<<std::endl;
    std::shared_ptr<spline1dinterpolant> s1 = spline_pair.first;
    std::shared_ptr<spline1dinterpolant> s2= spline_pair.second;
    for (size_t i = 0; i < num_vertices; i++)
    {   
        //counter-clockwise 
        double x = spline1dcalc(*s1,  1.0-i * 1.0/((double)num_vertices));
        double y = spline1dcalc(*s2,  1.0-i * 1.0/((double)num_vertices));
        Point_2 p(x,y);
        P.push_back(p);
    }

    return 0;
}
double ContourPairAnalyzer::computeArea(Polygon_2 &A)
{
    return static_cast<double>(A.area().exact());
}

//Compute A - B
double ContourPairAnalyzer::computeDifferenceArea(Polygon_2 &A, Polygon_2 &B, bool verbose = false)
{
     if ((CGAL::do_intersect (A, B))){
        std::cout << "The two polygons intersect in their interior." << std::endl;
            
        Pwh_list_2 symmR;
        Pwh_list_2::const_iterator it;
        //Segmentation fault when A is smaller B and fully contained within B. The other way around works
        //Check if A fully contained in B 
        
        
        CGAL::difference(A, B,std::back_inserter(symmR) );
        
        
        if (verbose){
            std::cout << "Difference: " <<  std::endl;
            for (it = symmR.begin(); it != symmR.end(); ++it) 
            {
                std::cout << "--> ";
                print_polygon_with_holes(*it);    
            }

        }
        
    if (! (*symmR.begin()).is_unbounded())
     {
        double diff_area = computeArea((*symmR.begin()).outer_boundary());
        
        return diff_area;
     }
    
    }
    
    // AREA: Returns the signed area of the polygon.
    //This means that the area is positive for counter clockwise polygons and negative for clockwise polygons.
    else{
        std::cout << "The two polygons do not intersect." << std::endl;
        return 0;
    }

}

double ContourPairAnalyzer::computeJoinArea(Polygon_2 &A, Polygon_2 &B, bool verbose = false)
{
     if ((CGAL::do_intersect (A, B))){
        std::cout << "The two polygons intersect in their interior." << std::endl;
            
        Polygon_with_holes_2 unionR;

        CGAL::join(A, B,unionR);
        if (verbose){
            std::cout << "Join: " <<  std::endl;
            
            print_polygon_with_holes(unionR);    
            

        }
        
    if (! unionR.is_unbounded())
     {
        double join_area = computeArea(unionR.outer_boundary());
        
        return join_area;
     }
    
    }
    
    // AREA: Returns the signed area of the polygon.
    //This means that the area is positive for counter clockwise polygons and negative for clockwise polygons.
    else{
        std::cout << "The two polygons do not intersect." << std::endl;
        return 0;
    }

}
double ContourPairAnalyzer::computeIntersectArea(Polygon_2 &A, Polygon_2 &B, bool verbose = false)
{  
    if ((CGAL::do_intersect (A, B)))
    {
        Pwh_list_2                  intR;
        Pwh_list_2::const_iterator  it;
        CGAL::intersection (A, B, std::back_inserter(intR));
        if (verbose)
        {
            std::cout << "The intersection:" << std::endl;
            for (it = intR.begin(); it != intR.end(); ++it) 
            {
                std::cout << "--> ";
                print_polygon_with_holes (*it);    
            }
        }
        if (! (*intR.begin()).is_unbounded())
        {
            double intersect_area = computeArea((*intR.begin()).outer_boundary());

            return intersect_area;
        }
    }
    else
    {
        std::cout << "The two polygons do not intersect." << std::endl;
        return 0;
    }
   
  
}
double ContourPairAnalyzer::computeFalseNegative(Polygon_2 &GroundTruth, Polygon_2 &Approximation)
{
    double join_area = computeJoinArea(GroundTruth,Approximation, false);
    double area_A = computeArea(GroundTruth);
    double area_B = computeArea(Approximation);
        //check if A contained in B or vice versa
    if (join_area == area_B )
    {
            //TODO: THINK ABOUT VERBOSITY
            std::cout << "GroundTruth is fully contained in Approximation and Joined area exactly equal to: "<< join_area<<  std::endl;
            return 0 ;
    }
    if (join_area == area_A )
    {
            //TODO: THINK ABOUT VERBOSITY
            std::cout << "Approximation is fully contained in GroundTruth and Joined area exactly equal to: "<< join_area<<  std::endl;
            return area_A-area_B ;
    }
    double area = computeDifferenceArea(Approximation,GroundTruth);
    std::cout << "FalseNegative Area: " << area << std::endl;
    return area;
}

double ContourPairAnalyzer::computeFalsePositive(Polygon_2 &GroundTruth, Polygon_2 &Approximation)
{
    double join_area = computeJoinArea(GroundTruth,Approximation, false);
    double area_A = computeArea(GroundTruth);
    double area_B = computeArea(Approximation);
        //check if A contained in B or vice versa
    if (join_area == area_B )
    {
            //TODO: THINK ABOUT VERBOSITY
            std::cout << "GroundTruth is fully contained in Approximation and Joined area exactly equal to: "<< join_area<<  std::endl;
            return area_B-area_A ;
    }
    if (join_area == area_A )
    {
            //TODO: THINK ABOUT VERBOSITY
            std::cout << "Approximation is fully contained in GroundTruth and Joined area exactly equal to: "<< join_area<<  std::endl;
            return 0.0;
    }
    double area = computeDifferenceArea(GroundTruth,Approximation);
    std::cout << "FalsePositive Area: " << area << std::endl;
    return area;
}
double ContourPairAnalyzer::computeTruePositive(Polygon_2 &GroundTruth, Polygon_2 &Approximation)
{

    double area = computeIntersectArea(GroundTruth,Approximation);
    std::cout << "Intersect Area: " << area << std::endl;
    return area;
}
double ContourPairAnalyzer::computeTrueNegative(double domainArea, Polygon_2 &GroundTruth, Polygon_2 &Approximation)
{
    double area = computeJoinArea(GroundTruth,Approximation);
    std::cout << "Join Area: " << area << std::endl;
    return domainArea - area;
    
}
double ContourPairAnalyzer::computeSpecificity(double domainArea, Polygon_2 &GroundTruth, Polygon_2 &Approximation)
{
    double true_negative = computeTrueNegative(domainArea, GroundTruth, Approximation);
    double specificity = true_negative/(true_negative+computeFalsePositive(GroundTruth,Approximation));
     std::cout<<"Specificity: "<< specificity<<std::endl;
    return specificity;
}
double ContourPairAnalyzer::computeSensitiviy(Polygon_2 &GroundTruth, Polygon_2 &Approximation)
{
    double true_positive = computeTruePositive(GroundTruth,Approximation);
    double sensitivity =  true_positive/(true_positive + computeFalseNegative(GroundTruth,Approximation));
    std::cout<<"Sensitivity: "<< sensitivity<<std::endl;
    return sensitivity;
}


double fx1(double x)
{
    return sin(x);
}
double fy1(double x)
{
    return cos(x);
}

double fx2(double x)
{
    return   0.5+0.1*sin(x);
}
double fy2(double x)
{
    return 0.5 +  0.1*cos(x);
}
double fx3(double x)
{
    return  0.9*sin(x);
}
double fy3(double x)
{
    return 0.9*cos(x);
}
std::pair<std::unique_ptr<alglib::spline1dinterpolant>,std::unique_ptr<alglib::spline1dinterpolant>> f_param(FunctionPtr &f_x, FunctionPtr &f_y, int n_samples_ = 10)
{
    double t[n_samples_];
    std::cout<<"sample size: "<<sizeof(t)/sizeof(double)<<std::endl;
    double f_1[n_samples_];
    double f_2[n_samples_];
    
    for (size_t i = 0; i < n_samples_; i++)
    {
        f_1[i] =f_x(1-(double) i/((double) n_samples_-1));
        

        f_2[i] =f_y(1-(double) i/((double) n_samples_-1));
        std::cout<< "VALUE" << f_1[i]<<", "<< f_2[i] <<"at: "<< (double) i/((double) n_samples_-1)<<std::endl;
        t[i] = i * 1.0/((double)n_samples_-1);
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
    
    return std::make_pair(std::make_unique<alglib::spline1dinterpolant>(s1),std::make_unique<alglib::spline1dinterpolant>(s2));
}



enum ShapeType {
    SHAPE_CIRCLE,
    SHAPE_TRIANGLE,
    SHAPE_RECTANGLE,
    SHAPE_TWOCIRCLES,
    SHAPE_UNKNOWN // for unrecognized strings
};
ShapeType getShapeType(const std::string& shape) {
    if (shape == "Circle") return SHAPE_CIRCLE;
    if (shape == "Triangle") return SHAPE_TRIANGLE;
    if (shape == "Rectangle") return SHAPE_RECTANGLE;
    if (shape == "TwoCircles") return SHAPE_TWOCIRCLES;
    return SHAPE_UNKNOWN;
}

int main(int argc, char *argv[]){
    //COMPUTE CONTOUR
    bayesopt::Parameters par;
 
  
  if (bayesopt::utils::ParamLoader::load("/home/raphael/robolab/displaygp/config/bo_parameters.txt", par))
  {
      std::cout << "Found bo_parameters.txt" << std::endl;
      
  }
  else
  {
    par = initialize_parameters_to_default();
    par.n_iterations = 60;
    par.n_init_samples = 10;
    par.crit_name = "cEI";
    par.epsilon = 3;
    par.random_seed = 10;
    par.init_method = 3;
    par.mean.name = "mZero";
    par.force_jump = 0;
    par.kernel.name = "kSEARD";
    par.kernel.hp_mean[0] = 0.08;
    par.kernel.hp_std[0] = 1.0;
    par.n_inner_iterations = 500;
    par.verbose_level = 1;
    
    }
    /*if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <Shape>\n";
        return 1;
    }*/
    
    //std::string arg = argv[1];
    std::string arg = "Triangle";
    ShapeType type = getShapeType(arg);
   
    std::vector< std::pair <std::function<double (const double&)> ,std::function<double (const double&)>> > groundTruths;
    std::unique_ptr<Shape> shape;

    switch (type) {

        case SHAPE_CIRCLE:
            std::cout << "Running Experiment on Circle" << std::endl;
            shape = std::make_unique<SmoothCircle>(par); 
            
            break;

        case SHAPE_TRIANGLE:
            std::cout << "Running Experiment on Triangle" << std::endl;
            shape = std::make_unique<Triangle>(par); 
            groundTruths.push_back(std::make_pair(shape->f_x(), shape->f_y()));
            break;

        case SHAPE_RECTANGLE:
            std::cout << "Running Experiment on Rectangle" << std::endl;
            shape = std::make_unique<Rectangle>(par); 
            break;

        case SHAPE_TWOCIRCLES:
            std::cout << "Running Experiment on Two Circles" << std::endl;
            shape = std::make_unique<TwoCircles>(par, 0.1,0.15,0.2,0.7,0.8,0.3,0.1); 
            break;
            
        default:
            std::cout << "Unknown Shape: " << arg << std::endl;
    }
    
    

    
    Contour contour(shape.get(),100);
    
    //run bayesian optimization
    contour.runGaussianProcess();
    contour.computeCluster();
    //dummy data
    
    contour.exploreContour();
    contour.approximateContour();
    
    SplineInterpolant_ptr_pair_vec spline_vec = contour.getSplineInterpolant();
    SplineInterpolant_ptr_pair  f_Groundtruth = static_cast<SplineInterpolant_ptr_pair>(f_param(groundTruths[0].first, groundTruths[0].second,1000));
    //SplineInterpolant_ptr_pair  test_spline = static_cast<SplineInterpolant_ptr_pair>(f_param(fx3,fy3,1000));
    /*TODO: Contour Analyzer works at the moment only for single Conotur
            -spline_vec may contour multiple countours -> correct ground truth must be assigned*/
    ContourPairAnalyzer analyzer(f_Groundtruth,spline_vec[0], 1.0);
    //ContourPairAnalyzer analyzer(f_Groundtruth,test_spline, 1.0);
    analyzer.analyzeContours();
    return 0;
}