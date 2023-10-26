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

using namespace alglib;
namespace ublas = boost::numeric::ublas;

evaluate::evaluate(ContourPairs contour_pairs):contour_pairs_(contour_pairs)
{   
    for (auto contourpair : contour_pairs_)
    {
        Contour *c1 = contourpair.first;
        Contour *c2 = contourpair.second;
        spline_pairs_.push_back(std::make_pair(c1->getSplineInterpolant()[0],c2->getSplineInterpolant()[1]));
    }
  


}   
evaluate::~evaluate()
{
}
bool evaluate::polygonizeSpline( std::pair<std::unique_ptr<spline1dinterpolant>,std::unique_ptr<spline1dinterpolant>>  &spline_pair, Polygon_2 &P, size_t num_vertices = 1000)
{   
    
    std::cout<<"Number of vertices: "<<num_vertices<<std::endl;
    std::unique_ptr<spline1dinterpolant> s1 = std::move(spline_pair.first);
    std::unique_ptr<spline1dinterpolant> s2= std::move(spline_pair.second);
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
double evaluate::computeArea(Polygon_2 &A)
{
    return static_cast<double>(A.area().exact());
}
double evaluate::computeDifferenceArea(Polygon_2 &A, Polygon_2 &B, bool verbose = false)
{
     if ((CGAL::do_intersect (A, B))){
        std::cout << "The two polygons intersect in their interior." << std::endl;
            
        Pwh_list_2 symmR;
        Pwh_list_2::const_iterator it;
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

double evaluate::computeJoinArea(Polygon_2 &A, Polygon_2 &B, bool verbose = false)
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
double  evaluate::computeIntersectArea(Polygon_2 &A, Polygon_2 &B, bool verbose = false)
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
double evaluate::computeFalseNegative(Polygon_2 &GroundTruth, Polygon_2 &Approximation)
{   
    double area = computeDifferenceArea(Approximation,GroundTruth);
    std::cout << "FalseNegative Area: " << area << std::endl;
    return area;
}

double evaluate::computeFalsePositive(Polygon_2 &GroundTruth, Polygon_2 &Approximation)
{
    double area = computeDifferenceArea(GroundTruth,Approximation);
    std::cout << "FalsePositive Area: " << area << std::endl;
    return area;
}
double evaluate::computeTruePositive(Polygon_2 &GroundTruth, Polygon_2 &Approximation)
{

    double area = computeIntersectArea(GroundTruth,Approximation);
    std::cout << "Intersect Area: " << area << std::endl;
    return area;
}
double evaluate::computeTrueNegative(double domainArea, Polygon_2 &GroundTruth, Polygon_2 &Approximation)
{
    double area = computeJoinArea(GroundTruth,Approximation);
    std::cout << "Join Area: " << area << std::endl;
    return domainArea - area;
    
}
double evaluate::computeSpecificity(double domainArea, Polygon_2 &GroundTruth, Polygon_2 &Approximation)
{
    double true_negative = computeTrueNegative(domainArea, GroundTruth, Approximation);
    double specificity = true_negative/(true_negative+computeFalsePositive(GroundTruth,Approximation));
     std::cout<<"Specificity: "<< specificity<<std::endl;
    return specificity;
}
double evaluate::computeSensitiviy(Polygon_2 &GroundTruth, Polygon_2 &Approximation)
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
    return 0.5+sin(x);
}
double fy2(double x)
{
    return cos(x);
}
std::pair<std::unique_ptr<spline1dinterpolant>,std::unique_ptr<spline1dinterpolant>> f_param(FunctionPtr f_x, FunctionPtr f_y, int n_samples_ = 10)
{
    double t[n_samples_];
    std::cout<<"sample size: "<<sizeof(t)/sizeof(double)<<std::endl;
    double f_1[n_samples_];
    double f_2[n_samples_];
    for (size_t i = 0; i < n_samples_; i++)
    {
        f_1[i] =f_x(2*M_PI*i/(n_samples_-1));
        f_2[i] =f_y(2*M_PI*i/(n_samples_-1));;
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

int main(){
    ContourPairs contourpairs;
    evaluate evaluation_(contourpairs);

    Polygon_2 A;
    std::pair<std::unique_ptr<spline1dinterpolant>,std::unique_ptr<spline1dinterpolant>>  f_A = f_param(fx1,fy1,1000);
    evaluation_.polygonizeSpline(f_A, A);
    std::cout << "AREA P:  = "<< A.area()<<std::endl;

    Polygon_2 B;
    std::pair<std::unique_ptr<spline1dinterpolant>,std::unique_ptr<spline1dinterpolant>>  f_B = f_param(fx2,fy2,1000);
    evaluation_.polygonizeSpline(f_B, B);
    std::cout << "AREA P:  = "<< B.area()<<std::endl;

    evaluation_.computeFalseNegative(A,B);
    evaluation_.computeFalsePositive(A,B);
    evaluation_.computeTrueNegative(9.0, A,B);
    evaluation_.computeTruePositive(A,B);

    evaluation_.computeSensitiviy(A,B);
    evaluation_.computeSpecificity(9.0,A,B);

    return 0;
}