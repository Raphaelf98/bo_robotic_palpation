#ifndef _EVALUATION_HPP_
#define _EVALUATION_HPP_
#include<boost/math/quadrature/trapezoidal.hpp>
#include "stdafx.h"
#include "interpolation.h"
#include "contour.hpp"
#include<list>
#include<iterator>
#include <CGAL/Polygon_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_rational.h>


typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                                   Point_2;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;

typedef CGAL::Polygon_with_holes_2<Kernel>                Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2>                   Pwh_list_2;

//Custom types
typedef std::vector<std::pair<Contour*, Contour*>> ContourPairs;
typedef double (*FunctionPtr)(double);

template<class Kernel, class Container>
void print_polygon (const CGAL::Polygon_2<Kernel, Container>& P)
{
  typename CGAL::Polygon_2<Kernel, Container>::Vertex_const_iterator vit;
  std::cout << "[ " << P.size() << " vertices:";
  for (vit = P.vertices_begin(); vit != P.vertices_end(); ++vit)
    std::cout << " (" << *vit << ')';
  std::cout << " ]" << std::endl;
}
template<class Kernel, class Container>

void print_polygon_with_holes(const CGAL::Polygon_with_holes_2<Kernel, Container> & pwh)
{
  if (! pwh.is_unbounded()) {
    std::cout << "{ Outer boundary = ";
    print_polygon (pwh.outer_boundary());
  } else
    std::cout << "{ Unbounded polygon." << std::endl;
  typename CGAL::Polygon_with_holes_2<Kernel,Container>::Hole_const_iterator hit;
  unsigned int k = 1;
  std::cout << " " << pwh.number_of_holes() << " holes:" << std::endl;
  for (hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit, ++k) {
    std::cout << " Hole #" << k << " = ";
    print_polygon (*hit);
  }
  std::cout << " }" << std::endl;
}
class evaluate
{

typedef std::vector<std::pair<alglib::spline1dinterpolant*, alglib::spline1dinterpolant*>> SplinePairs;
private: 
    ContourPairs contour_pairs_;
    SplinePairs spline_pairs_;
    
    /* data */
public:
    
    //Constructor to compare multiple Contour(ground truth) Contour(approximation) pairs
    evaluate(std::vector <std::pair<Contour*, Contour*>>);
    
    bool polygonizeSpline( std::pair<std::unique_ptr<alglib::spline1dinterpolant>,std::unique_ptr<alglib::spline1dinterpolant>>   &spline_pair, Polygon_2 &P,size_t num_vertices);
    
    double computeDifferenceArea(Polygon_2 &A, Polygon_2 &B, bool verbose);
    double computeJoinArea(Polygon_2 &A, Polygon_2 &B, bool verbose);
    double computeIntersectArea(Polygon_2 &A, Polygon_2 &B, bool verbose);

    double computeFalseNegative(Polygon_2 &GroundTruth, Polygon_2 &Approximation);
    double computeFalsePositive(Polygon_2 &GroundTruth, Polygon_2 &Approximation);
    double computeTruePositive(Polygon_2 &GroundTruth, Polygon_2 &Approximation);
    double computeTrueNegative(double domainArea, Polygon_2 &GroundTruth, Polygon_2 &Approximation);
    double computeArea(Polygon_2 &A);
    double computeSpecificity(double domainArea, Polygon_2 &GroundTruth, Polygon_2 &Approximation);
    double computeSensitiviy(Polygon_2 &GroundTruth, Polygon_2 &Approximation);
    ~evaluate();
};
#endif