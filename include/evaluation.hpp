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
#include <functional>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                                   Point_2;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;

typedef CGAL::Polygon_with_holes_2<Kernel>                Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2>                   Pwh_list_2;

//Custom types
typedef std::vector<std::pair<Contour*, Contour*>> ContourPairs;
typedef std::function<double(const double&)> FunctionPtr;

/**
 * @brief Plot information on nodes and vertices of CGAL Polygon_2.
 *
 * @param[in] P Polygon object. 
 *
 */

template<class Kernel, class Container>
void print_polygon (const CGAL::Polygon_2<Kernel, Container>& P)
{
  typename CGAL::Polygon_2<Kernel, Container>::Vertex_const_iterator vit;
  std::cout << "[ " << P.size() << " vertices:";
  for (vit = P.vertices_begin(); vit != P.vertices_end(); ++vit)
    std::cout << " (" << *vit << ')';
  std::cout << " ]" << std::endl;
}


/**
 * @brief Plot information on nodes and vertices of CGAL Polygon_with_holes_2.
 *
 * @param[in] P Polygon object. 
 *
 */
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

/**
* @brief This class analyses a pair of contours by polygonizing parametric (defined through spline functions) contours.
* Similarity of Contours is determined by scalars:
*   1. Specificity 
*   2. Sensitivity
* More information can be obtained in the IEEE paper: Fast Localization and Segmentation of Tissue
* Abnormalities by Autonomous Robotic Palpation by Youcan Yan and Jia Pan (2021)
*/
class ContourPairAnalyzer
{
private:
    SplineInterpolant_ptr_pair f_parametric_A_,f_parametric_B_; //parametric Spline functions of Contour A and B
    double domain_area_;
    std::string experiment_path_;
    
public:
    
    //Constructor to compare multiple Contour(ground truth) Contour(approximation) pairs
     /**
    * @brief ContourPairAnalyzer class constructor. Takes two contours represented as spline interpolant pair (parametric function representation) 
    *
    * @param[in] f_parametric_A Ground truth. 
    * @param[in] f_parametric_B Tumor contour approximated by algo.
    * @param[in] domain_area Area of ground truth domain (1.0x1.0).
    * @param[in] experiment_path_ Current working directory path as string.
    * 
    */
    ContourPairAnalyzer(SplineInterpolant_ptr_pair &f_parametric_A, SplineInterpolant_ptr_pair &f_parametric_B, double domain_area, std::string experiment_path_);
    /**
    * @brief Compare contours A and B and compute specificity and sensitivity. Results are logged and stored to file.
    *
    * @param[in] contour_number number of contour written to file (Needed for the case that multiple ground truths and multiple tumor contours are compared).
    * 
    * 
    */ 
    void analyzeContours(const size_t contour_number);
    /**
    * @brief Takes spline object and poligonizes it. 
    *
    * @param[in] spline_pair Parametric contour representation.
    * @param[in] P CGAL polygon object that stores the polygonized spline.
    * @param[in] num_vertices Number of points at which spline is sampled.
    * 
    * 
    */
    bool polygonizeSpline( SplineInterpolant_ptr_pair &spline_pair, Polygon_2 &P,size_t num_vertices);

    /**
    * @brief Compute area of difference A - B. 
    *
    * @param[in] A Polygon A.
    * @param[in] B Polygon B.
    * @param[in] verbose Print results.
    * 
    *  @retval floating point difference area.
    */
    double computeDifferenceArea(Polygon_2 &A, Polygon_2 &B, bool verbose);
    /**
    * @brief Computes area of union of two polygons A and B. 
    *
    * @param[in] A Polygon A.
    * @param[in] B Polygon B.
    * @param[in] verbose Print results.
    * 
    * @retval floating point join area.
    */
    double computeJoinArea(Polygon_2 &A, Polygon_2 &B, bool verbose);
    /**
    * @brief Computes area of intersection of two polygons A and B. 
    *
    * @param[in] A Polygon A.
    * @param[in] B Polygon B.
    * @param[in] verbose Print results.
    * 
    * @retval floating point intersect area.
    */
    double computeIntersectArea(Polygon_2 &A, Polygon_2 &B, bool verbose);
     /**
    * @brief Computes FalseNegative according to IEEE paper: Fast Localization and Segmentation of Tissue
    * Abnormalities by Autonomous Robotic Palpation by Youcan Yan and Jia Pan (2021) 
    *
    * @param[in] GroundTruth 
    * @param[in] Approximation
    * 
    * @retval floating point false-negative area.
    */
    double computeFalseNegative(Polygon_2 &GroundTruth, Polygon_2 &Approximation);
     /**
    * @brief Computes FalsePositive according to IEEE paper: Fast Localization and Segmentation of Tissue
    * Abnormalities by Autonomous Robotic Palpation by Youcan Yan and Jia Pan (2021) 
    *
    * @param[in] GroundTruth 
    * @param[in] Approximation
    * 
    * @retval floating point false-positive area.
    */
    double computeFalsePositive(Polygon_2 &GroundTruth, Polygon_2 &Approximation);
    /**
    * @brief Computes TruePositive according to IEEE paper: Fast Localization and Segmentation of Tissue
    * Abnormalities by Autonomous Robotic Palpation by Youcan Yan and Jia Pan (2021) 
    *
    * @param[in] GroundTruth 
    * @param[in] Approximation
    * 
    * @retval floating point true-positive area.
    */
    double computeTruePositive(Polygon_2 &GroundTruth, Polygon_2 &Approximation);
    /**
    * @brief Computes TrueNegative according to IEEE paper: Fast Localization and Segmentation of Tissue
    * Abnormalities by Autonomous Robotic Palpation by Youcan Yan and Jia Pan (2021) 
    * @param[in] domainArea  
    * @param[in] GroundTruth 
    * @param[in] Approximation
    * 
    * @retval floating point true-negative area.
    */
    double computeTrueNegative(double domainArea, Polygon_2 &GroundTruth, Polygon_2 &Approximation);
    /**
    * @brief Get area of polygon.
    * @param[in] A    
    * 
    * @retval floating point area.
    */
    double computeArea(Polygon_2 &A);
    /**
    * @brief Computes Specificity according to IEEE paper: Fast Localization and Segmentation of Tissue
    * Abnormalities by Autonomous Robotic Palpation by Youcan Yan and Jia Pan (2021) 
    * @param[in] domainArea  
    * @param[in] GroundTruth 
    * @param[in] Approximation
    * 
    * @retval floating point Specificity area.
    */
    double computeSpecificity(double domainArea, Polygon_2 &GroundTruth, Polygon_2 &Approximation);
    /**
    * @brief Computes Sensitivity according to IEEE paper: Fast Localization and Segmentation of Tissue
    * Abnormalities by Autonomous Robotic Palpation by Youcan Yan and Jia Pan (2021) 

    * @param[in] GroundTruth 
    * @param[in] Approximation
    * 
    * @retval floating point Sensitivity area.
    */
    double computeSensitivity(Polygon_2 &GroundTruth, Polygon_2 &Approximation);
    ~ContourPairAnalyzer();
};
#endif