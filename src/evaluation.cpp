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
#include "fileparser.hpp"
#include "param_loader.hpp"
#include<bayesopt/parameters.hpp>
#include<fstream>

using namespace alglib;
namespace ublas = boost::numeric::ublas;
//Constructor, receives and initializes two prametric spline functions (contours)
ContourPairAnalyzer::ContourPairAnalyzer(SplineInterpolant_ptr_pair &f_parametric_A, SplineInterpolant_ptr_pair &f_parametric_B, double domain_area, std::string experiment_path):domain_area_(domain_area), experiment_path_(experiment_path)
{   
    f_parametric_A_.first = f_parametric_A.first;
    f_parametric_A_.second = f_parametric_A.second;
    f_parametric_B_.first = f_parametric_B.first;
    f_parametric_B_.second = f_parametric_B.second;
}   
ContourPairAnalyzer::~ContourPairAnalyzer()
{
}
/*
Compute Sensitivity, Specificity of Contour A and B and store results to file
*/
void ContourPairAnalyzer::analyzeContours(const size_t contour_number){
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
    double sens = computeSensitivity(A,B);
    double spec = computeSpecificity(domain_area_,A,B);
    
    std::string posterior_path = FILE_METRICS;
    posterior_path = experiment_path_+posterior_path;
    std::cout<<"posterior_path:   " << posterior_path<<std::endl;
    saveMetricsToFile(contour_number, spec, sens ,posterior_path);
}
/*
Switch from parametric spline function representation to Polygon representation of contour
*/
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
/*
Computes Area of Polygon
*/
double ContourPairAnalyzer::computeArea(Polygon_2 &A)
{
    return static_cast<double>(A.area().exact());
}

/*
Compute area of difference A - B
*/
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
     else{
        return 0.0;
     }
    
    }
    
    // AREA: Returns the signed area of the polygon.
    //This means that the area is positive for counter clockwise polygons and negative for clockwise polygons.
    else{
        std::cout << "The two polygons do not intersect." << std::endl;
        return 0.0;
    }

}
/*
Computes area of union of two polygons
*/
double ContourPairAnalyzer::computeJoinArea(Polygon_2 &A, Polygon_2 &B, bool verbose = false)
{
    if ((CGAL::do_intersect (A, B)))
    {
        std::cout << "The two polygons intersect in their interior." << std::endl;
            
        Polygon_with_holes_2 unionR;

        CGAL::join(A, B,unionR);
        if (verbose)
        {
            std::cout << "Join: " <<  std::endl;
            
            print_polygon_with_holes(unionR);    
            

        }
        
        if (! unionR.is_unbounded())
        {
           double join_area = computeArea(unionR.outer_boundary());
           
           return join_area;
        }
        else
        {
            return 0.0;
        }
    
    }
    
    // AREA: Returns the signed area of the polygon.
    //This means that the area is positive for counter clockwise polygons and negative for clockwise polygons.
    else{
        std::cout << "The two polygons do not intersect." << std::endl;
        return 0.0;
    }

}
/*
Computes area of intersection of two polygons
*/
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
        else{
            return 0.0;
        }
    }
    else
    {
        std::cout << "The two polygons do not intersect." << std::endl;
        return 0.0;
    }
   
  
}
/*
Computes FalseNegative according to IEEE paper: Fast Localization and Segmentation of Tissue
Abnormalities by Autonomous Robotic Palpation by Youcan Yan and Jia Pan (2021)
*/
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
            return 0.0 ;
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
/*
Computes FalsePositive according to IEEE paper: Fast Localization and Segmentation of Tissue
Abnormalities by Autonomous Robotic Palpation by Youcan Yan and Jia Pan (2021)
*/
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
/*
Computes TruePositive according to IEEE paper: Fast Localization and Segmentation of Tissue
Abnormalities by Autonomous Robotic Palpation by Youcan Yan and Jia Pan (2021)
*/
double ContourPairAnalyzer::computeTruePositive(Polygon_2 &GroundTruth, Polygon_2 &Approximation)
{

    double area = computeIntersectArea(GroundTruth,Approximation);
    std::cout << "Intersect Area: " << area << std::endl;
    return area;
}
/*
Computes TrueNegative according to IEEE paper: Fast Localization and Segmentation of Tissue
Abnormalities by Autonomous Robotic Palpation by Youcan Yan and Jia Pan (2021)
*/
double ContourPairAnalyzer::computeTrueNegative(double domainArea, Polygon_2 &GroundTruth, Polygon_2 &Approximation)
{
    double area = computeJoinArea(GroundTruth,Approximation);
    std::cout << "Join Area: " << area << std::endl;
    return domainArea - area;
    
}
/*
Computes Specificity according to IEEE paper: Fast Localization and Segmentation of Tissue
Abnormalities by Autonomous Robotic Palpation by Youcan Yan and Jia Pan (2021)
*/
double ContourPairAnalyzer::computeSpecificity(double domainArea, Polygon_2 &GroundTruth, Polygon_2 &Approximation)
{
    double true_negative = computeTrueNegative(domainArea, GroundTruth, Approximation);
    double specificity = true_negative/(true_negative+computeFalsePositive(GroundTruth,Approximation));
     std::cout<<"Specificity: "<< specificity<<std::endl;
    return specificity;
}
/*
Computes Sensitivity according to IEEE paper: Fast Localization and Segmentation of Tissue
Abnormalities by Autonomous Robotic Palpation by Youcan Yan and Jia Pan (2021)
*/
double ContourPairAnalyzer::computeSensitivity(Polygon_2 &GroundTruth, Polygon_2 &Approximation)
{
    double true_positive = computeTruePositive(GroundTruth,Approximation);
    double sensitivity =  true_positive/(true_positive + computeFalseNegative(GroundTruth,Approximation));
    std::cout<<"Sensitivity: "<< sensitivity<<std::endl;
    return sensitivity;
}


/*
Approximates two dimensional parametric function defined analytical expressions through a pair of splines. 

*/
std::pair<std::unique_ptr<alglib::spline1dinterpolant>,std::unique_ptr<alglib::spline1dinterpolant>> f_param(std::ofstream &file,std::string file_path, FunctionPtr &f_x, FunctionPtr &f_y, int n_samples_ = 1000)
{
    double t[n_samples_];
    //std::cout<<"sample size: "<<sizeof(t)/sizeof(double)<<std::endl;
    double f_1[n_samples_];
    double f_2[n_samples_];
    
    for (size_t i = 0; i < n_samples_; i++)
    {
        f_1[i] =f_x(1-(double) i/((double) n_samples_-1));
       
        f_2[i] =f_y(1-(double) i/((double) n_samples_-1));
        //std::cout<< "VALUE" << f_1[i]<<", "<< f_2[i] <<"at: "<< (double) i/((double) n_samples_-1)<<std::endl;
        t[i] = i * 1.0/((double)n_samples_-1);
    }
    
    // Write data to file

    file.open(file_path, std::ofstream::app);
    for (int i = 0; i < n_samples_; ++i) 
    {
        file << f_1[i] << "," << f_2[i] << "\n";
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
    
    return std::make_pair(std::make_unique<alglib::spline1dinterpolant>(s1),std::make_unique<alglib::spline1dinterpolant>(s2));
}

/*
enum to carry tumor model shapes
*/
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





int main(int argc, char *argv[])
{

  
  bayesopt::Parameters par;
  TumorModelParameters model_parameters;
  ContourParamters contour_parameters;
  //Either load bayesian optimization parameters from file or set them manually
  //TODO adjust file parsing
  std::string config_path = generateFilePath(CONFIG_PATH,"");
  std::string bo_parameters = FILE_BO_PARAMETERS;
  bo_parameters = config_path+bo_parameters;
  if (bayesopt::utils::ParamLoader::load(bo_parameters, par))
  {
      std::cout << "Found "+ bo_parameters << std::endl;
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
     std::cout << "Using default paramters" << std::endl;
  }
  std::string contour_params_path = FILE_CONTOUR_PARAMETERS;
  contour_params_path = config_path+contour_params_path;
  if (contour_parameters.loadContourParameters(contour_params_path, contour_parameters))
  {
      std::cout << "Found contour_parameters.txt" << std::endl;
      contour_parameters.PrintParameters();
      
  }
  else
  {
    std::cout<<"Could not load contour_parameters.txt"<< std::endl;
  }
  std::string tumor_params_path = FILE_TUMOR_MODEL_PARAMETERS;
  tumor_params_path = config_path+tumor_params_path;
  if (model_parameters.loadModelParameters(tumor_params_path, model_parameters))
  {
      std::cout << "Found tumor_model_parameters.txt" << std::endl;
       model_parameters.printParameters();
  }
  else
  {
    std::cout<<"Could not load tumor_model_parameters"<< std::endl;
  }

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <Shape>\n";
        return 1;
    }

    std::string arg = argv[1];
    std::string experiment_path = generateFilePath(DATA_PATH,"");
    std::string dirName = createShapeDirectory(experiment_path,arg);
    std::string log_dir = LOG_PATH;
    log_dir = dirName + log_dir;
    createOrOverwriteDirectory(log_dir);
    std::string results_dir = RESULTS_PATH;
    results_dir = dirName + results_dir;
    createOrOverwriteDirectory(results_dir);
    //std::string arg = "TwoCircles";


    std::string contour_params_file = FILE_CONTOUR_PARAMETERS;
    
    std::string model_params_file = FILE_TUMOR_MODEL_PARAMETERS;
    std::string params_dir = PARAMS_PATH;
    params_dir = dirName + params_dir;
    createOrOverwriteDirectory(params_dir);
    std::cout<<"params_dir: "<<params_dir<<std::endl;
    //model_params_file = log_dir + model_params_file;
    std::string bo_paramters = FILE_BO_PARAMETERS;
    bo_parameters = params_dir+bo_paramters;
    bayesopt::utils::ParamLoader::save(bo_parameters,par);
    copyFileToDirectory(contour_params_path, params_dir, contour_params_file);
    copyFileToDirectory(tumor_params_path, params_dir, model_params_file);

    ShapeType type = getShapeType(arg);
   
    
    std::unique_ptr<Shape> shape;
    std::vector< std::pair <std::function<double (const double&)> ,std::function<double (const double&)>> > groundTruths;
    switch (type) {

        case SHAPE_CIRCLE:
            std::cout << "Running Experiment on Circle" << std::endl;
            shape = std::make_unique<SmoothCircle>(par,model_parameters.circle_low,model_parameters.circle_high,model_parameters.circle_radius,
                                                    model_parameters.circle_x_trans,model_parameters.circle_y_trans,model_parameters.circle_epsilon,model_parameters.circle_noise); 
            groundTruths.push_back(std::make_pair(shape->f_x(), shape->f_y()));
            break;

        case SHAPE_TRIANGLE:
            std::cout << "Running Experiment on Triangle" << std::endl;
            shape = std::make_unique<Triangle>(par,model_parameters.triangle_low,model_parameters.triangle_high,model_parameters.triangle_radius,
                                                    model_parameters.triangle_x_trans,model_parameters.triangle_y_trans,model_parameters.triangle_epsilon, model_parameters.triangle_noise); 
            groundTruths.push_back(std::make_pair(shape->f_x(), shape->f_y()));
            break;

        case SHAPE_RECTANGLE:
            std::cout << "Running Experiment on Rectangle" << std::endl;
            shape = std::make_unique<Rectangle>(par,model_parameters.rectangle_low,model_parameters.rectangle_high,model_parameters.rectangle_radius,
                                                    model_parameters.rectangle_x_trans,model_parameters.rectangle_y_trans,model_parameters.rectangle_epsilon, model_parameters.rectangle_noise); 
            groundTruths.push_back(std::make_pair(shape->f_x(), shape->f_y()));
            break;

        case SHAPE_TWOCIRCLES:
            std::cout << "Running Experiment on Two Circles" << std::endl;
            shape = std::make_unique<TwoCircles>(par,model_parameters.two_circles_low,model_parameters.two_circles_high,model_parameters.two_circles_radius_1,model_parameters.two_circles_radius_2,
                                                    model_parameters.two_circles_x_trans_1,model_parameters.two_circles_x_trans_2,  model_parameters.two_circles_y_trans_1, model_parameters.two_circles_y_trans_2,
                                                            model_parameters.rectangle_epsilon, model_parameters.two_circles_noise);
            groundTruths.push_back(std::make_pair(shape->f_x(), shape->f_y()));
            groundTruths.push_back(std::make_pair(shape->f_x(), shape->f_y()));
            break;
            
        default:
            std::cout << "Unknown Shape: " << arg << std::endl;
    }
    
    shape->saveGroundTruth(contour_parameters.c_points,log_dir);
    Contour contour(shape.get(),contour_parameters, dirName);
    //run bayesian optimization
    contour.runGaussianProcess();
    contour.computeCluster();
    contour.exploreContour();
    contour.approximateContour();
    
    SplineInterpolant_ptr_pair_vec spline_vec = contour.getSplineInterpolant();
    
    std::string path = generateExperimentFilePath(dirName, LOG_PATH,FILE_EVAL_GROUND_TRUTH);
    std::ofstream file(path);
    file << "X,Y\n";
    file.close();
    std::vector<int> idx_list ={0};
    //selects the correct index based on list of centroids provided in ms_centroids.csv for list iterating over over gorund truth contours
    if (groundTruths.size() > 1)

    {   idx_list.clear();
        std::string centroids_file = FILE_MS_CENTROIDS;
        centroids_file = log_dir + centroids_file;
        std::vector<std::pair<double, double>> centroids = readCoordinatesFromCSV(centroids_file);
        std::pair<double, double> ground_truth_centroid_1= std::make_pair(model_parameters.two_circles_x_trans_1,model_parameters.two_circles_y_trans_1);
        std::pair<double, double> ground_truth_centroid_2= std::make_pair(model_parameters.two_circles_x_trans_2,model_parameters.two_circles_y_trans_2);
        std::vector<std::pair<double, double>> gts = {ground_truth_centroid_1,ground_truth_centroid_2};
        
        int idx;
        double dist = 10e6;
        
        for(size_t i = 0; i < centroids.size(); i++)
        {   
            for(size_t j = 0; j < gts.size(); j++)
            {   
                double norm = sqrt((centroids[i].first-gts[j].first)* (centroids[i].first-gts[j].first) + (centroids[i].second-gts[j].second) * (centroids[i].second-gts[j].second) );
                if(norm < dist){
                    idx = j;
                }
                dist = norm;
            }
            idx_list.push_back(idx);
            dist = 10e6;
        }
    }
    

    //Iterate over contours approximated by Contour class, number of contours depends on how many centers could be detected by clustering algorithm,
    //NOTE: When comparing against ground truth, number of detected centers has to match groundTruths vector
    for (size_t i = 0; i< groundTruths.size(); i++)
    {
    
        std::cout << "\nCOMPARE #"<< i+1<<"st  PAIR\n" << std::endl;

        SplineInterpolant_ptr_pair  f_Groundtruth = static_cast<SplineInterpolant_ptr_pair>(f_param(file, path, groundTruths[idx_list[i-1+groundTruths.size()]].first, groundTruths[idx_list[i-1+groundTruths.size()]].second,1000));
        
        ContourPairAnalyzer analyzer(f_Groundtruth,spline_vec[i], 1.0,results_dir);
        
        analyzer.analyzeContours(i+1);

    }
    file.close();
    return 0;
}