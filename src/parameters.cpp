#include "parameters.hpp"




void TumorModelParameters::loadModel(bayesopt::utils::FileParser &fp, TumorModelParameters &cp)
{
    fp.readOrWrite("triangle_low", cp.triangle_low);
    fp.readOrWrite("triangle_high", cp.triangle_high);
    fp.readOrWrite("triangle_radius", cp.triangle_radius);
    fp.readOrWrite("triangle_x_trans", cp.triangle_x_trans);
    fp.readOrWrite("triangle_y_trans", cp.triangle_y_trans);
    fp.readOrWrite("triangle_epsilon", cp.triangle_epsilon);
    fp.readOrWrite("triangle_noise", cp.triangle_noise);

    fp.readOrWrite("rectangle_low", cp.rectangle_low);
    fp.readOrWrite("rectangle_high", cp.rectangle_high);
    fp.readOrWrite("rectangle_radius", cp.rectangle_radius);
    fp.readOrWrite("rectangle_x_trans", cp.rectangle_x_trans);
    fp.readOrWrite("rectangle_y_trans", cp.rectangle_y_trans);
    fp.readOrWrite("rectangle_epsilon", cp.rectangle_epsilon);
    fp.readOrWrite("rectangle_noise", cp.rectangle_noise);

    fp.readOrWrite("circle_low", cp.circle_low);
    fp.readOrWrite("circle_high", cp.circle_high);
    fp.readOrWrite("circle_radius", cp.circle_radius);
    fp.readOrWrite("circle_x_trans", cp.circle_x_trans);
    fp.readOrWrite("circle_y_trans", cp.circle_y_trans);
    fp.readOrWrite("circle_epsilon", cp.circle_epsilon);
    fp.readOrWrite("circle_noise", cp.circle_noise);

    fp.readOrWrite("two_circles_low", cp.two_circles_low);
    fp.readOrWrite("two_circles_high", cp.two_circles_high);
    fp.readOrWrite("two_circles_radius_1", cp.two_circles_radius_1);
    fp.readOrWrite("two_circles_radius_2", cp.two_circles_radius_2);
    fp.readOrWrite("two_circles_x_trans_1", cp.two_circles_x_trans_1);
    fp.readOrWrite("two_circles_x_trans_2", cp.two_circles_x_trans_2);
    fp.readOrWrite("two_circles_y_trans_1", cp.two_circles_y_trans_1);
    fp.readOrWrite("two_circles_y_trans_2", cp.two_circles_y_trans_2);
    fp.readOrWrite("two_circles_epsilon", cp.two_circles_epsilon);
    fp.readOrWrite("two_circles_noise", cp.two_circles_noise);

    

}
 void TumorModelParameters::printParameters() {
        std::cout << "Triangle Parameters:" << std::endl;
        std::cout << "Low: " << triangle_low << ", High: " << triangle_high << ", Radius: " << triangle_radius 
                  << ", X Trans: " << triangle_x_trans << ", Y Trans: " << triangle_y_trans 
                  << ", Epsilon: " << triangle_epsilon << ", Noise: " << triangle_noise << std::endl;

        std::cout << "Rectangle Parameters:" << std::endl;
        std::cout << "Low: " << rectangle_low << ", High: " << rectangle_high << ", Radius: " << rectangle_radius 
                  << ", X Trans: " << rectangle_x_trans << ", Y Trans: " << rectangle_y_trans 
                  << ", Epsilon: " << rectangle_epsilon << ", Noise:" << rectangle_noise << std::endl;

    std::cout << "Circle Parameters:" << std::endl;
    std::cout << "Low: " << circle_low << ", High: " << circle_high << ", Radius: " << circle_radius 
              << ", X Trans: " << circle_x_trans << ", Y Trans: " << circle_y_trans 
              << ", Epsilon: " << circle_epsilon << ", Noise: " << circle_noise << std::endl;

    std::cout << "Two Circles Parameters:" << std::endl;
    std::cout << "Low: " << two_circles_low << ", High: " << two_circles_high 
              << ", Radius 1: " << two_circles_radius_1 << ", Radius 2: " << two_circles_radius_2
              << ", X Trans 1: " << two_circles_x_trans_1 << ", X Trans 2: " << two_circles_x_trans_2 
              << ", Y Trans 1: " << two_circles_y_trans_1 << ", Y Trans 2: " << two_circles_y_trans_2
              << ", Epsilon: " << two_circles_epsilon << ", Noise: " << two_circles_noise << std::endl;
}
void ContourParamters::loadContour(bayesopt::utils::FileParser &fp, ContourParamters &cp)
{
    fp.readOrWrite("n_exploration_directions", cp.n_exploration_directions);
    fp.readOrWrite("c_points", cp.c_points);
    fp.readOrWrite("means_shift_bandwidth", cp.means_shift_bandwidth);
    fp.readOrWrite("lim_steps", cp.lim_steps);
    fp.readOrWrite("threshold_multiplier", cp.threshold_multiplier);
    fp.readOrWrite("tumor_stiffness_guess_low", cp.tumor_stiffness_guess_low);
    fp.readOrWrite("tumor_stiffness_guess_high", cp.tumor_stiffness_guess_high);
}

bool TumorModelParameters::loadModelParameters(std::string filename, TumorModelParameters &cp)
{
    bayesopt::utils::FileParser fp(filename);
    if(!fp.fileExists())
    {
        return false;
    }
    
    fp.openInput();
    loadModel(fp,cp);
    return true;

}

bool ContourParamters::loadContourParameters(std::string filename, ContourParamters &cp)
{
    bayesopt::utils::FileParser fp(filename);
    if(!fp.fileExists())
    {
        return false;
    }
    
    fp.openInput();
    loadContour(fp,cp);
    return true;
}
  void ContourParamters::PrintParameters() {
        std::cout << "Contour Parameters:" << std::endl;
        std::cout << "Number of exploration directions: " << n_exploration_directions << std::endl;
        std::cout << "Number of contour points (c_points): " << c_points << std::endl;
        std::cout << "Mean shift bandwidth: " << means_shift_bandwidth << std::endl;
        std::cout << "Limit steps (lim_steps): " << lim_steps << std::endl;
        std::cout << "Threshold multiplier: " << threshold_multiplier << std::endl;
        std::cout << "Tumor stiffness guess low: " << tumor_stiffness_guess_low << std::endl;
        std::cout << "Tumor stiffness guess high: " << tumor_stiffness_guess_high << std::endl;

    }