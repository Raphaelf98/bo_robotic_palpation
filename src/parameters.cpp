#include "parameters.hpp"




void TumorModelParameters::loadModel(bayesopt::utils::FileParser &fp, TumorModelParameters &cp)
{
    fp.readOrWrite("triangle_low", cp.triangle_low);
    fp.readOrWrite("triangle_high", cp.triangle_high);
    fp.readOrWrite("triangle_radius", cp.triangle_radius);
    fp.readOrWrite("triangle_x_trans", cp.triangle_x_trans);
    fp.readOrWrite("triangle_y_trans", cp.triangle_y_trans);
    fp.readOrWrite("triangle_epsilon", cp.triangle_epsilon);

    fp.readOrWrite("rectangle_low", cp.rectangle_low);
    fp.readOrWrite("rectangle_high", cp.rectangle_high);
    fp.readOrWrite("rectangle_radius", cp.rectangle_radius);
    fp.readOrWrite("rectangle_x_trans", cp.rectangle_x_trans);
    fp.readOrWrite("rectangle_y_trans", cp.rectangle_y_trans);
    fp.readOrWrite("rectangle_epsilon", cp.rectangle_epsilon);

    fp.readOrWrite("circle_low", cp.circle_low);
    fp.readOrWrite("circle_high", cp.circle_high);
    fp.readOrWrite("circle_radius", cp.circle_radius);
    fp.readOrWrite("circle_x_trans", cp.circle_x_trans);
    fp.readOrWrite("circle_y_trans", cp.circle_y_trans);
    fp.readOrWrite("circle_epsilon", cp.circle_epsilon);

    fp.readOrWrite("two_circles_low", cp.two_circles_low);
    fp.readOrWrite("two_circles_high", cp.two_circles_high);
    fp.readOrWrite("two_circles_radius_1", cp.two_circles_radius_1);
    fp.readOrWrite("two_circles_radius_2", cp.two_circles_radius_2);
    fp.readOrWrite("two_circles_x_trans_1", cp.two_circles_x_trans_1);
    fp.readOrWrite("two_circles_x_trans_2", cp.two_circles_x_trans_2);
    fp.readOrWrite("two_circles_y_trans_1", cp.two_circles_y_trans_1);
    fp.readOrWrite("two_circles_y_trans_2", cp.two_circles_y_trans_2);
    fp.readOrWrite("two_circles_epsilon", cp.two_circles_epsilon);

    

}
void ContourParamters::loadContour(bayesopt::utils::FileParser &fp, ContourParamters &cp)
{
    fp.readOrWrite("n_exploration_directions", cp.n_exploration_directions);
    fp.readOrWrite("c_points", cp.c_points);
    fp.readOrWrite("means_shift_bandwidth", cp.means_shift_bandwidth);
    fp.readOrWrite("lim_steps", cp.lim_steps);
    fp.readOrWrite("threshold_multiplier", cp.threshold_multiplier);
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
