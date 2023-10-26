class F_x{
    spline1dinterpolant *s_;
    public:
    F_x(){}
    F_x(spline1dinterpolant *s) : s_(s){}
    double operator()(double x) const
    {
        return spline1dcalc(*s_, x);
    }
};
class D_f_x{
    spline1dinterpolant *s_;
    public:
    D_f_x(){}
    D_f_x(spline1dinterpolant *s) : s_(s){}
    double operator()(double x) const
    {
        double s, ds, dds;
        spline1ddiff(*s_, x, s, ds, dds);
        return ds;
    }
    
};

class Spline_function{
    public:
    Spline_function(){}
    Spline_function(spline1dinterpolant *s) :s_(s){}

    F_x f_s_x_()
    {
        return F_x(s_);
    }
    
    D_f_x d_f_s_x_()
    {   
        return D_f_x(s_);
    }
    spline1dinterpolant* s_;
};
class Contour_Integrand{
    
    F_x f_x_1_,f_x_2_;
    D_f_x d_f_x_1_, d_f_x_2_;
    public:
        Contour_Integrand(Spline_function *spline_1, Spline_function *spline_2 )
        {   
    
            f_x_1_ = spline_1->f_s_x_();
            f_x_2_ = spline_2->f_s_x_();
            d_f_x_1_ = spline_1->d_f_s_x_();
            d_f_x_2_ = spline_2->d_f_s_x_();


        }
        double operator()(double x)
        {
            return -0.5*(f_x_1_(x) * d_f_x_2_(x) - f_x_2_(x) * d_f_x_1_(x));
        }
};








// Newton-Raphson method for a 2D system
std::pair<double, double> newtonRaphson2D(double t0, double s0, double t_min, double t_max, double s_min, double s_max, int maxIter = 100, double tol = 1e-6) {
    double t = t0;
    double s = s0;

    for (int i = 0; i < maxIter; ++i) {
        ublas::vector<double> F(2);
        F(0) = f1(t, s);
        F(1) = f2(t, s);

        if (norm_2(F) < tol) {
            break;
        }

        ublas::matrix<double> J = jacobian(t, s);
        ublas::permutation_matrix<std::size_t> pm(J.size1());
        lu_factorize(J, pm);
        lu_substitute(J, pm, F);

        // Update and clip to bounds
        t = std::min(t_max, std::max(t_min, t - F(0)));
        s = std::min(s_max, std::max(s_min, s - F(1)));
    }

    return {t, s};
}

/*
    int n_samples_ = 10;
    double t[n_samples_];
    std::cout<<"sample size: "<<sizeof(t)/sizeof(double)<<std::endl;
    double f_1[n_samples_];
    double f_2[n_samples_];
    for (size_t i = 0; i < n_samples_; i++)
    {
        f_1[i] =0.5+sin(2*M_PI*i/(n_samples_-1));
        f_2[i] =cos(2*M_PI*i/(n_samples_-1));;
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
    
   
    
     using boost::math::quadrature::trapezoidal;
     // This function has an analytic form for its integral over a period: 2pi/3.
     //auto f_x = [](double x) { return sin(x)*sin(x)+cos(x)*cos(x); };
     
    Spline_function spline_1(&s1);
    Spline_function spline_2(&s2);

    Contour_Integrand c_i(&spline_1,&spline_2);

    int N = 1000;
    
    // Size of each sub-interval
    double a = 0.0;
    double b = 1.0;
    double step = (b - a) / N;
    
    std::vector<double> total_integral_vec;
    std::vector<double> t_val;
    double total_integral = 0;
    for (int i = 0; i < N; ++i) {
        double start = a + i * step;
        double end = start + step;
        
        // Integrate over the sub-interval [start, end]
        double ground_truth = boost::math::quadrature::trapezoidal(f, start, end);
        double approximation= trapezoidal(c_i, start, end);
        total_integral += ground_truth-approximation;
        total_integral_vec.push_back(total_integral);
        t_val.push_back((start+end)/2);
    }
    
    writeVectorsToCSV(t_val,total_integral_vec,"integral_vals.csv");
    return 0;
    */
   template <typename T1, typename T2>
void writeVectorsToCSV(const std::vector<T1>& vec1, const std::vector<T2>& vec2, const std::string& filename) {
    if (vec1.size() != vec2.size()) {
        std::cerr << "Vectors have different sizes!" << std::endl;
        return;
    }

    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Failed to open the file: " << filename << std::endl;
        return;
    }

    for (size_t i = 0; i < vec1.size(); ++i) {
        outfile << vec1[i] << "," << vec2[i] << "\n";
    }

    outfile.close();
}


// Define your functions
double f1(double t, double s) {
    return sin(t) + 0.5 -sin(s);
}
double d_ds_f1(double t, double s) {
    return -cos(s);
}
double d_dt_f1(double t, double s) {
    return -cos(t);
}

double f2(double t, double s) 
{
    return cos(t)-cos (s);
}
double d_ds_f2(double t, double s) 
{
    return sin(s);
}
double d_dt_f2(double t, double s) 
{
    return -sin(t);
}

// Define the Jacobian matrix
ublas::matrix<double> jacobian(double t, double s) {
    ublas::matrix<double> j(2, 2);
    // Fill in the partial derivatives
    // For this example, we're just putting in dummy values
    j(0, 0) = d_dt_f1(t,s);  // df1/dt
    j(0, 1) = d_ds_f1(t,s);  // df1/ds
    j(1, 0) = d_dt_f2(t,s);  // df2/dt
    j(1, 1) = d_ds_f2(t,s);  // df2/ds
    return j;
}
double f(double x)
{
    return sin(x)*sin(x)+cos(x)*cos(x);
}