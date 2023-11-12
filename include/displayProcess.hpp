#include<thread>
#include<mutex>
#include <condition_variable>
#include <chrono>
#include "matplotpp.h"
#include "contour.hpp"

class DisplayProcess : public MatPlot
{
private:
    
public:
    DisplayProcess(std::mutex &mtx);
    ~DisplayProcess();
    void displayProcess()
};

