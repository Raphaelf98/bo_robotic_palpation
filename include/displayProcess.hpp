#include<thread>
#include<mutex>
#include <condition_variable>
#include <chrono>
#include "matplotpp.h"
#include "contour.hpp"

class DisplayProcess : public MatPlot
{
private:
    std::mutex mtx_;
    
    std::condition_variable cv_;
    bool received_;
    bool sent_;
public:
    DisplayProcess(/* args */);
    ~DisplayProcess();
    void displayProcess()
};

