# Active Exploration Strategy for Tumor Localization in Robot-Assisted Surgery Using Bayesian Optimization
Robot-assisted minimally invasive surgery (RMIS) has become increasingly popular
in the resection of cancers. However, the lack of tactile feedback in clinical RMIS limits the
surgeon’s haptic understanding of tissue mechanics, making it hard to detect tissue abnormalities
(e.g., tumors) efficiently. This project will use Gaussian Processes (GP) to model the stiffness
distribution and develop an active exploration strategy to autonomously guide the robot to localize
the tumor profile via robot palpation using Bayesian Optimization (BO). Compared to exhaustively
palpating the entire organ to achieve tumor localization, the proposed framework aims at improving
efficiency with minimum times of exploration. Simulation and comparison experiments will be
conducted to verify the performance of the proposed strategy.

## Installation
### 1. Install dependencies 
For Ubuntu/Debian, the minimum dependencies (C/C++) can be optained by running: 
```
sudo apt install libboost-dev cmake g++
sudo apt install freeglut3-dev

```
The required BayesOpt C++ library to perform Bayesian optimizaiton needs to be cloned to a convenient location and build. For this run:
```
git clone https://github.com/rmcantin/bayesopt
cd bayesopt
cmake . 
make
sudo make install
```

### 2. Clone repository
Clone repository and build the package.
```
git clone https://git.tu-berlin.de/raphael/bo_robotic_palpation
cd bo_robotic_palpation
mkdir build && cd build
cmake ..
make
```
## Using the package 
### 1. Configure Bayesian Optimization parameters
Adjust optimization parameters by entering desired values under /config/bo_parameters.txt . 
### 2. Run the optimization
```
cd build
./display_gp ../config/bo_parameters.txt
```


