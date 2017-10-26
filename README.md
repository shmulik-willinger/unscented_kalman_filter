# Unscented Kalman Filter

[![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)

Overview
---

In this project the goal is to implement an `Unscented Kalman Filter` with C++. The kalman filter estimate the state of a moving object of interest with noisy LIDAR and RADAR measurements, while obtaining RMSE values that are very low. The communication between the project and the simulator is done using WebSocket.

Compared to the Extended kalman filter where we used a constant velocity model, In the Unscented kalman filter we are using a `constant turn rate and velocity magnitude model` (CTRV), which is a great model to accurately track non-linear motion.

The CTRV model:
![]( https://github.com/shmulik-willinger/unscented_kalman_filter/blob/master/readme_img/ctrv_model.jpg?raw=true)

The three most important probabilities of the UKF as a central fusion tool for self driving cars:
* The UKF is able to take noisy measurement data as input and provide a smooth position and velocity estimation of dynamic objects around it, without introducing a delay.
* We provide an estimation of the orientation and the yaw rate of other vehicles using sensors that can't even directly observe these things.
* The UKF give information on how precise the result is, because it always provides a covariance matrix for every estimation, and we know that this covariance matrix is realistic since the UKF performs a consistency check.


Prerequisites and Executing
---

This project requires the following dependencies:

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)


  This project involves the Udacity Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

  The WebSocket and other initial info regarding dependencies  installation can be found  [here](https://github.com/udacity/CarND-Extended-Kalman-Filter-Project)

Build Instructions:

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
   * On windows, you may need to run: `cmake .. -G "Unix Makefiles" && make`

Running the project:
1. Run the project main file: `./UnscentedKF `
2. You will get the following output:
 `Listening to port 4567. Connected!!! `
 3. Run the Simulator on `project 1` and start the process on one of the datasets.

Process results
---

The simulator final state after running the EKL:

Dataset 1              |  Dataset 2
:---------------------:|:---------------------:
![]( https://github.com/shmulik-willinger/unscented_kalman_filter/blob/master/readme_img/dataset_1.jpg?raw=true)  |  ![]( https://github.com/shmulik-willinger/unscented_kalman_filter/blob/master/readme_img/dataset_2.jpg?raw=true)

The Unscented Kalman Filter Residual error is calculated by mean squared error, and the accuracy is very good (low) for both datasets:

| Input |   MSE for dataset 1   | MSE for dataset 2   |
| ----- | ------- | ------- |
|  px   | 0.0737 | 0.0738 |
|  py   | 0.0846 | 0.0752 |
|  vx   | 0.3236 | 0.3292 |
|  vy   | 0.2402 | 0.2922 |

The results are considerably better than the [Extended Kalman Filter project.](https://github.com/shmulik-willinger/extended_kalman_filter)

The Kalman Filter algorithm handles the first measurements appropriately. The measurement is handled at src/ukf.cpp from line 102 to line 138

The estimation part is being done by two steps: the algorithm first predicts then updates. Predict process can be found at src/ukf.cpp from line 450. Update process can be found at lines 153-182.
![]( https://github.com/shmulik-willinger/extended_kalman_filter/blob/master/readme_img/loop.jpg?raw=true)

The project can handle both RADAR and LIDAR measurements. For experience, I ran the process each time for only one of the inputs (lidar / radar), observing how the estimations change when running against a single sensor type. Here are the result:

only LIDAR               |  only RADAR
:---------------------:|:---------------------:
![]( https://github.com/shmulik-willinger/unscented_kalman_filter/blob/master/readme_img/radar_only.jpg?raw=true)  |  ![]( https://github.com/shmulik-willinger/unscented_kalman_filter/blob/master/readme_img/lidar_only.jpg?raw=true)
![]( https://github.com/shmulik-willinger/unscented_kalman_filter/blob/master/readme_img/radar_only_2.jpg?raw=true)  |  ![]( https://github.com/shmulik-willinger/unscented_kalman_filter/blob/master/readme_img/lidar_only_2.jpg?raw=true)

Its easy to see that both results have higher RMSE comparing to the process with both inputs together. The RADAR screw up a large part of the way, giving the worst result.

The algorithm avoid unnecessary calculations. For example: The calculation of the weights matrix is being done only at the constructor, skipping the 'prediction' step if more than one measurment came in the same time, and more.
The code is clean and readable, including verification for dividing by Zero and wrong inputs.

The video below shows how the simulator performs when running with this Kalman filter project to track the object.
* LIDAR measurements are red circles
* RADAR measurements are blue circles with an arrow pointing in the direction of the observed angle
* Estimation markers are green triangles

[![video output](https://github.com/shmulik-willinger/unscented_kalman_filter/blob/master/readme_img/dataset_1.gif)]

I'm curently trying to finish the extra challenge where We need to catch a run away car with the implementation of an Unscented Kalman Filter. The capture car and the run away car have the same max velocity, so in order for the capture vehicle to catch the run away car, it will need to predict where the car will be ahead of time. curently this is how the process looks when running in the simulator:

[![video output](https://github.com/shmulik-willinger/unscented_kalman_filter/blob/master/readme_img/Catch-Run-Away-Car.gif)]

comment:
 In order for the project to run faster, I canceled the printing of the data during the running process.
