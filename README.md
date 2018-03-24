# Unscented Kalman Filter

In this project utilize an Unscented Kalman Filter to estimate the state of a moving object of interest with noisy lidar and radar measurements. 

This project involves the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

This repository includes two files that can be used to set up and install [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) for either Linux or Mac systems. For windows you can use either Docker, VMware, or even [Windows 10 Bash on Ubuntu](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/) to install uWebSocketIO. Please see [this concept in the classroom](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/16cf4a78-4fc7-49e1-8621-3450ca938b77) for the required version and installation scripts.

Once the install for uWebSocketIO is complete, the main program can be built and ran by doing the following from the project top directory.

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./UnscentedKF

## Steps

There are 2 main steps to take

### Predition

* Generate Sigma Points
* Predict Sigma Points
* Predict Mean and Covariance

### Update

* Predict Measurement
* Update State

The purpose of this project is to measure RSME values of the calculated value of position x, position y, velocity x, and velocity y compare with ground truth measurement sent from the simulated program and as result, you can see we can get RMSE <= [.09, .10, 0.40, 0.30] which is lower if compared [Extended Kalman Filters](https://github.com/anonymint/extended-kalman-filter) I worked on earlier ([0.11, 0.11, 0.52, 0.52]) or simply say performs better in term of accuracy.

#### Result from dataset 1

##### Radar only

![Radar only](asset/radar_only_1.png)

##### Lidar only

![Lidar only](asset/lidar_only_1.png)

##### Combine both Radar and Lidar

![combine both Radar and Lidar](asset/lidar_radar_1.png)

#### Result from dataset 2

##### Radar only

![Radar only](asset/radar_only_2.png)

##### Lidar only

![Lidar only](asset/lidar_only_2.png)

##### Combine both Radar and Lidar

![combine both Radar and Lidar](asset/lidar_radar_2.png)

#### Results

As we can see from result captured above, either Radar and Lidar can beat combined sensors as result of RSME value. In addition to that, Unscented Kalman Filter can give a better result than [Extended Kalman Filters](https://github.com/anonymint/extended-kalman-filter).