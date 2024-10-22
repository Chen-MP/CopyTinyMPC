#ifndef PENDULUM_DATA_H
#define PENDULUM_DATA_H

#include "../common/types_constants.h"

#ifndef USE_PENDULUM_CONFIG
#define USE_PENDULUM_CONFIG
#endif

const double M = 0.5;
const double m = 0.5;
const double l = 0.3;
const double g = 9.81;

double A_data[NUM_STATES * NUM_STATES] = 
    {0,                         0, 1, 0,
     0,                         0, 0, 1,
     0,                 m * g / M, 0, 0,
     0, (m * g + M * g) / (M * l), 0, 0};
     
double B_data[NUM_STATES * NUM_INPUTS] = 
    {          0, 
               0, 
           1 / M, 
     1 / (M * l)};

#endif