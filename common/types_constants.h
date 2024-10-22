#ifndef TYPES_CONSTANTS_H
#define TYPES_CONSTANTS_H

#include "../include/Eigen/Eigen/Dense"

#if defined(USE_PENDULUM_CONFIG) && defined(USE_TESTMAT_CONFIG)
#error "Cannot define both USE_PENDULUM_CONFIG and USE_TESTMAT_CONFIG"
#elif defined(USE_PENDULUM_CONFIG)
namespace Config {
    constexpr int NUM_STATES = 4;
    constexpr int NUM_INPUTS = 1;
    constexpr int CONTROL_HORIZON = 8;
}
#elif defined(USE_TESTMAT_CONFIG)
namespace Config {
    constexpr int NUM_STATES = 2;
    constexpr int NUM_INPUTS = 1;
    constexpr int CONTROL_HORIZON = 4;
}
#else
namespace Config {
    // NOTE: "No configuration selected, setting default values."
    constexpr int NUM_STATES = 5; 
    constexpr int NUM_INPUTS = 2;
    constexpr int CONTROL_HORIZON = 8;
}
#endif

constexpr int NUM_STATES = Config::NUM_STATES;
constexpr int NUM_INPUTS = Config::NUM_INPUTS;
constexpr int HORIZON = Config::CONTROL_HORIZON;

constexpr int HORIZON_STATES = NUM_STATES * HORIZON;
constexpr int HORIZON_INPUTS = NUM_INPUTS * HORIZON;

using StateVector = Eigen::Matrix<double, NUM_STATES, 1>;
using InputVector = Eigen::Matrix<double, NUM_INPUTS, 1>;

using StateMatrix = Eigen::Matrix<double, NUM_STATES, NUM_STATES>;
using InputMatrix = Eigen::Matrix<double, NUM_STATES, NUM_INPUTS>;

using StateWeightMatrix = Eigen::Matrix<double, NUM_STATES, NUM_STATES>;
using InputWeightMatrix = Eigen::Matrix<double, NUM_INPUTS, NUM_INPUTS>;
using StateWeightDiagonalMatrix = Eigen::DiagonalMatrix<double, NUM_STATES>;
using InputWeightDiagonalMatrix = Eigen::DiagonalMatrix<double, NUM_INPUTS>;

using ExtendStateWeightVector = Eigen::Matrix<double, HORIZON_STATES, 1>;
using ExtendInputWeightVector = Eigen::Matrix<double, HORIZON_INPUTS, 1>;
using ExtendStateWeightMatrix = Eigen::Matrix<double, HORIZON_STATES, HORIZON_STATES>;
using ExtendInputWeightMatrix = Eigen::Matrix<double, HORIZON_INPUTS, HORIZON_INPUTS>;

using ExtendStateBoundVector = ExtendStateWeightVector;
using ExtendInputBoundVector = ExtendInputWeightVector;
// using ExtendStateConstrainMatrix = ExtendStateWeightMatrix;
// using ExtendInputConstrainMatrix = ExtendInputWeightMatrix;

using LqrGainMatrix = Eigen::Matrix<double, NUM_INPUTS, NUM_STATES>;
using LqrRiccatiMatrix = Eigen::Matrix<double, NUM_STATES, NUM_STATES>;    

extern Eigen::IOFormat MatrixApiFmt;

#endif