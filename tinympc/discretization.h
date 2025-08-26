#ifndef DISCRETIZATION_H
#define DISCRETIZATION_H

// #include <eigen3/unsupported/Eigen/MatrixFunctions>

#include "../common/types_constants.h"

class Discretization {
public:
    Discretization(const double &sampleTime);
    virtual ~Discretization() = default;

    void FowardEulerDiscretization(const StateMatrix &continuousA, const InputMatrix &continuousB);
    // void MatrixExponentialDiscretization(const StateMatrix &continuousA, const InputMatrix &continuousB);
    
    StateMatrix GetDiscreteA() const { return discreteA_; }
    InputMatrix GetDiscreteB() const { return discreteB_; }
    
private:
    StateMatrix discreteA_;
    InputMatrix discreteB_;
    double sampleTime_;
    
    Eigen::Matrix<double, NUM_STATES + NUM_INPUTS, NUM_STATES + NUM_INPUTS> expMatrix_;
    Eigen::Matrix<double, NUM_STATES + NUM_INPUTS, NUM_STATES + NUM_INPUTS> expDiscreteMatrix_;
};

#endif
