#include "discretization.h"

Discretization::Discretization(const double &sampleTime) : sampleTime_(sampleTime) {}

void Discretization::FowardEulerDiscretization(const StateMatrix &continuousA, const InputMatrix &continuousB)
{
    discreteA_ = StateMatrix::Identity() + continuousA * sampleTime_;
    discreteB_ = continuousB * sampleTime_;
}

// void Discretization::MatrixExponentialDiscretization(const StateMatrix &continuousA, const InputMatrix &continuousB)
// {
//     expMatrix_.setZero();
//     expMatrix_.block(0, 0, NUM_STATES, NUM_STATES) = continuousA;
//     expMatrix_.block(0, NUM_STATES, NUM_STATES, NUM_INPUTS) = continuousB;
//     expMatrix_ *= sampleTime_;
//     expDiscreteMatrix_ = expMatrix_.exp();
//     discreteA_ = expDiscreteMatrix_.block(0, 0, NUM_STATES, NUM_STATES);
//     discreteB_ = expDiscreteMatrix_.block(0, NUM_STATES, NUM_STATES, NUM_INPUTS);
// }