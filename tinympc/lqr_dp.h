#ifndef LQR_DP_H
#define LQR_DP_H

#include "../common/types_constants.h"

class LqrDp {
public:
    LqrDp(const StateMatrix &discreteA, const InputMatrix &discreteB, 
          const StateVector &QdiagVec, const InputVector &RdiagVec, 
          const double &rho);
    virtual ~LqrDp() = default ;

    void SetWeightMatrix();
    void SetWeightMatrixBarAdmm();
    void PrecomputeGainAndCache(); 
    void ComputeLqrDp();

    StateVector GetNextState(const StateVector &xInit);
    StateVector GetNextState(const StateVector &xRefCurr, const StateVector &xInit);
    
    StateMatrix GetDiscreteStateMatrix() const { return discreteA_; }
    InputMatrix GetDiscreteInputMatrix() const { return discreteB_; }
    StateWeightMatrix GetStateWeightMatrix() const { return Qmat_; }
    InputWeightMatrix GetInputWeightMatrix() const { return Rmat_; }
    LqrGainMatrix GetLqrGainMatrixKk() const { return K_k_; }
    LqrRiccatiMatrix GetRiccatiMatrixPk() const { return P_k_; }
    Eigen::Matrix<double, NUM_INPUTS, NUM_INPUTS> GetCacheMatQuuInv() const { return QuuInv_; }
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> GetCacheMatAmBKMatTran() const { return AmBKMatTranspose_; }

private:
    StateMatrix discreteA_;
    InputMatrix discreteB_;
    StateVector QdiagVec_;
    InputVector RdiagVec_;
    double rho_;
    bool verbose_;
    bool equ_verbose_ = false;

    StateWeightMatrix Qmat_;
    InputWeightMatrix Rmat_; 
    StateWeightMatrix QmatBar_;
    InputWeightMatrix RmatBar_; 
    bool weightMatBarVerbose_ =false;

    Eigen::Matrix<double, NUM_INPUTS, NUM_INPUTS> tempSymMat_;
    Eigen::Matrix<double, NUM_INPUTS, NUM_INPUTS> tempSymMatInv_;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> AmBKMat_;
    LqrGainMatrix K_minus_, K_k_;
    LqrRiccatiMatrix P_kplus_, P_k_ ;
    bool convergenceVerbose_ = false;

    Eigen::Matrix<double, NUM_INPUTS, NUM_INPUTS> tempQuuMat_;
    Eigen::Matrix<double, NUM_INPUTS, NUM_INPUTS> QuuInv_;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> AmBKMatTranspose_;
    bool gainVerbose_ = false;
    bool cacheVerbose_ = false;

    StateVector xNextVec_;
    InputVector uLqrDpVec_;
};

#endif