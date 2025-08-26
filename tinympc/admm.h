#ifndef ADMM_H
#define ADMM_H

#include "lqr_dp.h"

class Admm : public LqrDp {
public:
    Admm(const StateMatrix &discreteA, const InputMatrix &discreteB, 
         const StateVector &QdiagVec, const InputVector &RdiagVec, 
         const double &rho);
    virtual ~Admm() = default ;

    void GetLqrDpParams();

    void SetXrefVecMat();
    void SetXrefVecMat(const StateVector &xRef);

    void SetUrefVecMat();
    void SetUrefVecMat(const InputVector &uRef);

    void SetXmatTinyMpcInit(const StateVector &xInit) { xTinyMpcMat_.col(0) = xInit; }

    void DisableStateConstrain() { this->enStateBound_ = false; }
    void EnableStateConstrain(const StateVector &xMax, const StateVector &xMin);

    void DisableInputConstrain() { this->enInputBound_ = false; }
    void EnableInputConstrain(const InputVector &uMax, const InputVector &uMin);

    void UpdateLinearCost(); 
    void BackwardPassGrad();
    void ForwardPass(); 
    void UpdateSlack(); 
    void UpdateDual(); 
    bool TerminationCondition();

    int SolveTinyMpc() ; 

    InputVector GetUvec() const { return this->uTinyMpcMat_.col(0); }
    StateVector GetXvecNext(const StateVector &xInit);

private:
    double rho_;
    StateMatrix discreteA_;
    InputMatrix discreteB_;
    StateWeightMatrix Qmat_;
    InputWeightMatrix Rmat_;

    // StateVector xRefVec_;
    // InputVector uRefVec_; 
    Eigen::Matrix<double, NUM_STATES, HORIZON> xRefMat_;
    Eigen::Matrix<double, NUM_INPUTS, HORIZON - 1> uRefMat_;

    Eigen::Matrix<double, NUM_STATES, HORIZON> xTinyMpcMat_;
    Eigen::Matrix<double, NUM_INPUTS, HORIZON - 1> uTinyMpcMat_;

    Eigen::Matrix<double, NUM_STATES, HORIZON> xMinMat_, xMaxMat_;

    Eigen::Matrix<double, HORIZON_STATES, 1> xMinHorizonVec_, xMaxHorizonVec_;
    Eigen::Matrix<double, HORIZON_INPUTS, 1> uMinHorizonVec_, uMaxHorizonVec_;

    Eigen::Matrix<double, NUM_INPUTS, HORIZON - 1> uMinMat_, uMaxMat_;

    LqrGainMatrix K_k_;
    LqrRiccatiMatrix P_k_;
    
    Eigen::Matrix<double, NUM_INPUTS, NUM_INPUTS> QuuInv_;
    Eigen::Matrix<double, NUM_STATES, NUM_STATES> AmBKMatTranspose_;

    bool enStateBound_ = false; 
    bool enInputBound_ = false;

    int iter_;
    int maxIter_ = 100;
    int checkTermination_ = 1;
    double absPriTol_ = 1e-3;
    double absDuaTol_ = 1e-3;

    double primalResidualState_ ; 
    double dualResidualState_ ; 
    double primalResidualInput_ ; 
    double dualResidualInput_ ;

    bool solved_;

    Eigen::Matrix<double, NUM_STATES, HORIZON> qLinearStateCost_;    // x_size x N
    Eigen::Matrix<double, NUM_INPUTS, HORIZON - 1> rLinearInputCost_;    // u_size x N-1

    Eigen::Matrix<double, NUM_STATES, HORIZON> p;    // x_size x N
    Eigen::Matrix<double, NUM_INPUTS, HORIZON - 1> d;    // u_size x N-1

    Eigen::Matrix<double, NUM_STATES, HORIZON> vSlackStateVar_;    // x_size x N
    Eigen::Matrix<double, NUM_STATES, HORIZON> vNewSlack_; // x_size x N
    Eigen::Matrix<double, NUM_INPUTS, HORIZON - 1> zSlackInputVar_;    // u_size x N-1
    Eigen::Matrix<double, NUM_INPUTS, HORIZON - 1> zNewSlack_; // u_size x N-1
 
    Eigen::Matrix<double, NUM_STATES, HORIZON> gDualStateVar_;    // x_size x N
    Eigen::Matrix<double, NUM_INPUTS, HORIZON - 1> yDualInputVar_;    // u_size x N-1

    StateVector xTinyMpcNextVec_;
};

#endif






