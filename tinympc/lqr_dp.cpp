#include <iostream>
#include <chrono> 

#include "lqr_dp.h"

LqrDp::LqrDp(const StateMatrix &discreteA, const InputMatrix &discreteB, 
             const StateVector &QdiagVec, const InputVector &RdiagVec, 
             const double &rho)
        : discreteA_(discreteA), discreteB_(discreteB), QdiagVec_(QdiagVec), RdiagVec_(RdiagVec), 
          rho_(rho)
{
    // this->P_kplus = rho * tinyMatrix::Ones(x_size, 1).array().matrix().asDiagonal();
}

void LqrDp::SetWeightMatrix()
{
    this->Qmat_ = QdiagVec_.asDiagonal();
    this->Rmat_ = RdiagVec_.asDiagonal();
    this->P_kplus_ = Qmat_;
}

void LqrDp::SetWeightMatrixBarAdmm()
{
    this->QmatBar_ = Qmat_ + rho_ * StateWeightMatrix::Identity();
    this->RmatBar_ = Rmat_ + rho_ * InputWeightMatrix::Identity();
    
    if (weightMatBarVerbose_) {
        std::cout << "QmatBar = " << QmatBar_.format(MatrixApiFmt) << std::endl;
        std::cout << "RmatBar = " << RmatBar_.format(MatrixApiFmt) << std::endl;
        printf("rho = %f. \n", rho_);
    }
}

void LqrDp::PrecomputeGainAndCache()
{
    for (int i = 0; i < 1000; ++i) {
        this->tempSymMat_.noalias() = discreteB_.transpose() * P_kplus_ * discreteB_;
        tempSymMat_.diagonal() += RmatBar_.diagonal();
        tempSymMatInv_ = tempSymMat_.inverse();
        // tempSymMatInv_ = tempSymMat_.ldlt().solve(Eigen::MatrixXd::Identity(tempSymMat_.rows(), tempSymMat_.cols()));
        this->K_k_.noalias() = tempSymMatInv_ * discreteB_.transpose();
        K_k_ *= P_kplus_;
        K_k_ *= discreteA_;
        // this->K_k_.noalias() = (RmatBar_ + discreteB_.transpose() * P_kplus_ * discreteB_).inverse() * discreteB_.transpose() * P_kplus_ * discreteA_;

        this->AmBKMat_.noalias() = -(discreteB_ * K_k_);
        AmBKMat_.noalias() += discreteA_;
        this->P_k_.noalias() = discreteA_.transpose() * P_kplus_;
        P_k_ *= AmBKMat_;
        P_k_.diagonal() += QmatBar_.diagonal();

        if ((K_k_ - K_minus_).cwiseAbs().maxCoeff() < 1e-5) {
            if (convergenceVerbose_) {
                std::cout << "Kinf converged after " << i + 1 << " iterations" << std::endl;
            }
            break;
        }

        K_minus_ = K_k_;
        P_kplus_ = P_k_;
    }

    this->tempQuuMat_.noalias() = discreteB_.transpose() * P_k_ * discreteB_;
    tempQuuMat_.diagonal() += RmatBar_.diagonal();
    this->QuuInv_ = tempQuuMat_.inverse();
    // this->QuuInv_  = tempQuuMat_.ldlt().solve(Eigen::MatrixXd::Identity(tempQuuMat_.rows(), tempQuuMat_.cols()));
    
    this->AmBKMatTranspose_.noalias() = -(discreteB_ * K_k_);
    AmBKMatTranspose_.noalias() += discreteA_;
    AmBKMatTranspose_.transposeInPlace();
    // this->AmBKMatTranspose_ = (discreteA_ - discreteB_ * K_k_).transpose();

    if (gainVerbose_) {
        std::cout << "K_k_ = " << K_k_.format(MatrixApiFmt) << std::endl;
        std::cout << "P_k_ = " << P_k_.format(MatrixApiFmt) << std::endl;
        std::cout << "\n Precomputation finished! \n" << std::endl;
    } else if (cacheVerbose_) {
        std::cout << "QuuInv_ = " << QuuInv_.format(MatrixApiFmt) << std::endl;
        std::cout << "AmBKMatTranspose = " << AmBKMatTranspose_.format(MatrixApiFmt) << std::endl;
        std::cout << "\n Precomputation finished! \n" << std::endl;
    } else {
        std::cout << "\n Precomputation finished! \n" << std::endl;
    }
}

void LqrDp::ComputeLqrDp()
{
    SetWeightMatrix();
    SetWeightMatrixBarAdmm();
    // auto start = std::chrono::high_resolution_clock::now();
    PrecomputeGainAndCache();
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> duration_ms = end - start;
    // printf("Execution time: %f milliseconds. \n", duration_ms.count());
}


// LQR
StateVector LqrDp::GetNextState(const StateVector &xInit)
{   
    this->uLqrDpVec_ = K_k_ * (StateVector::Zero() -  xInit);
    this->xNextVec_ = discreteA_ * xInit + discreteB_ * uLqrDpVec_;
    return xNextVec_;
}

StateVector LqrDp::GetNextState(const StateVector &xRefCurr, const StateVector &xInit)
{   
    this->uLqrDpVec_ = K_k_ * (xRefCurr -  xInit);
    this->xNextVec_ = discreteA_ * xInit + discreteB_ * uLqrDpVec_;
    return xNextVec_;
}
