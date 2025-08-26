#include <iostream>

#include "admm.h"

Admm::Admm(const StateMatrix &discreteA, const InputMatrix &discreteB, 
           const StateVector &QdiagVec, const InputVector &RdiagVec, 
           const double &rho)
        : LqrDp(discreteA, discreteB, QdiagVec, RdiagVec, rho), rho_(rho)
{
    // settings->check_termination = 1 ; // 能整除任何数
    // settings->abs_pri_tol = 1e-3;
    // settings->abs_dua_tol = 1e-3;
    // settings->max_iter = 100;

    // this->x = tinyMatrix::Zero(x_size, N);
    // this->u = tinyMatrix::Zero(u_size, N-1);

    this->qLinearStateCost_.setZero();
    this->rLinearInputCost_.setZero();

    this->p.setZero();
    this->d.setZero();

    this->vSlackStateVar_.setZero();
    this->vNewSlack_.setZero();
    this->zSlackInputVar_.setZero();
    this->zNewSlack_.setZero();
    
    this->gDualStateVar_.setZero();
    this->yDualInputVar_.setZero();
    // solver->settings->check_termination = 1; 

}

void Admm::GetLqrDpParams()
{
    ComputeLqrDp();
    this->discreteA_ = GetDiscreteStateMatrix();
    this->discreteB_ = GetDiscreteInputMatrix();
    this->Qmat_ = GetStateWeightMatrix();
    this->Rmat_ = GetInputWeightMatrix();
    this->K_k_ = GetLqrGainMatrixKk();
    this->P_k_ = GetRiccatiMatrixPk();
    this->AmBKMatTranspose_ = GetCacheMatAmBKMatTran();
}

void Admm::SetXrefVecMat()
{
    // this->xRefVec_.setZero();
    this->xRefMat_.setZero();
}
void Admm::SetXrefVecMat(const StateVector &xRef)
{
    // this->xRefVec_ = xRef;
    this->xRefMat_ = xRef.replicate(1, HORIZON);
}
void Admm::SetUrefVecMat()
{
    // this->uRefVec_.setZero();
    this->uRefMat_.setZero();
}
void Admm::SetUrefVecMat(const InputVector &uRefVec)
{
    // this->uRefVec_ = uRefVec;
    this->uRefMat_ = uRefVec.replicate(1, HORIZON - 1);
}

// TODO
void Admm::EnableStateConstrain(const StateVector &xMax, const StateVector &xMin)
{
    /* tinympc */
    this->enStateBound_ = true;
    this->xMaxMat_ = xMax.replicate(1, HORIZON);
    this->xMinMat_ = xMin.replicate(1, HORIZON);
    // std::cout << "xMinMat_ = " << xMinMat_.format(MatrixApiFmt) << std::endl;
}
void Admm::EnableInputConstrain(const InputVector &uMax, const InputVector &uMin)
{
    this->enInputBound_ = true;
    this->uMaxMat_ = uMax.replicate(1, HORIZON - 1);
    this->uMinMat_ = uMin.replicate(1, HORIZON - 1);
}


/* tinympc */
/* 更新 prime variables */
void Admm::ForwardPass()
{
    for (int i = 0; i < HORIZON - 1; ++i) {
        (this->uTinyMpcMat_.col(i)).noalias() = -K_k_.lazyProduct(xTinyMpcMat_.col(i)) - d.col(i);
        (this->xTinyMpcMat_.col(i + 1)).noalias() = discreteA_.lazyProduct(xTinyMpcMat_.col(i)) + discreteB_.lazyProduct(uTinyMpcMat_.col(i));
    }
}

void Admm::UpdateSlack()
{
    this->zNewSlack_ = uTinyMpcMat_ + yDualInputVar_;
    this->vNewSlack_ = xTinyMpcMat_ + gDualStateVar_;

    if (this->enStateBound_) {
        vNewSlack_ = this->xMaxMat_.cwiseMin(this->xMinMat_.cwiseMax(vNewSlack_));
        // std::cout << "vNewSlack_ = " << vNewSlack_(0,0) << std::endl;
    }

    if (this->enInputBound_) {
        zNewSlack_ = this->uMaxMat_.cwiseMin(this->uMinMat_.cwiseMax(zNewSlack_));
        // std::cout << "zNewSlack_ = " << zNewSlack_(0,0) << std::endl;
    }
}

void Admm::UpdateDual()
{
    this->yDualInputVar_ = yDualInputVar_ + uTinyMpcMat_ - zNewSlack_;
    this->gDualStateVar_ = gDualStateVar_ + xTinyMpcMat_ - vNewSlack_;
}

void Admm::UpdateLinearCost()
{
    // std::cout << "xRefMat_ = " << xRefMat_.format(MatrixApiFmt) << std::endl;
    this->rLinearInputCost_ = -(uRefMat_.array().colwise() * Rmat_.diagonal().array()); 
    rLinearInputCost_.noalias() -= rho_ * (zNewSlack_ - yDualInputVar_);
    this->qLinearStateCost_ = -(xRefMat_.array().colwise() * Qmat_.diagonal().array());
    qLinearStateCost_.noalias() -= rho_ * (vNewSlack_ - gDualStateVar_);
    // 更新最后一项
    p.col(HORIZON - 1) = (-(xRefMat_.col(HORIZON - 1).transpose().lazyProduct(P_k_)));
    (p.col(HORIZON - 1)).noalias() -= rho_ * (vNewSlack_.col(HORIZON - 1) - gDualStateVar_.col(HORIZON - 1));
}

bool Admm::TerminationCondition()
{
    if (this->iter_ % this->checkTermination_ == 0) {
        // cwiseAbs(): 用于对向量或矩阵的每个元素取绝对值。
        // .maxCoeff(): 用于计算向量或矩阵中的最大元素值。
        this->primalResidualState_ = (this->xTinyMpcMat_ - this->vNewSlack_).cwiseAbs().maxCoeff();
        this->dualResidualState_ = ((this->vSlackStateVar_ - this->vNewSlack_).cwiseAbs().maxCoeff()) * this->rho_; // 要注意乘以 rho
        this->primalResidualInput_ = (this->uTinyMpcMat_ - this->zNewSlack_).cwiseAbs().maxCoeff();
        this->dualResidualInput_ = ((this->zSlackInputVar_ - this->zNewSlack_).cwiseAbs().maxCoeff()) * this->rho_;

        if (primalResidualState_ < this->absPriTol_ && primalResidualInput_ < absPriTol_ &&
            dualResidualState_ < this->absDuaTol_ && dualResidualInput_ < absDuaTol_) {
            return true;                 
        }
    }
    return false;
}

void Admm::BackwardPassGrad()
{
    // 从后往前更新
    for (int i = HORIZON - 2; i >= 0; --i) {
        (d.col(i)).noalias() = QuuInv_ * (discreteB_.transpose() * p.col(i + 1) 
                               + rLinearInputCost_.col(i));
        (p.col(i)).noalias() = qLinearStateCost_.col(i) + 
                               AmBKMatTranspose_.lazyProduct(p.col(i + 1)) - 
                               (K_k_.transpose()).lazyProduct(rLinearInputCost_.col(i)); 
        // std::cout << "p = " << p.col(i).format(MatrixApiFmt) << std::endl;
    }
}

int Admm::SolveTinyMpc()
{
    // Initialize variables
    this->solved_ = false;
    this->iter_ = 0;

    for (int i = 0; i < this->maxIter_; ++i) {
        // Solve linear system with Riccati and roll out to get new trajectory. u[i] x[i+1]
        ForwardPass();

        // Project slack variables into feasible domain. zNewSlack_ - u, vNewSlack_ - x
        UpdateSlack();

        // Compute next iteration of dual variables. 拉格朗日乘子:yDualInputVar_ - u, gDualStateVar_ - x
        UpdateDual();

        // Update linear control cost terms using reference trajectory, duals, and slack variables
        // rLinearInputCost_ - u, p qLinearStateCost_ - x
        UpdateLinearCost();

        iter_ += 1;

        // Check for whether cost is minimized by calculating residuals. 残差 < 阈值
        if (TerminationCondition()) {
            solved_ = true;
            // this->xTinyMpcMat_ = this->vNewSlack_;
            // this->uTinyMpcMat_ = this->zNewSlack_;
            return 0; 
        }

        // Save previous slack variables
        this->vSlackStateVar_ = this->vNewSlack_;
        this->zSlackInputVar_ = this->zNewSlack_;

        BackwardPassGrad();
        // std::cout<<"true"<<std::endl ;

    }
    solved_ = false;
    this->xTinyMpcMat_ = this->vNewSlack_;
    this->uTinyMpcMat_ = this->zNewSlack_;
    // std::cout<<"true"<<std::endl ;
    return 1;
}

// ADMM Update
StateVector Admm::GetXvecNext(const StateVector &xInit)
{
    this->xTinyMpcNextVec_ = discreteA_ * xInit + discreteB_ * this->uTinyMpcMat_.col(0);
    return xTinyMpcNextVec_; 
}

