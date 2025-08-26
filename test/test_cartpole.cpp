#include "../include/matplotlibcpp/matplotlibcpp.h" 

#include "../tinympc/discretization.h"
#include "../tinympc/lqr_dp.h"
#include "../tinympc/admm.h"

int main()
{
    float rho_value = 1.0; // penalty term

    // states - discrete
    double A_data[NUM_STATES * NUM_STATES] = {1.0, 0.01, 0.0, 0.0, 
                                             0.0, 1.0, 0.039, 0.0, 
                                             0.0, 0.0, 1.002, 0.01, 
                                             0.0, 0.0, 0.458, 1.002};
                                             
    double B_data[NUM_STATES * NUM_INPUTS] = {0.0, 0.02, 0.0, 0.067};
    // weight
    double Q_data[NUM_STATES] = {10.0, 1.0, 10.0, 1.0};
    double R_data[NUM_INPUTS] = {1.0};

    StateMatrix discreteA = Eigen::Map<Eigen::Matrix<double, NUM_STATES, NUM_STATES, Eigen::RowMajor>>(A_data);
    InputMatrix discreteB = Eigen::Map<Eigen::Matrix<double, NUM_STATES, NUM_INPUTS>>(B_data);

    double xInit_data[NUM_STATES] = {0.5, 0.0, 0.0, 0.0}; 
    double rho = 1000.0; 
    // double sampleTime = 0.01;
    bool verbose = true;

    StateVector Q_diag = Eigen::Map<StateVector>(Q_data);
    InputVector R_diag = Eigen::Map<InputVector>(R_data);

    StateVector xInit = Eigen::Map<StateVector>(xInit_data);

    StateVector xMax;
    xMax << 0.0, 0.0, 0.0, 10.0;
    StateVector xMin;
    xMin << -10.0, -10.0, -0.028, -0.0;

    InputVector uMax;
    uMax << 10.0;
    InputVector uMin;
    uMin << -10.0;

    std::vector<double> x;
    std::vector<double> q;
    std::vector<double> xDot;
    std::vector<double> qDot;
    std::vector<double> u;
    
    // Discretization DiscreteSystem(sampleTime);
    // DiscreteSystem.FowardEulerDiscretization(continuousA, continuousB);
    // StateMatrix discreteA = DiscreteSystem.GetDiscreteA();
    // InputMatrix discreteB = DiscreteSystem.GetDiscreteB();
    
    Admm AdmmAlg(discreteA, discreteB, Q_diag, R_diag, rho);
    
    AdmmAlg.GetLqrDpParams();
    
    // AdmmAlg.SetXrefVecMat();
    StateVector xRef;
    xRef << 1.0, 0, 0, 0 ;
    AdmmAlg.SetXrefVecMat(xRef); // 
    
    AdmmAlg.SetUrefVecMat();

    /************** Disable Constrain *******************/
    AdmmAlg.DisableStateConstrain();
    AdmmAlg.DisableInputConstrain();

    /********** Only Enable Input Constrain *************/
    // AdmmAlg.DisableStateConstrain(); // Only Enable State Constrain is not use
    // AdmmAlg.EnableInputConstrain(uMax, uMin);

    /********** Enable State && Input Constrain *************/
    // TODO： TinyMpc Enable State Constrain fail
    AdmmAlg.EnableStateConstrain(xMax, xMin);
    AdmmAlg.EnableInputConstrain(uMax, uMin);

    
    for (int i = 0; i < 500; ++i) {
        AdmmAlg.SetXmatTinyMpcInit(xInit);
        // AdmmAlg.SetNewStateConstrain(xInit);
        x.push_back(xInit(0));
        q.push_back(xInit(1));
        xDot.push_back(xInit(2));
        qDot.push_back(xInit(3));
        AdmmAlg.SolveTinyMpc();
        StateVector xNext = AdmmAlg.GetXvecNext(xInit);
        // StateVector xNext = AdmmAlg.GetNextState(xRef, xInit);
        u.push_back(AdmmAlg.GetUvec()(0));
        // std::cout << "u = " << AdmmAlg.GetUvec()(0) << std::endl;
        xInit = xNext;
    }

    matplotlibcpp::subplot(5, 1, 1);
    matplotlibcpp::plot(x);
    matplotlibcpp::ylabel("x");
    matplotlibcpp::subplot(5, 1, 2);
    matplotlibcpp::plot(q);
    matplotlibcpp::ylabel("q");
    matplotlibcpp::subplot(5, 1, 3);
    matplotlibcpp::plot(xDot);
    matplotlibcpp::ylabel("xDot");
    matplotlibcpp::subplot(5, 1, 4);
    matplotlibcpp::plot(qDot);
    matplotlibcpp::ylabel("qDot");
    matplotlibcpp::subplot(5, 1, 5);
    matplotlibcpp::plot(u);
    matplotlibcpp::ylabel("u");
    matplotlibcpp::show();
    return 0;
}
