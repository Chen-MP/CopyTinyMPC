#include "../include/matplotlibcpp/matplotlibcpp.h" 

#include "../data_store/pendulum_data.h"
#include "../tinympc/discretization.h"
#include "../tinympc/lqr_dp.h"
#include "../tinympc/admm.h"

int main()
{
    StateMatrix continuousA = Eigen::Map<Eigen::Matrix<double, NUM_STATES, NUM_STATES, Eigen::RowMajor>>(A_data);
    InputMatrix continuousB = Eigen::Map<Eigen::Matrix<double, NUM_STATES, NUM_INPUTS>>(B_data);

    double Q_data[NUM_STATES] = {10, 10, 1, 1};
    double R_data[NUM_INPUTS] = {1};
    double xInit_data[NUM_STATES] = {0.3, 3.14/40.0, 0, 0}; // x, q, \dot{x}, \dot{q}
    double rho = 100.0; 
    double sampleTime = 0.01;
    bool verbose = true;

    StateVector Q_diag = Eigen::Map<StateVector>(Q_data);
    InputVector R_diag = Eigen::Map<InputVector>(R_data);

    StateVector xInit = Eigen::Map<StateVector>(xInit_data);

    StateVector xMax;
    xMax << 10.0, 9.0, 8.0, 7.0;
    StateVector xMin;
    xMin << -3.0, -2.0, -2.1, -0.22;

    InputVector uMax;
    uMax << 10.0;
    InputVector uMin;
    uMin << -10.0;

    std::vector<double> x;
    std::vector<double> q;
    std::vector<double> xDot;
    std::vector<double> qDot;
    std::vector<double> u;
    
    Discretization DiscreteSystem(sampleTime);
    DiscreteSystem.FowardEulerDiscretization(continuousA, continuousB);
    StateMatrix discreteA = DiscreteSystem.GetDiscreteA();
    InputMatrix discreteB = DiscreteSystem.GetDiscreteB();
    
    Admm AdmmAlg(discreteA, discreteB, Q_diag, R_diag, rho);
    
    AdmmAlg.GetLqrDpParams();

    AdmmAlg.SetXrefVecMat();
    AdmmAlg.SetUrefVecMat();
    /************** Disable Constrain *******************/
    // AdmmAlg.DisableStateConstrain();
    // AdmmAlg.DisableInputConstrain();

    /********** Only Enable Input Constrain *************/
    // AdmmAlg.DisableStateConstrain(); // Only Enable State Constrain is not use
    // AdmmAlg.EnableInputConstrain(uMax, uMin);

    /********** Enable State && Input Constrain *************/
    // TODO： TinyMpc Enable State Constrain 失败
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


