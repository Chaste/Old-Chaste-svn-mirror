#include "TysonNovak2001OdeSystem.hpp"
#include <cmath>
/**
 * Constructor
 */
TysonNovak2001OdeSystem::TysonNovak2001OdeSystem() : AbstractOdeSystem()
{
    /*
     * State variable
     */
    mVariableNames.push_back("CycB");
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(0.1-2.0*1.0e-2);
    
    mVariableNames.push_back("Cdh1");
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(9.8770e-01);
    
    mVariableNames.push_back("Cdc20T");
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(1.5011e+00);
    
    mVariableNames.push_back("Cdc20A");
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(1.2924e+00);
    
    mVariableNames.push_back("IEP");
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(6.5405e-01);
    
    mVariableNames.push_back("mass");
    mVariableUnits.push_back("");
    mInitialConditions.push_back(4.7039e-01);
    
    mNumberOfStateVariables=6;
    
    Init();
}


/**
 * Destructor
 */
TysonNovak2001OdeSystem::~TysonNovak2001OdeSystem(void)
{
    // Do nothing
}


void TysonNovak2001OdeSystem::Init()
{
    // Initialize model parameters
    mK1 = 0.04;
    mK2d = 0.04;
    mK2dd = 1.0;
    mK2ddd = 1.0;
    mCycB_threshold = 0.1;
    mK3d = 1.0;
    mK3dd = 10.0;
    mK4d = 2.0;
    mK4 = 35;
    mJ3 = 0.04;
    mJ4 = 0.04;
    mK5d = 0.005;
    mK5dd = 0.2;
    mK6 = 0.1;
    mJ5 = 0.3;
    mN = 4;
    mK7 = 1.0;
    mK8 = 0.5;
    mJ7 = 1e-3;
    mJ8 = 1e-3;
    mMad = 1.0;
    mK9 = 0.1;
    mK10 = 0.02;
    mK11 = 1.0;
    mK12d = 0.2;
    mK12dd = 50.0;
    mK12ddd = 100.0;
    mKeq = 1e3;
    mK13 = 1.0;
    mK14 = 1.0;
    mK15d = 1.5;
    mK15dd = 0.05;
    mK16d = 1.0;
    mK16dd = 3.0;
    mJ15 = 0.01;
    mJ16 = 0.01;
    mMu = 0.01;
    mMstar = 10.0;
}

/**
 * Returns a vector representing the RHS of the Tyson & Novak system of Odes at each time step, y' = [y1' ... yn'].
 * Some ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
 *
 * @return std::vector<double> RHS of Tyson & Novak system of equations
 */
std::vector<double> TysonNovak2001OdeSystem::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
    double x1 = rY[0];
    double x2 = rY[1];
    double x3 = rY[2];
    double x4 = rY[3];
    double x5 = rY[4];
    double x6 = rY[5];
    
    double dx1 = 0.0;
    double dx2 = 0.0;
    double dx3 = 0.0;
    double dx4 = 0.0;
    double dx5 = 0.0;
    double dx6 = 0.0;
    double temp1 = 0.0;
    double temp2 = 0.0;
    double temp3 = 0.0;
    
    /**
    % 1. [CycB]
    % 2. [Cdh1]
    % 3. [Cdc20T]
    % 4. [Cdc20A]
    % 5. [IEP]
    % 6. m - mass of the cell
    */
    
    std::vector<double> RHS;
    
    dx1 = mK1-(mK2d+mK2dd*x2)*x1;
    
    // This line just models the start transition, no cycling, without Cdc20A
    //temp1 = ((mK3d)*(1.0-x2))/(mJ3+1.0-x2);
    
    temp1 = ((mK3d+mK3dd*x4)*(1.0-x2))/(mJ3+1.0-x2);
    temp2 = (mK4*x6*x1*x2)/(mJ4+x2);
    dx2 = temp1-temp2;
    
    temp1 = mK5dd*(pow(x1*x6/mJ5,mN)/(1+pow(x1*x6/mJ5,mN)));
    temp2 = mK6*x3;
    dx3 = mK5d + temp1 - temp2;
    
    temp1 = (mK7*x5*(x3-x4))/(mJ7+x3-x4);
    temp2 = (mK8*mMad*x4)/(mJ8+x4);
    temp3 = mK6*x4;
    dx4 = temp1 - temp2 - temp3;
    
    dx5 = mK9*x6*x1*(1.0-x5) - mK10*x5;
    
    dx6 = mMu*x6*(1.0-x6/mMstar);
    
    RHS.push_back(dx1);
    RHS.push_back(dx2);
    RHS.push_back(dx3);
    RHS.push_back(dx4);
    RHS.push_back(dx5);
    RHS.push_back(dx6);
    
    return RHS;
}


