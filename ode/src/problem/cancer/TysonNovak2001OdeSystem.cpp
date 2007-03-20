#include "TysonNovak2001OdeSystem.hpp"
#include <cmath>
#include <cassert>
#include "ReplicatableVector.hpp"
/**
 * Constructor
 */
TysonNovak2001OdeSystem::TysonNovak2001OdeSystem()
        : AbstractOdeSystemWithAnalyticJacobian(6)
{
    /*
     * State variables
     * 
     * These Initial conditions arethe steady state (approximate) 
     * solutions and the commented out conditions are from tyson 
     * novak paper we think.
     */
    mVariableNames.push_back("CycB");
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(0.1);
    
    mVariableNames.push_back("Cdh1");
    mVariableUnits.push_back("nM");
    //mInitialConditions.push_back(9.8770e-01);
    mInitialConditions.push_back(0.98913);
    
    mVariableNames.push_back("Cdc20T");
    mVariableUnits.push_back("nM");
    //mInitialConditions.push_back(1.5011e+00);
    mInitialConditions.push_back(1.54217);
    
    mVariableNames.push_back("Cdc20A");
    mVariableUnits.push_back("nM");
    //mInitialConditions.push_back(1.2924e+00);
    mInitialConditions.push_back(1.40563);
    
    mVariableNames.push_back("IEP");
    mVariableUnits.push_back("nM");
    //mInitialConditions.push_back(6.5405e-01);
    mInitialConditions.push_back(0.67083);
    
    mVariableNames.push_back("mass");
    mVariableUnits.push_back("");
    //mInitialConditions.push_back(4.7039e-01);
    mInitialConditions.push_back(0.95328/2.0);
    
    
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
 * @param rDy filled in with the resulting derivatives using Tyson & Novak system of equations
 */
void TysonNovak2001OdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
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
    // multiplying by 60 beacuase the paper has time in minutes wheras we have hours.
    rDY[0] = dx1*60.0;
    rDY[1] = dx2*60.0;
    rDY[2] = dx3*60.0;
    rDY[3] = dx4*60.0;
    rDY[4] = dx5*60.0;
    rDY[5] = dx6*60.0;
}



void TysonNovak2001OdeSystem::AnalyticJacobian(std::vector<double> &solutionGuess, double** jacobian, double time, double timeStep)
{

    timeStep *=60.0; // to scale Jacobian so in hours not minutes
    double x1 = solutionGuess[0];
    double x2 = solutionGuess[1];
    double x3 = solutionGuess[2];
    double x4 = solutionGuess[3];
    double x5 = solutionGuess[4];
    double x6 = solutionGuess[5];
    
    // f1
    double df1_dx1 = -mK2d - mK2dd*x2;
    double df1_dx2 = -mK2dd*x1;
    
    jacobian[0][0] =  1-timeStep*df1_dx1;
    jacobian[0][1] = -timeStep*df1_dx2;
    
    //f2
    double df2_dx1 = -mK4*x6*x2/(mJ4+x2);
    double df2_dx2 = -mJ3*(mK3d + mK3dd*x4)/(pow((mJ3 + 1 - x2),2))
                     -mJ4*mK4*x6*x1/(pow((mJ4+x2),2));
    double df2_dx4 =  mK3dd*(1-x2)/(mJ3+1-x2);
    double df2_dx6 = -mK4*x1*x2/(mJ4+x2);
    
    jacobian[1][0] = -timeStep*df2_dx1;
    jacobian[1][1] =  1-timeStep*df2_dx2;
    jacobian[1][3] = -timeStep*df2_dx4;
    jacobian[1][5] = -timeStep*df2_dx6;
    
    
    //f3
    double z = x1*x6/mJ5;
    double df3_dx1 = (mK5dd*x6/mJ5)*mN*pow(z,mN-1)/(pow((1-pow(z,mN)),2));
    double df3_dx3 = -mK6;
    double df3_dx6 = (mK5dd*x1/mJ5)*mN*pow(z,mN-1)/(pow((1-pow(z,mN)),2));
    
    jacobian[2][0] = -timeStep*df3_dx1;
    jacobian[2][2] = 1-timeStep*df3_dx3;
    jacobian[2][5] = -timeStep*df3_dx6;
    
    //f4
    double df4_dx3 =  mJ7*mK7*x5/(pow(mJ7+x3-x4,2));
    double df4_dx4 = -mJ7*mK7*x5/(pow(mJ7+x3-x4,2)) - mK6 - mJ8*mK8*mMad/(pow(mJ8+x4,2));
    double df4_dx5 =  mK7*(x3-x4)/(mJ7+x3-x4);
    
    jacobian[3][2] = -timeStep*df4_dx3;
    jacobian[3][3] = 1-timeStep*df4_dx4;
    jacobian[3][4] = -timeStep*df4_dx5;
    
    //f5
    double df5_dx1 =  mK9*x6*(1-x5);
    double df5_dx5 = -mK10 - mK9*x6*x1;
    double df5_dx6 =  mK9*x1*(1-x5);
    
    jacobian[4][0] = -timeStep*df5_dx1;
    jacobian[4][4] = 1-timeStep*df5_dx5;
    jacobian[4][5] = -timeStep*df5_dx6;
    
    
    //f6
    double df6_dx6 = mMu - 2*mMu*x6/mMstar;
    
    jacobian[5][5] = 1-timeStep*df6_dx6;
    
    
    
}



