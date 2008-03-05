#include "Lee2003WntSignallingOdeSystem.hpp"
#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>

#define COVERAGE_IGNORE
/**
 * Constructor.
 *
 * @param WntLevel is a non-dimensional Wnt value between 0 and 1. This sets up the Wnt pathway in its steady state.
 * 
 */
Lee2003WntSignallingOdeSystem::Lee2003WntSignallingOdeSystem(double WntLevel)
        : AbstractOdeSystem(8)
{
    //
    // State variables
    // 
    // % The variables are
    // % 0. X2 Dsh_active
    // % 1. X3 APC*/axin*/GSK3
    // % 2. X4 APC/axin/GSK3
    // % 3. X9 beta-cat*/APC*/axin*/GSK3
    // % 4. X10 beta-cat*
    // % 5. X11 beta-cat
    // % 6. X12 axin
    // % 7. WntLevel
    //
    Init(); //Set up parameters
    // unstimulated state first of all...
    
    mVariableNames.push_back("Dsh_active");
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(0.0);
    
    mVariableNames.push_back("APC_axin_GSK3");
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(4.83e-3);
    
    mVariableNames.push_back("beta_cat_P_APC_P_axin_P_GSK3");
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(2.02e-3);
    
    mVariableNames.push_back("beta_cat_P");
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(1.00);
    
    mVariableNames.push_back("APC_P_axin_P_GSK3");
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(9.66e-3);
    
    mVariableNames.push_back("beta_cat");
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(25.1);
    
    mVariableNames.push_back("axin");
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(4.93e-4);
    
    mVariableNames.push_back("Wnt");
    mVariableUnits.push_back("non_dim");
    mInitialConditions.push_back(WntLevel);
    
    mNumberOfStateVariables=8;
}


/**
 * Destructor
 */
Lee2003WntSignallingOdeSystem::~Lee2003WntSignallingOdeSystem(void)
{
    // Do nothing
}


void Lee2003WntSignallingOdeSystem::Init()
{
    // Initialize model parameters
    // Lee (2003) Parameters
    mDsh0 = 100.0;
    mAPC0 = 100.0;
    mTCF0 = 15.0;
    mGSK0 = 50.0;
    mK7 = 50.0;
    mK8 = 120.0;
    mK16 = 30.0;
    mK17 = 1200.0;
    mk1 = 0.182;
    mk2 = 1.82e-2;
    mk3 = 5.0e-2;
    mk4 = 0.267;
    mk5 = 0.133;
    mk6 = 9.09e-2;
    mk_6 = 0.909;
    mk9 = 206.0;
    mk10 = 206.0;
    mk11 = 0.417;
    mk13 = 2.57e-4;
    mk15 = 0.167;
    mv12 = 0.423;
    mv14 = 8.22e-5;
}

/**
 * Returns a vector representing the RHS of the odes at each time step, y' = [y1' ... yn'].
 * Some ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
 *
 * @param rDY filled in with the resulting derivatives (using Lee et al. (2003) system of equations)
 */
void Lee2003WntSignallingOdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    double X5 = mGSK0;
    double X2 = rY[0];
    double X4 = rY[1];
    double X9 = rY[2];
    double X10 = rY[3];
    double X3 = rY[4];
    double X11 = rY[5];
    double X12 = rY[6];
    double WntLevel = rY[7];  

    for(unsigned i=0 ; i<rY.size() ; i++)
    {
        assert( rY[i] >= 0.0 ); // all protein concentrations are positive...
    }

    // Easy Ones - A.32, A.37, A.34, A.35
    double dX2 = mk1*WntLevel*(mDsh0-X2) - mk2*X2;
    double dX4 = -(mk3*X2+mk4+mk_6)*X4 + mk5*X3 + mk6*X5*((mK17*X12*mAPC0)/(mK7*(mK17+X11)));
    double dX9 = (mk9*X3*X11)/mK8 - mk10*X9;
    double dX10 = mk10*X9 - mk11*X10;
    
    // Bit of rearranging of A.43 and A.44 simultaneous equations
    double a = 1 + X11/mK8;
    double b = X3/mK8;
    double c = X11/mK8;
    double d = 1.0 + X3/mK8 + (mTCF0*mK16)/((mK16+X11)*(mK16+X11)) + (mAPC0*mK17)/((mK17+X11)*(mK17+X11));
    double e = mk4*X4 - mk5*X3 - (mk9*X3*X11)/mK8 + mk10*X9;
    double f = mv12 - ((mk9*X3)/mK8 +mk13)*X11;
    double dX3 = (e-(b/d)*f)/(a-(b/d)*c);
    double dX11 = (e-(a/c)*f)/(b-(a/c)*d);
    
    // And a bit more for A.42
    double temp1 = mk3*X2*X4 - (mk6*mGSK0*mAPC0*mK17*X12)/(mK7*(mK17+X11)) + mk_6*X4+mv14-mk15*X12;
    double temp2 = (dX11*mAPC0*mK17*X12)/(mK7*(mK17+X11)*(mK17+X11));
    double temp3 = 1 + (mAPC0*mK17)/(mK7*(mK17+X11));
    double dX12 = (temp1 + temp2)/temp3;
    
    double factor = 60.0;  // Convert d/dt in minutes to d/dt in hours.
        
    rDY[0] = dX2*factor;
    rDY[1] = dX4*factor;
    rDY[2] = dX9*factor;
    rDY[3] = dX10*factor;
    rDY[4] = dX3*factor;
    rDY[5] = dX11*factor;
    rDY[6] = dX12*factor; 
    rDY[7] = 0.0; // Do not change the Wnt level.
}

#undef COVERAGE_IGNORE

