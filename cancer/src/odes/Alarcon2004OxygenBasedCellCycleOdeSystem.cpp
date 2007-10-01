#include "Alarcon2004OxygenBasedCellCycleOdeSystem.hpp"

#include <cmath>
#include <cassert>
#include <vector>

/**
 * Constructor.
 *
 * @param oxygenConcentration is a non-dimensional oxygen concentration between 0 and 1.
 * @param mutationState affects the ODE system and is given by CryptCellMutationStates.hpp
 */
Alarcon2004OxygenBasedCellCycleOdeSystem::Alarcon2004OxygenBasedCellCycleOdeSystem(double oxygenConcentration, const bool &rIsCancerCell)
        : AbstractOdeSystem(6)
{
    /*
     * State variables
     * 
     % 0. x = Cdh1-APC complexes
     % 1. y = cyclin-CDK
     % 2. z = p27
     % 3. m = mass
     % 4. u = RBNP
     % 5. oxygenConcentration
    */
    Init(); // set up parameters
    
    mIsCancerCell = rIsCancerCell;

    // parameter values taken from the Alarcon et al. (2004) paper
    if (mIsCancerCell)
    {           
        ma1 = 0.04;
        mc1 = 0.007;
        mxThreshold = 0.04; // should this be 0.004??
        myThreshold = 0.05;
    }
    else
    {
        ma1 = 0.05;
        mc1 = 0.1;
        mxThreshold = 0.004;
        myThreshold = 0.2;
    }  

    mVariableNames.push_back("Cdh1_APC_complexes");
    mVariableUnits.push_back("non_dim");
    mInitialConditions.push_back(0.9);
    
    mVariableNames.push_back("cyclin_CDK");
    mVariableUnits.push_back("non_dim");
    mInitialConditions.push_back(0.01);
    
    mVariableNames.push_back("p27");
    mVariableUnits.push_back("non_dim");
    mInitialConditions.push_back(0.0);    
    
    mVariableNames.push_back("mass");
    mVariableUnits.push_back("non_dim");
    mInitialConditions.push_back(mMstar/2.0);
    
    mVariableNames.push_back("RBNP");
    mVariableUnits.push_back("non_dim");
    mInitialConditions.push_back(1.0);
        
    mVariableNames.push_back("O2");
    mVariableUnits.push_back("non_dim");
    mInitialConditions.push_back(oxygenConcentration);
}

/**
 * This should be called by the relevant cell cycle model before any solving
 * of the ODE system (as it is used to evaluate the Y derivatives).
 */
void Alarcon2004OxygenBasedCellCycleOdeSystem::SetIsCancerCell(const bool& rIsCancerCell)
{
    mIsCancerCell = rIsCancerCell;
}

/**
 * Destructor
 */
Alarcon2004OxygenBasedCellCycleOdeSystem::~Alarcon2004OxygenBasedCellCycleOdeSystem(void)
{
    // Do nothing
}

void Alarcon2004OxygenBasedCellCycleOdeSystem::Init()
{
    // parameters values taken from the Alarcon et al. (2004) paper
    ma2 = 1.0;
    ma3 = 0.25;
    ma4 = 0.04;
    mb3 = 10.0;
    mb4 = 5.5;
    mc2 = 0.01;
    md1 = 0.01;
    md2 = 0.1;
    mJ3 = 0.04;
    mJ4 = 0.04;
    mEta = 0.01;
    mMstar = 10.0;
    mB = 0.01;
}

/**
 * Returns a vector representing the RHS of the odes at each time step, y' = [y1' ... yn'].
 * Some ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
 *
 * @param rDY filled in with the resulting derivatives (using Alarcons et al. (2004) system of equations)
 */
void Alarcon2004OxygenBasedCellCycleOdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    double x = rY[0];
    double y = rY[1];
    double z = rY[2];
    double mass = rY[3];
    double u = rY[4];
    double oxygen_concentration = rY[5];
    
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    double dmass = 0.0;
    double du = 0.0;
    
    /*
     % The variables are
     % 1. x = Cdh1-APC complexes
     % 2. y = cyclin-CDK
     % 3. z = p27
     % 4. m = mass
     % 5. u = RBNP
    */

    dx = ((1 + mb3*u)*(1-x))/(mJ3 + 1 - x) - (mb4*mass*x*y)/(mJ4 + x);
    
    dy = ma4 -(ma1 + ma2*x + ma3*z)*y;
    
    if (mIsCancerCell)
    {
        dz = mc1 - mc2*oxygen_concentration*z/(mB + oxygen_concentration);   
    }
    else
    {
        dz = mc1*(1 - mass/mMstar) - mc2*oxygen_concentration*z/(mB + oxygen_concentration);   
    }
    
    dmass = mEta*mass*(1 - mass/mMstar);
        
    du = md1 - (md2 + md1*y)*u;

    rDY[0] = dx;
    rDY[1] = dy;
    rDY[2] = dz;
    rDY[3] = dmass;
    rDY[4] = du;
    rDY[5] = 0.0; // do not change the oxygen concentration
}

