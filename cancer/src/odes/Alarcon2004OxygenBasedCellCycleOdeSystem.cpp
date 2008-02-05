#include <cmath>
#include <cassert>
#include <vector>

#include "Alarcon2004OxygenBasedCellCycleOdeSystem.hpp"

/**
 * Constructor.
 *
 * @param oxygenConcentration is a non-dimensional oxygen concentration between 0 and 1.
 * @param rIsCancerCell affects the ODE system 
 */
Alarcon2004OxygenBasedCellCycleOdeSystem::Alarcon2004OxygenBasedCellCycleOdeSystem(double oxygenConcentration, const CellMutationState& rMutationState)
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
    
    assert(rMutationState == ALARCON_NORMAL || rMutationState == ALARCON_CANCER);
       
    mMutationState = rMutationState;
    
    // parameter values taken from the Alarcon et al. (2004) paper        
    if (mMutationState == ALARCON_NORMAL)    // normal cells
    {
        ma1 = 0.05;
        mc1 = 0.1;
        mxThreshold = 0.004;
        myThreshold = 0.2;
    }
    else // cancer cells
    {
        ma1 = 0.04;
        mc1 = 0.007;
        mxThreshold = 0.04; // should this be 0.004??
        myThreshold = 0.05;
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
void Alarcon2004OxygenBasedCellCycleOdeSystem::SetMutationState(const CellMutationState& rMutationState)
{
    mMutationState = rMutationState;
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
    
    assert(mMutationState == ALARCON_NORMAL || mMutationState == ALARCON_CANCER);
    
    // parameter values taken from the Alarcon et al. (2004) paper        
    if (mMutationState == ALARCON_NORMAL)    // normal cells
    {
        dz = mc1*(1 - mass/mMstar) - mc2*oxygen_concentration*z/(mB + oxygen_concentration);
    }
    else // cancer cells
    {
        dz = mc1 - mc2*oxygen_concentration*z/(mB + oxygen_concentration);   
    } 
        
    dmass = mEta*mass*(1 - mass/mMstar);        
    du = md1 - (md2 + md1*y)*u;

    // rescale time to be in hours
    rDY[0] = 60.0*dx;
    rDY[1] = 60.0*dy;
    rDY[2] = 60.0*dz;
    rDY[3] = 60.0*dmass;
    rDY[4] = 60.0*du;
    rDY[5] = 0.0; // do not change the oxygen concentration
}

