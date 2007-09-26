#include "GaryWntOdeSystem.hpp"
#include "CryptCellMutationStates.hpp"

#include <cmath>
#include <cassert>
#include <vector>

/**
 * Constructor.
 *
 * @param WntLevel is a non-dimensional Wnt value between 0 and 1. This sets up the Wnt pathway in its steady state.
 * @param mutationState affects the ODE system and is given by CryptCellMutationStates.hpp
 */
GaryWntOdeSystem::GaryWntOdeSystem(double WntLevel, const CryptCellMutationState& rMutationState)
        : AbstractOdeSystem(4)
{
    /*
     * State variables
    % 
    % 0. c = APC (Active)
    % 1. b1 = Beta-Catenin (1st allele's copy)
    % 2. b2 = Beta-Catenin (2nd allele's copy)
    % 3. WntLevel 
    */
    Init(); //Set up parameters
    
    double destruction_level = ma5d/(ma4d*WntLevel+ma5d);
    double beta_cat_level_1 = -1.0;
    double beta_cat_level_2 = -1.0;
    
    mMutationState = rMutationState;
    
    // These three lines set up a wnt signalling pathway in a steady state
    if (mMutationState == HEALTHY || mMutationState == LABELLED)	// healthy cells
    {
        beta_cat_level_1 = 0.5*ma2d/(ma2d+ma3d*destruction_level);
        beta_cat_level_2 = 0.5*ma2d/(ma2d+ma3d*destruction_level);
    }
    else if (mMutationState == APC_ONE_HIT) // APC +/-
    {
        beta_cat_level_1 = 0.5*ma2d/(ma2d+0.5*ma3d*destruction_level); // only half are active
        beta_cat_level_2 = 0.5*ma2d/(ma2d+0.5*ma3d*destruction_level);
    }
    else if (mMutationState == BETA_CATENIN_ONE_HIT) // Beta-cat delta 45
    {
        beta_cat_level_1 = 0.5*ma2d/(ma2d+ma3d*destruction_level);
        beta_cat_level_2 = 0.5;
    }
    else if (mMutationState == APC_TWO_HIT) // APC -/-
    {
        destruction_level = 0.0; // no active destruction complex
        beta_cat_level_1 = 0.5; // fully active beta-catenin
        beta_cat_level_2 = 0.5; // fully active beta-catenin
    }
    else
    {
        // can't get here until new mutation states are added to CryptCellMutationState
        #define COVERAGE_IGNORE
        assert(0);
        #undef COVERAGE_IGNORE
    }
    
    mVariableNames.push_back("APC");
    mVariableUnits.push_back("non_dim");
    mInitialConditions.push_back(destruction_level);
    
    mVariableNames.push_back("Beta_Cat1");
    mVariableUnits.push_back("non_dim");
    mInitialConditions.push_back(beta_cat_level_1);
    
    mVariableNames.push_back("Beta_Cat2");
    mVariableUnits.push_back("non_dim");
    mInitialConditions.push_back(beta_cat_level_2);
    
    mVariableNames.push_back("Wnt");
    mVariableUnits.push_back("non_dim");
    mInitialConditions.push_back(WntLevel);
}

/**
 * This should be called by the relevant cell cycle model before any solving
 * of the ODE system (as it is used to evaluate the Y derivatives).
 */
void GaryWntOdeSystem::SetMutationState(const CryptCellMutationState& rMutationState)
{
    mMutationState = rMutationState;
}


/**
 * Destructor
 */
GaryWntOdeSystem::~GaryWntOdeSystem(void)
{
    // Do nothing
}


void GaryWntOdeSystem::Init()
{
    double phi_E2F1 = 0.1;
    // Gary Parameters
    //double a1 = 0.423;
    double a2 = 2.57e-4;
    double a3 = 1.72;
    double a4 = 10.0;
    double a5 = 0.5;
    double WntMax = 10.0;

    //\todo change this without breaking the build
    double APC_Total = 0.02;
    
//  Non-dimensionalise...
//    mk2d = k2/(Km2*phi_E2F1);
//    mk3d = k3*a1*mitogenic_factorF/(Km4*phi_E2F1*a2);
//    mk34d = k34/phi_E2F1;
//    mk43d = k43/phi_E2F1;
//    mk23d = k23*Km2/(Km4*phi_E2F1);
//    mad = a/Km2;
//    mJ11d = J11*phi_E2F1/k1;
//    mJ12d = J12*phi_E2F1/k1;
//    mJ13d = J13*phi_E2F1/k1;
//    mJ61d = J61*phi_E2F1/k1;
//    mJ62d = J62*phi_E2F1/k1;
//    mJ63d = J63*phi_E2F1/k1;
//    mKm1d = Km1/Km2;
//    mkpd = kp/(Km2*phi_E2F1);
//    mphi_r = phi_pRb/phi_E2F1;
//    mphi_i = phi_CycDi/phi_E2F1;
//    mphi_j = phi_CycDa/phi_E2F1;
//    mphi_p = phi_pRbp/phi_E2F1;
    ma2d = a2/phi_E2F1;
    ma3d = a3*APC_Total/phi_E2F1;
    ma4d = a4*WntMax/phi_E2F1;
    ma5d = a5/phi_E2F1;
//    mk16d = k16*Km4/phi_E2F1;
//    mk61d = k61/phi_E2F1;
    mPhiE2F1 = phi_E2F1;
}

/**
 * Returns a vector representing the RHS of the odes at each time step, y' = [y1' ... yn'].
 * Some ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
 *
 * @param rDY filled in with the resulting derivatives (using Mirams et al. (2007?) system of equations)
 */
void GaryWntOdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    double c = rY[0];
    double b1 = rY[1];
    double b2 = rY[2];
    double WntLevel = rY[3];
    
    double dx6 = 0.0;
    double dx7 = 0.0;
    double dx8 = 0.0;
    
    /*
    % The variables are
    % 1. c = APC (Active)
    % 2. b1 = Beta-Catenin allele 1
    % 3. b2 = Beta-Catenin allele 2
    % 4. Wnt level
    */
    
    
    // Mutations take effect by altering the level of beta-catenin
    if (mMutationState==HEALTHY || mMutationState==LABELLED)	// HEALTHY CELL
    {
        // da
        dx6 = ma5d*(1.0-c) - ma4d*WntLevel*c;
        //db
        dx7 = ma2d*(0.5-b1) - ma3d*b1*c;
        dx8 = ma2d*(0.5-b2) - ma3d*b2*c;
    }
    else if (mMutationState==APC_ONE_HIT) // APC +/-
    {
        dx6 = ma5d*(1.0-c) - ma4d*WntLevel*c;
        dx7 = ma2d*(0.5-b1) - 0.5*ma3d*b1*c;
        dx8 = ma2d*(0.5-b2) - 0.5*ma3d*b2*c;
    }
    else if (mMutationState==BETA_CATENIN_ONE_HIT) // Beta-Cat D45
    {
        dx6 = ma5d*(1.0-c) - ma4d*WntLevel*c;
        dx7 = ma2d*(0.5-b1) - ma3d*b1*c;
        dx8 = ma2d*(0.5-b2);
    }
    else if (mMutationState==APC_TWO_HIT) // APC -/-
    {
        dx6 = 0.0;
        dx7 = ma2d*(0.5-b1);
        dx8 = ma2d*(0.5-b2);
    }
    else
    {
        // can't get here until new mutation states are added to CryptCellMutationState
        #define COVERAGE_IGNORE
        assert(0);
        #undef COVERAGE_IGNORE
    }
   
    double factor = mPhiE2F1*60.0;  // Convert non-dimensional d/dt s to d/dt in hours.
    

    rDY[0] = dx6*factor;
    rDY[1] = dx7*factor; // beta-cat allele 1
    rDY[2] = dx8*factor; // beta-cat allele 2
    rDY[3] = 0.0; // Do not change the Wnt level.
}

