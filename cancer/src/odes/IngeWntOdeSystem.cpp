#include "IngeWntOdeSystem.hpp"
#include "CellMutationStates.hpp"

#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>

/**
 * Constructor.
 *
 * @param WntLevel is a non-dimensional Wnt value between 0 and 1. This sets up the Wnt pathway in its steady state.
 * @param mutationState affects the ODE system and is given by CellMutationStates.hpp
 */
IngeWntOdeSystem::IngeWntOdeSystem(double wntLevel, const CellMutationState& rMutationState)
        : AbstractOdeSystem(17),
          mMutationState(rMutationState)
        
{
    Init(); //Set up parameters
    
    double d_d_hat = mDd + mXiD*wntLevel;
    double d_d_x_hat = mDdx + mXiDx*wntLevel;
    double d_x_hat = mDx + mXiX*wntLevel;
    double p_c_hat = mPc + mXiC*wntLevel;
    
    double sigma_D = 0.0;   // for healthy cells
    double sigma_B = 0.0;   // for healthy cells
    
    switch(mMutationState)
    {
        case HEALTHY:
        {
            break;
        }   
        case LABELLED:
        {
            break;
        }
        case APC_ONE_HIT:
        {
            sigma_D = 0.5;
            break;
        }
        case APC_TWO_HIT:
        {
            sigma_D = 1.0;
            break;   
        }
        case BETA_CATENIN_ONE_HIT:
        {
            sigma_B = 0.5;
            break;
        }
        default:
        #define COVERAGE_IGNORE
            assert(0);  // this can't happen if all mutation states are catered for.
        #undef COVERAGE_IGNORE
    }
    
    double steady_D = ((1.0-sigma_D)*mSd*mSx)/((1.0-sigma_D)*mSd*d_d_hat + d_x_hat*(d_d_hat + d_d_x_hat));
    
    mVariableNames.push_back("D");  //  Destruction complex (APC/Axin/GSK3B)
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(steady_D);
    
    double temp = (mSx*(d_d_hat+d_d_x_hat))/((1.0-sigma_D)*mSd*d_d_hat+d_x_hat*(d_d_hat+d_d_x_hat));
    
    mVariableNames.push_back("X");  //  Axin
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(temp);
    
    double steady_Cf = ((mSc-mDc*mKd - mPu*steady_D)+sqrt(pow((mSc-mDc*mKd - mPu*steady_D),2) + (4.0*mSc*mDc*mKd)))/(2.0*mDc);
    temp = (mPu*steady_D*steady_Cf)/(mDu*(steady_Cf+mKd));
    
    mVariableNames.push_back("Cu"); //  beta-catenin to be ubiquitinated
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(temp);
    
    double theta = mDc+ (mPu*steady_D)/(steady_Cf + mKd);
    
    double steady_Co = ( mSc - p_c_hat - theta*mKc + sqrt(4.0*mSc*theta*mKc + pow((mSc - p_c_hat - theta*mKc),2)) )/(2.0*theta);
    
    mVariableNames.push_back("Co"); //  Open form beta-catenin
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(steady_Co);
    
    double steady_Cc = steady_Cf - steady_Co;
    
    mVariableNames.push_back("Cc"); //  Closed form beta-catenin
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(steady_Cc);
    
    mVariableNames.push_back("Mo"); //  Open form mutant beta-catenin
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(0.0);
    
    mVariableNames.push_back("Mc"); //  Closed form mutant beta-catenin
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(0.0);
    
    mVariableNames.push_back("A");  //  `Free' adhesion molecules
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(mSa/mDa);
    
    mVariableNames.push_back("Ca"); //  Co-A    Adhesion complex
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(mSa*mSca*steady_Co/(mDa*mDca));
    
    mVariableNames.push_back("Ma"); //  Mo-A    Mutant adhesion complex
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(0.0);
    
    mVariableNames.push_back("T");  //  `free' transcription molecules (TCF)
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(mSt/mDt);
    
    mVariableNames.push_back("Cot");//  Co-T open form beta-catenin/TCF
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(mSct*mSt*steady_Co/(mDt*mDct));
    
    mVariableNames.push_back("Cct");//  Cc-T closed beta-catenin/TCF
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(mSct*mSt*steady_Cc/(mDt*mDct));
    
    mVariableNames.push_back("Mot");//  Mo-T open form mutant beta-catenin/TCF
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(0.0);
    
    mVariableNames.push_back("Mct");//  Mc-T closed form mutant beta-catenin/TCF
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(0.0);
    
    temp = (mSct*mSt*mSy*steady_Cf)/(mDy*(mSct*mSt*steady_Cf + mDct*mDt*mKt));
    
    mVariableNames.push_back("Y");  //  Wnt target protein
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(temp);
    
    mVariableNames.push_back("Sw");  //  Wnt stimulus
    mVariableUnits.push_back("nM");
    mInitialConditions.push_back(wntLevel);
}

/**
 * This should be called by the relevant cell cycle model before any solving
 * of the ODE system (as it is used to evaluate the Y derivatives).
 */
void IngeWntOdeSystem::SetMutationState(const CellMutationState& rMutationState)
{
    mMutationState = rMutationState;
}

/** 
 * Called by the archive function on the wnt cell cycle model.
 * @return mMutationState the mutation state of the cell defined by 
 * CellMutationStates.hpp
 */
CellMutationState& IngeWntOdeSystem::rGetMutationState()
{
    return mMutationState;
}

/**
 * Destructor
 */
IngeWntOdeSystem::~IngeWntOdeSystem(void)
{
    // Do nothing
}


void IngeWntOdeSystem::Init()
{
    // Initialize model parameters
    mSa = 20;   //  nM/h
    mSca = 250; //  (nMh)^-1
    mSc = 25;   //  nM/h
    mSct = 30;  //  (nMh)^-1
    mSd = 100;  //  h^-1
    mSt = 10;   //  nM/h
    mSx = 10;   //  nM/h
    mSy = 10;   //  h^-1
    mDa = 2;    //  h^-1
    mDca = 350; //  h^-1
    mDc = 1;    //  h^-1
    mDct = 750; //  h^-1
    mDd = 5;    //  h^-1
    mDdx = 5;   //  h^-1
    mDt = 0.4;  //  h^-1
    mDu = 50;   //  h^-1
    mDx = 100;  //  h^-1
    mDy = 1;    //  h^-1
    mKc = 200;  //  nM
    mKd = 5;    //  nM
    mKt = 50;   //  nM
    mPc = 0.0;  //  h^-1
    mPu = 100;  //  h^-1
    mXiD = 5;   //  h^-1
    mXiDx = 5;  //  h^-1
    mXiX = 200; //  h^-1
    mXiC = 0.0; //  h^-1 (FOR HYPOTHESIS ONE)
}

/**
 * Returns a vector representing the RHS of the odes at each time step, y' = [y1' ... yn'].
 * Some ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
 *
 * @param rDY filled in with the resulting derivatives (using Mirams et al. (2007?) system of equations)
 */
void IngeWntOdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    // variables
    double D = rY[0];
    double X = rY[1];
    double Cu = rY[2];
    double Co = rY[3];
    double Cc = rY[4];
    double Mo = rY[5];
    double Mc = rY[6];
    double A = rY[7];
    double Ca = rY[8];
    double Ma = rY[9];
    double T = rY[10];
    double Cot = rY[11];
    double Cct = rY[12];
    double Mot = rY[13];
    double Mct = rY[14];
    double Y = rY[15];
    double stimulus_wnt = rY[16];
    
    // Totals
    double Cf = Cc+Co;
    double Ct = Cct+Cot;
    double Mf = Mc+Mo;
    double Mt = Mct+Mot;
    
    double d_d_hat = mDd + mXiD*stimulus_wnt;
    double d_d_x_hat = mDdx + mXiDx*stimulus_wnt;
    double d_x_hat = mDx + mXiX*stimulus_wnt;
    double p_c_hat = mPc + mXiC*stimulus_wnt;
    
    double sigma_D = 0.0;   // for healthy cells
    double sigma_B = 0.0;   // for healthy cells
    
    switch(mMutationState)
    {
        case HEALTHY:
        {
            break;
        }   
        case LABELLED:
        {
            break;
        }
        case APC_ONE_HIT:
        {
            sigma_D = 0.5;
            break;
        }
        case APC_TWO_HIT:
        {
            sigma_D = 1.0;
            break;   
        }
        case BETA_CATENIN_ONE_HIT:
        {
            sigma_B = 0.5;
            break;
        }
        default:
        #define COVERAGE_IGNORE
            assert(0);  // this can't happen if all mutation states are catered for.
        #undef COVERAGE_IGNORE
    }
      
    rDY[0] = (1.0-sigma_D)*mSd*X - (d_d_hat + d_d_x_hat)*D;
    rDY[1] = mSx - (1.0-sigma_D)*mSd*X - d_x_hat*X + d_d_x_hat*D;
    rDY[2] = (mPu*D*Cf)/(Cf+mKd) - mDu*Cu;

    rDY[3] = (1.0-sigma_B)*mSc + mDca*Ca + mDct*Cot - (mSca*A + mSct*T + mDc)*Co
             - (p_c_hat*Co)/(Co + Mo + mKc) - (mPu*D*Co)/(Cf+mKd);
    

    rDY[4] = (p_c_hat*Co)/(Co + Mo + mKc) + mDct*Cct - (mSct*T + mDc)*Cc
             - (mPu*D*Cc)/(Cf+mKd);
    rDY[5] = sigma_B*mSc + mDca*Ma + mDct*Mot - (mSca*A + mSct*T + mDc)*Mo
             - (p_c_hat*Mo)/(Co + Mo + mKc);
    rDY[6] = (p_c_hat*Mo)/(Co + Mo + mKc) + mDct*Mct - (mSct*T + mDc)*Mc;    
    rDY[7] = mSa + mDca*(Ca+Ma) - (mSca*(Co+Mo) + mDa)*A;
    rDY[8] = mSca*Co*A - mDca*Ca; 

    rDY[9] = mSca*Mo*A - mDca*Ma; 

    rDY[10] = mSt + mDct*(Ct+Mt) - mSct*(Cf+Mf)*T - mDt*T;
    rDY[11] = mSct*Co*T - mDct*Cot; 
    rDY[12] = mSct*Cc*T - mDct*Cct; 
    rDY[13] = mSct*Mo*T - mDct*Mot; 
    rDY[14] = mSct*Mc*T - mDct*Mct; 
    rDY[15] = (mSy*(Ct+Mt))/(Ct + Mt + mKt) - mDy*Y;
    rDY[16] = 0.0;  // don't interfere with Wnt stimulus.
}

void IngeWntOdeSystem::SetUseHypothesisTwo(bool hypothesisTwo)
{
    if (hypothesisTwo)
    {
        mXiC = 5000.0;  // hypothesis two
    }
    else
    {
        mXiC = 0.0;     // hypothesis one
    }
}   

