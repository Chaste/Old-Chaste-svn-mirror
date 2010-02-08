/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#include "IngeWntSwatCellCycleOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"

IngeWntSwatCellCycleOdeSystem::IngeWntSwatCellCycleOdeSystem(unsigned hypothesis, double wntLevel, const CryptCellMutationState& rMutationState)
    : AbstractOdeSystem(22),
      mMutationState(rMutationState),
      mHypothesis(hypothesis)
{
    if (hypothesis!=1u && hypothesis!=2u)
    {
        EXCEPTION("You must set up this cell cycle ODE system with hypothesis one or two.");
    }

    mpSystemInfo.reset(new CellwiseOdeSystemInformation<IngeWntSwatCellCycleOdeSystem>);

    /**
     * State variables are
     *
     *  0. r = pRb
     *  1. e = E2F1 (This is the S-phase indicator)
     *  2. i = CycD (inactive)
     *  3. j = CycD (active)
     *  4. p = pRb-p
     *  5. D = APC destruction complex
     *  6. X = Axin
     *  7. Cu = Beta Cat marked for ubiquitination
     *  8. Co = Open form Beta Cat
     *  9. Cc = Closed form Beta Cat
     *  10. Mo = Open form Mutant Beta Cat
     *  11. Mc = Closed form Mutant Beta Cat
     *  12. A = Free Adhesion molecules
     *  13. Ca = BetaCat/Adhesion
     *  14. Ma = Mutant BetaCat/Adhesion
     *  15. T = free TCF
     *  16. Cot = Open BetaCat/TCF
     *  17. Cct = Closed BetaCat/TCF
     *  18. Mot = Open Mutant BetaCat/TCF
     *  19. Mct = Closed Mutant BetaCat/TCF
     *  20. Y = Wnt Target protein
     *  21. Wnt level
     */
    Init(); // set up parameter values

    double d_d_hat = mDd + mXiD*wntLevel;
    double d_d_x_hat = mDdx + mXiDx*wntLevel;
    double d_x_hat = mDx + mXiX*wntLevel;
    double p_c_hat = mPc + mXiC*wntLevel;

    double sigma_D = 0.0; // for healthy cells
    double sigma_B = 0.0; // for healthy cells

    switch (mMutationState)
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
            // This can't happen if all mutation states are catered for
            NEVER_REACHED;
    }

    // Cell-specific initial conditions
    double steady_D = ((1.0-sigma_D)*mSd*mSx)/((1.0-sigma_D)*mSd*d_d_hat + d_x_hat*(d_d_hat + d_d_x_hat));
    SetInitialConditionsComponent(5, steady_D); // Destruction complex (APC/Axin/GSK3B)

    double temp = (mSx*(d_d_hat+d_d_x_hat))/((1.0-sigma_D)*mSd*d_d_hat+d_x_hat*(d_d_hat+d_d_x_hat));
    SetInitialConditionsComponent(6, temp);  // Axin

    double steady_Cf = ((mSc-mDc*mKd - mPu*steady_D)+sqrt(SmallPow((mSc-mDc*mKd - mPu*steady_D),2) + (4.0*mSc*mDc*mKd)))/(2.0*mDc);
    temp = (mPu*steady_D*steady_Cf)/(mDu*(steady_Cf+mKd));
    SetInitialConditionsComponent(7, temp); // beta-catenin to be ubiquitinated

    double theta = mDc + (mPu*steady_D)/(steady_Cf + mKd);

    double steady_Co = ( mSc - p_c_hat - theta*mKc + sqrt(4.0*mSc*theta*mKc + SmallPow((mSc - p_c_hat - theta*mKc),2)) )/(2.0*theta);
    SetInitialConditionsComponent(8, steady_Co); // Open form beta-catenin

    double steady_Cc = steady_Cf - steady_Co;

    if ((steady_Cc < 0) && (steady_Cc+100*DBL_EPSILON > 0) ) // stop protein values going -ve
    {
        steady_Cc = 0.0;
    }
    SetInitialConditionsComponent(9, steady_Cc); // Closed form beta-catenin

    SetInitialConditionsComponent(12, mSa/mDa);  // 'Free' adhesion molecules

    SetInitialConditionsComponent(13, mSa*mSca*steady_Co/(mDa*mDca)); // Co-A Adhesion complex

    SetInitialConditionsComponent(15, mSt/mDt); // `Free' transcription molecules (TCF)

    SetInitialConditionsComponent(16, mSct*mSt*steady_Co/(mDt*mDct)); // Co-T open form beta-catenin/TCF

    SetInitialConditionsComponent(17, mSct*mSt*steady_Cc/(mDt*mDct)); // Cc-T closed beta-catenin/TCF

    temp = (mSct*mSt*mSy*steady_Cf)/(mDy*(mSct*mSt*steady_Cf + mDct*mDt*mKt));
    SetInitialConditionsComponent(20, temp); // Wnt target protein

    SetInitialConditionsComponent(21, wntLevel); // Wnt stimulus
}

void IngeWntSwatCellCycleOdeSystem::SetMutationState(const CryptCellMutationState& rMutationState)
{
    mMutationState = rMutationState;
}

IngeWntSwatCellCycleOdeSystem::~IngeWntSwatCellCycleOdeSystem()
{
    // Do nothing
}

void IngeWntSwatCellCycleOdeSystem::Init()
{
    // Swat (2004) parameters
    double k1 = 1.0;
    double k2 = 1.6;
    double k3 = 0.05;
    double k16 = 0.4;
    double k34 = 0.04;
    double k43 = 0.01;
    double k61 = 0.3;
    double k23 = 0.3;
    double a = 0.04;
    double J11 = 0.5;
    double J12 = 5.0;
    double J61 = 5.0;
    double J62 = 8.0;
    double J13 = 0.002;
    double J63 = 2.0;
    double Km1 = 0.5;
    double Km2 = 4.0;
    double Km4 = 0.3;
    double kp = 0.05;
    double phi_pRb = 0.005;
    double phi_E2F1 = 0.1;
    double phi_CycDi = 0.023;
    double phi_CycDa = 0.03;
    double phi_pRbp = 0.06;

    // Value of the mitogenic factor to make the van Leeuwen model influence cell cycle just the same
    double mitogenic_factorF = 1.0/25.0;

    // Non-dimensionalise parameters
    mk2d = k2/(Km2*phi_E2F1);
    mk3d = k3*mitogenic_factorF/(Km4*phi_E2F1);
    mk34d = k34/phi_E2F1;
    mk43d = k43/phi_E2F1;
    mk23d = k23*Km2/(Km4*phi_E2F1);
    mad = a/Km2;
    mJ11d = J11*phi_E2F1/k1;
    mJ12d = J12*phi_E2F1/k1;
    mJ13d = J13*phi_E2F1/k1;
    mJ61d = J61*phi_E2F1/k1;
    mJ62d = J62*phi_E2F1/k1;
    mJ63d = J63*phi_E2F1/k1;
    mKm1d = Km1/Km2;
    mkpd = kp/(Km2*phi_E2F1);
    mphi_r = phi_pRb/phi_E2F1;
    mphi_i = phi_CycDi/phi_E2F1;
    mphi_j = phi_CycDa/phi_E2F1;
    mphi_p = phi_pRbp/phi_E2F1;
    mk16d = k16*Km4/phi_E2F1;
    mk61d = k61/phi_E2F1;
    mPhiE2F1 = phi_E2F1;

    // Initialize van Leeuwen model parameters
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
    if (mHypothesis==1)
    {
        mXiC = 0.0; //  h^-1 (FOR HYPOTHESIS ONE)
    }
    else
    {
        mXiC = 5000.0;   //  h^-1 (FOR HYPOTHESIS TWO)
    }
}

void IngeWntSwatCellCycleOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    double r = rY[0];
    double e = rY[1];
    double i = rY[2];
    double j = rY[3];
    double p = rY[4];

    double dx1 = 0.0;
    double dx2 = 0.0;
    double dx3 = 0.0;
    double dx4 = 0.0;
    double dx5 = 0.0;

    // Bit back-to-front, but work out the Wnt section first...

    // Variables
    double D = rY[5];
    double X = rY[6];
    double Cu = rY[7];
    double Co = rY[8];
    double Cc = rY[9];
    double Mo = rY[10];
    double Mc = rY[11];
    double A = rY[12];
    double Ca = rY[13];
    double Ma = rY[14];
    double T = rY[15];
    double Cot = rY[16];
    double Cct = rY[17];
    double Mot = rY[18];
    double Mct = rY[19];
    double Y = rY[20];
    double stimulus_wnt = rY[21];

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

    switch (mMutationState)
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
            // This can't happen if all mutation states are catered for
            NEVER_REACHED;
    }

    // Now the cell cycle stuff...

    // dr
    dx1 = e/(mKm1d+e)*mJ11d/(mJ11d+r)*mJ61d/(mJ61d+p) - mk16d*r*j+mk61d*p-mphi_r*r;
    // de
    dx2 = mkpd+mk2d*(mad*mad+e*e)/(1+e*e)*mJ12d/(mJ12d+r)*mJ62d/(mJ62d+p) - e;
    // di - changed to include Ct+Mt - transcriptional beta-catenin
    dx3 = mk3d*(Ct+Mt) + mk23d*e*mJ13d/(mJ13d+r)*mJ63d/(mJ63d+p) + mk43d*j - mk34d*i*j/(1+j) - mphi_i*i;
    // dj
    dx4 = mk34d*i*j/(1+j) - (mk43d+mphi_j)*j;
    // dp
    dx5 = mk16d*r*j - mk61d*p - mphi_p*p;

    double factor = mPhiE2F1*60.0;  // Convert non-dimensional d/dt s to d/dt in hours.

    rDY[0] = dx1*factor;
    rDY[1] = dx2*factor;
    rDY[2] = dx3*factor;
    rDY[3] = dx4*factor;
    rDY[4] = dx5*factor;

    // The van Leeuwen ODE system
    rDY[5] = (1.0-sigma_D)*mSd*X - (d_d_hat + d_d_x_hat)*D;
    rDY[6] = mSx - (1.0-sigma_D)*mSd*X - d_x_hat*X + d_d_x_hat*D;
    rDY[7] = (mPu*D*Cf)/(Cf+mKd) - mDu*Cu;

    rDY[8] = (1.0-sigma_B)*mSc + mDca*Ca + mDct*Cot - (mSca*A + mSct*T + mDc)*Co
             - (p_c_hat*Co)/(Co + Mo + mKc) - (mPu*D*Co)/(Cf+mKd);

    rDY[9] = (p_c_hat*Co)/(Co + Mo + mKc) + mDct*Cct - (mSct*T + mDc)*Cc
             - (mPu*D*Cc)/(Cf+mKd);

    rDY[10] = sigma_B*mSc + mDca*Ma + mDct*Mot - (mSca*A + mSct*T + mDc)*Mo
             - (p_c_hat*Mo)/(Co + Mo + mKc);

    rDY[11] = (p_c_hat*Mo)/(Co + Mo + mKc) + mDct*Mct - (mSct*T + mDc)*Mc;
    rDY[12] = mSa + mDca*(Ca+Ma) - (mSca*(Co+Mo) + mDa)*A;
    rDY[13] = mSca*Co*A - mDca*Ca;
    rDY[14] = mSca*Mo*A - mDca*Ma;
    rDY[15] = mSt + mDct*(Ct+Mt) - mSct*(Cf+Mf)*T - mDt*T;
    rDY[16] = mSct*Co*T - mDct*Cot;
    rDY[17] = mSct*Cc*T - mDct*Cct;
    rDY[18] = mSct*Mo*T - mDct*Mot;
    rDY[19] = mSct*Mc*T - mDct*Mct;
    rDY[20] = (mSy*(Ct+Mt))/(Ct + Mt + mKt) - mDy*Y;
    rDY[21] = 0.0;  // don't interfere with Wnt stimulus
}

CryptCellMutationState& IngeWntSwatCellCycleOdeSystem::rGetMutationState()
{
    return mMutationState;
}

bool IngeWntSwatCellCycleOdeSystem::CalculateStoppingEvent(double time, const std::vector<double>& rY)
{
    std::vector<double> dy(rY.size());
    EvaluateYDerivatives(time, rY, dy);

    return (rY[1] > 1.0 && dy[1] > 0.0);
}

double IngeWntSwatCellCycleOdeSystem::CalculateRootFunction(double time, const std::vector<double>& rY)
{
    return rY[1]-1.0;
}


template<>
void CellwiseOdeSystemInformation<IngeWntSwatCellCycleOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("pRb");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(7.357000000000000e-01);

    this->mVariableNames.push_back("E2F1");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(1.713000000000000e-01);

    this->mVariableNames.push_back("CycD_i");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(6.900000000000001e-02);

    this->mVariableNames.push_back("CycD_a");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(3.333333333333334e-03);

    this->mVariableNames.push_back("pRb_p");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(1.000000000000000e-04);

    this->mVariableNames.push_back("D");  // Destruction complex (APC/Axin/GSK3B)
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("X");  // Axin
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Cu"); // beta-catenin to be ubiquitinated
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Co"); // Open form beta-catenin
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Cc"); // Closed form beta-catenin
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Mo"); // Open form mutant beta-catenin
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("Mc"); // Closed form mutant beta-catenin
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("A");  // 'Free' adhesion molecules
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Ca"); // Co-A Adhesion complex
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Ma"); // Mo-A Mutant adhesion complex
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("T"); // `Free' transcription molecules (TCF)
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Cot"); // Co-T open form beta-catenin/TCF
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Cct"); // Cc-T closed beta-catenin/TCF
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Mot"); // Mo-T open form mutant beta-catenin/TCF
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("Mct"); // Mc-T closed form mutant beta-catenin/TCF
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("Y"); // Wnt target protein
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Sw"); // Wnt stimulus
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mInitialised = true;
}
