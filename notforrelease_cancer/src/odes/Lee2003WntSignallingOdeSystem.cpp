/*

Copyright (C) University of Oxford, 2005-2009

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
#include "Lee2003WntSignallingOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"

Lee2003WntSignallingOdeSystem::Lee2003WntSignallingOdeSystem(double wntStimulus)
    : AbstractOdeSystem(8)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Lee2003WntSignallingOdeSystem>);

    /**
     * The state variables are
     *
     *  0. X2 Dsh_active
     *  1. X3 APC/axin/GSK3
     *  2. X4 APC/axin/GSK3
     *  3. X9 beta-cat/APC/axin/GSK3
     *  4. X10 beta-cat
     *  5. X11 beta-cat
     *  6. X12 axin
     *  7. WntLevel
     */

    Init(); // set up parameters values

    // Cell-specific initial conditions
    SetInitialConditionsComponent(7u, wntStimulus);
}

Lee2003WntSignallingOdeSystem::~Lee2003WntSignallingOdeSystem()
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

void Lee2003WntSignallingOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
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

    for (unsigned i=0; i<rY.size(); i++)
    {
        assert( rY[i] >= 0.0 ); // all protein concentrations are positive...
    }

    // Easy ones - A.32, A.37, A.34, A.35
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

    double factor = 60.0;  // convert d/dt in minutes to d/dt in hours

    rDY[0] = dX2*factor;
    rDY[1] = dX4*factor;
    rDY[2] = dX9*factor;
    rDY[3] = dX10*factor;
    rDY[4] = dX3*factor;
    rDY[5] = dX11*factor;
    rDY[6] = dX12*factor;
    rDY[7] = 0.0; // do not change the Wnt level
}


template<>
void CellwiseOdeSystemInformation<Lee2003WntSignallingOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Dsh_active");
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("APC_axin_GSK3");
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(4.83e-3);

    this->mVariableNames.push_back("beta_cat_P_APC_P_axin_P_GSK3");
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(2.02e-3);

    this->mVariableNames.push_back("beta_cat_P");
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(1.00);

    this->mVariableNames.push_back("APC_P_axin_P_GSK3");
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(9.66e-3);

    this->mVariableNames.push_back("beta_cat");
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(25.1);

    this->mVariableNames.push_back("axin");
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(4.93e-4);

    this->mVariableNames.push_back("Wnt");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mInitialised = true;
}
