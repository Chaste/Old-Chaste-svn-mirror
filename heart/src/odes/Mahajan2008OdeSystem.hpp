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

#ifndef _MAHAJAN2008ODESYSTEM_
#define _MAHAJAN2008ODESYSTEM_

#include <cmath>
#include <cassert>
#include "AbstractCardiacCell.hpp"
#include "Exception.hpp"
#include "OdeSystemInformation.hpp"
#include "AbstractStimulusFunction.hpp"

class Mahajan2008OdeSystem : public AbstractCardiacCell
{
private:
    /** Scale factor for Gks*/
    double mScaleFactorGks;
    /** Scale factor for Gto*/
    double mScaleFactorIto;
    /** Scale factor for Gkr*/
    double mScaleFactorGkr;
public:
    Mahajan2008OdeSystem(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                         boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    ~Mahajan2008OdeSystem();

    /**
     * verify the state variables
     */
    void VerifyGatingVariables();

     /**
     * Set the scale factor for Gks in order to differentiate epi M and endo cells
     */
    void SetScaleFactorGks(double sfgks);

    /**
     * Set the scale factor for Gks in order to differentiate epi M and endo cells
     */
    void SetScaleFactorIto(double sfito);

    /**
     * Set the scale factor for Gkr (might be useful for simulating LQT2)
     */
    void SetScaleFactorGkr(double sfgkr);

    /**
     * Calculates the ionic current
     */
    double GetIIonic();

    /**
     * Evaluate the derivatives of the state variables
     */
    void EvaluateYDerivatives(double var_Environment__time, const std::vector<double> &rY, std::vector<double> &rDY);

};


#endif
