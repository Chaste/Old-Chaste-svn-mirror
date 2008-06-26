/*

Copyright (C) University of Oxford, 2008

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
#ifndef _FITZHUGHNAGUMO1961ODESYSTEM_HPP_
#define _FITZHUGHNAGUMO1961ODESYSTEM_HPP_

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * Represents the FitzHugh-Nagumo system of ODEs.
 */
class FitzHughNagumo1961OdeSystem : public AbstractCardiacCell
{
private:
    /** 
     *  Constants for the FitzHugh-Nagumo model
     */
    static const double mAlpha = -0.08; 
    static const double mGamma = 3.00;
    static const double mEpsilon = 0.005;    

public:
    // Constructor
    FitzHughNagumo1961OdeSystem(AbstractIvpOdeSolver *pOdeSolver,
                                AbstractStimulusFunction *pIntracelullarStimulus,
                                AbstractStimulusFunction *pExtracelullarStimulus=NULL);
                                
    // Destructor
    ~FitzHughNagumo1961OdeSystem();
    
    // Compute the RHS of the FitHugh-Nagumo system of ODEs
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY);
    double GetIIonic();
};

#endif //_FITZHUGHNAGUMO1961ODESYSTEM_HPP_
