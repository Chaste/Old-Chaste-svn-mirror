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


#ifndef GROWTHBYCONSTANTMASSODESYSTEM_HPP_
#define GROWTHBYCONSTANTMASSODESYSTEM_HPP_

#include <vector>
#include <cmath>
#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "AbstractGrowingTumourSourceModel.hpp"


/**
 */
template<unsigned DIM>
class GrowthByConstantMassOdeSystem : public AbstractOdeSystem
{
private:
    double mRho;
    unsigned mSourceModelIndex;

    AbstractGrowingTumourSourceModel<DIM>* mpSourceModel;

public:
    // Constructor
    GrowthByConstantMassOdeSystem(double rho,
                                  unsigned sourceModelIndex,
                                  AbstractGrowingTumourSourceModel<DIM>* pSourceModel)
            : AbstractOdeSystem(1)
    {
        mSourceModelIndex = sourceModelIndex;
        assert(rho>0);
        mRho = rho;
        assert(pSourceModel!=NULL);
        mpSourceModel = pSourceModel;
        
        mpSystemInfo = OdeSystemInformation<GrowthByConstantMassOdeSystem<DIM> >::Instance(); 
         
        SetStateVariables(GetInitialConditions());
      }

    ~GrowthByConstantMassOdeSystem()
    {}


    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
    {
        double s = mpSourceModel->GetSourceValue(mSourceModelIndex);
        rDY[0] = (1.0/DIM)*rY[0]*mRho*s;
    }
    

};

///\todo Don't know how to make this template generic
template<>
void OdeSystemInformation< GrowthByConstantMassOdeSystem<2> >::Initialise(void)
{
    this->mVariableNames.push_back("Dunno");
    this->mVariableUnits.push_back("Dunno");
    this->mInitialConditions.push_back(1.0);
    this->mInitialised = true;
}
template<>
void OdeSystemInformation< GrowthByConstantMassOdeSystem<3> >::Initialise(void)
{
    this->mVariableNames.push_back("Dunno");
    this->mVariableUnits.push_back("Dunno");
    this->mInitialConditions.push_back(1.0);
    this->mInitialised = true;
}

#endif /*GROWTHBYCONSTANTMASSODESYSTEM_HPP_*/
