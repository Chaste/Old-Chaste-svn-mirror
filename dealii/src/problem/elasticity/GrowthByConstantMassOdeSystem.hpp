/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GROWTHBYCONSTANTMASSODESYSTEM_HPP_
#define GROWTHBYCONSTANTMASSODESYSTEM_HPP_

#include <vector>
#include <cmath>
#include "AbstractOdeSystem.hpp"
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
        
        mInitialConditions.push_back(1.0);
        
        SetStateVariables(mInitialConditions);
    }
    
    ~GrowthByConstantMassOdeSystem()
    {}
    
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
    {
        double s = mpSourceModel->GetSourceValue(mSourceModelIndex);
        rDY[0] = (1.0/DIM)*rY[0]*mRho*s;
    }
    
};
#endif /*GROWTHBYCONSTANTMASSODESYSTEM_HPP_*/
