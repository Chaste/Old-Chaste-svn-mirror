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


#ifndef SIMPLETUMOURSOURCEMODEL_HPP_
#define SIMPLETUMOURSOURCEMODEL_HPP_

#include "AbstractGrowingTumourSourceModel.hpp"

/**
 * A simple tumour source model (for testing) that just returns
 * s = n, at an evaluation point, where n is the index of the evaluation
 * point
 */
template<unsigned DIM>
class SimpleTumourSourceModel : public AbstractGrowingTumourSourceModel<DIM>
{
public :
    SimpleTumourSourceModel()
       : AbstractGrowingTumourSourceModel<DIM>()
    {
    }

    void Run(double tStart, double tEnd, FiniteElasticityAssembler<DIM>* pFiniteElasticity)
    {
        typename std::map<unsigned,EvaluationPointInfo<DIM> >::iterator iter
        = this->mEvaluationPoints.begin();
        while (iter!=this->mEvaluationPoints.end())
        {
            unsigned index = iter->first;
            iter->second.SourceValue = index;
            iter++;
        }
    }
};

#endif /*SIMPLETUMOURSOURCEMODEL_HPP_*/
