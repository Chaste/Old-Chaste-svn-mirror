/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef STOKESFLOWPROBLEMDEFINITION_HPP_
#define STOKESFLOWPROBLEMDEFINITION_HPP_

#include "ContinuumMechanicsProblemDefinition.hpp"

template<unsigned DIM>
class StokesFlowProblemDefinition : public ContinuumMechanicsProblemDefinition<DIM>
{
private:
    double mMu;

public:
    StokesFlowProblemDefinition(QuadraticMesh<DIM>& rMesh)
        : ContinuumMechanicsProblemDefinition<DIM>(rMesh),
          mMu(-1.0)
    {
    }

    void SetViscosity(double mu)
    {
        assert(mu > 0.0);
        mMu = mu;
    }

    double GetViscosity()
    {
        if(mMu < 0.0)
        {
            EXCEPTION("Viscosity hasn't been set yet (for the Stokes' flow problem)");
        }
        return mMu;
    }

    void SetZeroFlowNodes(std::vector<unsigned>& rZeroFlowNodes)
    {
        this->SetZeroDirichletNodes(rZeroFlowNodes);
    }

    void SetPrescribedFlowNodes(std::vector<unsigned>& rPrescribedFlowNodes, std::vector<c_vector<double,DIM> >& rPrescribedFlow)
    {
        assert(rPrescribedFlowNodes.size()==rPrescribedFlow.size());

        this->mDirichletNodes = rPrescribedFlowNodes;
        this->mDirichletNodeValues = rPrescribedFlow;
    }
};


#endif // STOKESFLOWPROBLEMDEFINITION_HPP_
