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

#ifdef CHASTE_CVODE

#include <sstream>
#include <cassert>

#include "AbstractCvodeSystem.hpp"
#include "Exception.hpp"
#include "VectorHelperFunctions.hpp"

AbstractCvodeSystem::AbstractCvodeSystem(unsigned numberOfStateVariables)
    : AbstractParameterisedSystem<N_Vector>(numberOfStateVariables),
      mUseAnalyticJacobian(false)
{
}

void AbstractCvodeSystem::Init()
{
    DeleteVector(mStateVariables);
    mStateVariables = GetInitialConditions();
    DeleteVector(mParameters);
    mParameters = N_VNew_Serial(rGetParameterNames().size());
    for (int i=0; i<NV_LENGTH_S(mParameters); i++)
    {
        NV_Ith_S(mParameters, i) = 0.0;
    }
}

AbstractCvodeSystem::~AbstractCvodeSystem()
{
}

//
//double AbstractCvodeSystem::CalculateRootFunction(double time, const std::vector<double>& rY)
//{
//    bool stop = CalculateStoppingEvent(time, rY);
//    return stop ? 0.0 : 1.0;
//}
//
//bool AbstractCvodeSystem::GetUseAnalyticJacobian()
//{
//    return mUseAnalyticJacobian;
//}

//
//std::string AbstractCvodeSystem::DumpState(const std::string& message,
//                                           N_Vector Y)
//{
//    std::stringstream res;
//    res << message << "\nState:\n";
//    if (Y == NULL)
//    {
//        Y = mStateVariables;
//    }
//    for (int i=0; i<NV_LENGTH_S(Y); i++)
//    {
//        res << "\t" << rGetStateVariableNames()[i] << ":" << NV_Ith_S(Y, i) << "\n";
//    }
//    return res.str();
//}

#endif // CHASTE_CVODE
