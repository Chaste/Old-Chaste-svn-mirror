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


#ifndef TESTCARDIACELECTROMECHANICSPROBLEM2LONG_HPP_
#define TESTCARDIACELECTROMECHANICSPROBLEM2LONG_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechanicsProblem2.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

class TestCardiacElectroMechanicsProblem2Long : public CxxTest::TestSuite
{
public:
    void Test2dHardcodedResult() throw(Exception)
    {
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 2> cell_factory(-1000*1000);

        CardiacElectroMechanicsProblem2<2> problem(&cell_factory,
                                                   10, /* end time */
                                                   5, /*mech mesh size*/
                                                   100, /* 100*0.01ms mech dt */
                                                   0.01, /*NHS ode dt */
                                                   "TestCardiacElectroMechImplicit2");
        problem.SetNoElectricsOutput();
        problem.Solve();

        // test by checking the length of the tissue against hardcoded value
        std::vector<c_vector<double,2> >& r_deformed_position = problem.rGetDeformedPosition();
        TS_ASSERT_DELTA(r_deformed_position[5](0), 0.998313, 1e-4);
    }
};
#endif /*TESTCARDIACELECTROMECHANICSPROBLEM2LONG_HPP_*/
