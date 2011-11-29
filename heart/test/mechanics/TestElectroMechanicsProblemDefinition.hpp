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

#ifndef TESTELECTROMECHANICSPROBLEMDEFINITION_HPP_
#define TESTELECTROMECHANICSPROBLEMDEFINITION_HPP_

#include <cxxtest/TestSuite.h>
#include "ElectroMechanicsProblemDefinition.hpp"

class TestElectroMechanicsProblemDefinition : public CxxTest::TestSuite
{
public:
    void TestInterface() throw (Exception)
    {
        QuadraticMesh<2> mesh(0.1, 1.0, 1.0);

        ElectroMechanicsProblemDefinition<2> problem_defn(mesh);

        TS_ASSERT_THROWS_THIS(problem_defn.GetContractionModel(), "Contraction model hasn't been set yet");
        TS_ASSERT_THROWS_THIS(problem_defn.GetContractionModelOdeTimestep(), "Contraction model hasn't been set yet");
        TS_ASSERT_THROWS_THIS(problem_defn.GetMechanicsSolveTimestep(), "Timestep for mechanics solve hasn't been set yet");

        problem_defn.SetContractionModel(NASH2004, 0.01);
        TS_ASSERT_EQUALS(problem_defn.GetContractionModel(), NASH2004);
        TS_ASSERT_EQUALS(problem_defn.GetContractionModelOdeTimestep(), 0.01);

        problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
        NashHunterPoleZeroLaw<2>* p_law = dynamic_cast<NashHunterPoleZeroLaw<2>*>(problem_defn.GetIncompressibleMaterialLaw(0));
        TS_ASSERT(p_law != NULL);

        problem_defn.SetUseDefaultCardiacMaterialLaw(COMPRESSIBLE);
        CompressibleExponentialLaw<2>* p_law_2 = dynamic_cast<CompressibleExponentialLaw<2>*>(problem_defn.GetCompressibleMaterialLaw(0));
        TS_ASSERT(p_law_2 != NULL);

        problem_defn.SetDeformationAffectsElectrophysiology(true,false);
        TS_ASSERT_EQUALS(problem_defn.GetDeformationAffectsConductivity(), true);
        TS_ASSERT_EQUALS(problem_defn.GetDeformationAffectsCellModels(), false);
        problem_defn.SetDeformationAffectsElectrophysiology(false,true);
        TS_ASSERT_EQUALS(problem_defn.GetDeformationAffectsConductivity(), false);
        TS_ASSERT_EQUALS(problem_defn.GetDeformationAffectsCellModels(), true);

        problem_defn.SetMechanicsSolveTimestep(1.0);
        TS_ASSERT_EQUALS(problem_defn.GetMechanicsSolveTimestep(), 1.0);
    }
};

#endif // TESTELECTROMECHANICSPROBLEMDEFINITION_HPP_
