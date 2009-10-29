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


#ifndef TESTEXPLICITCARDIACMECHANICSASSEMBLER_HPP_
#define TESTEXPLICITCARDIACMECHANICSASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "ExplicitCardiacMechanicsAssembler.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "QuadraturePointsGroup.hpp"
#include "NonlinearElasticityTools.hpp"
#include "ReplicatableVector.hpp"

class TestExplicitCardiacMechanicsAssembler : public CxxTest::TestSuite
{
public:
    void TestWithSimpleContractionModel() throw(Exception)
    {
        EXIT_IF_PARALLEL; // unlike above tests, this one doesn't pass in parallel (only written for sequential)

        QuadraticMesh<2> mesh(1.0, 1.0, 4, 4);
        MooneyRivlinMaterialLaw<2> law(1);


        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0.0);

        // coverage (bad contraction model)
        TS_ASSERT_THROWS_THIS(ExplicitCardiacMechanicsAssembler<2> assembler(NHS,&mesh,"",fixed_nodes,&law), "Unknown or stretch-rate-dependent contraction model");
          
        // TEST1 => NonphysiologicalContractionModel 1
        ExplicitCardiacMechanicsAssembler<2> assembler(TEST1,&mesh,"TestExplicitCardiacMech",fixed_nodes,&law);
        QuadraturePointsGroup<2> quad_points(mesh, *(assembler.GetQuadratureRule()));

        std::vector<double> calcium_conc(assembler.GetTotalNumQuadPoints(), 0.0);
        std::vector<double> voltages(assembler.GetTotalNumQuadPoints(), 0.0);

        assembler.SetCalciumAndVoltage(calcium_conc, voltages);

        // solve UP TO t=0. So Ta(lam_n,t_{n+1})=5*sin(0)=0, ie no deformation
        assembler.Solve(-0.01,0.0,0.01);        
        TS_ASSERT_EQUALS(assembler.GetNumNewtonIterations(),0u);

        assembler.Solve(0.99,1.0,0.01);
        
        TS_ASSERT_DELTA(assembler.rGetDeformedPosition()[4](0),  0.8730, 1e-2);
        TS_ASSERT_DELTA(assembler.rGetDeformedPosition()[4](1), -0.0867, 1e-2);
    }
    
    //EMTODO1: compare with implicit when stretch independent (also when stretch dependent?)
};

#endif /*TESTEXPLICITCARDIACMECHANICSASSEMBLER_HPP_*/
