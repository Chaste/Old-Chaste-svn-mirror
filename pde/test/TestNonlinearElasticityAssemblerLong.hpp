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


#ifndef TESTNONLINEARELASTICITYASSEMBLERLONG_HPP_
#define TESTNONLINEARELASTICITYASSEMBLERLONG_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "NonlinearElasticityAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ExponentialMaterialLaw.hpp"
#include "MooneyRivlinMaterialLaw.hpp"

class TestNonlinearElasticityAssemblerLong : public CxxTest::TestSuite
{
public:
    void TestSolve3d() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools
        
        QuadraticMesh<3> mesh("mesh/test/data/cube_136_elements_quadratic");

        MooneyRivlinMaterialLaw<3> law(0.02, 0.0);
        c_vector<double,3> body_force;
        body_force(0) = 0.06;
        body_force(1) = 0.0;
        body_force(2) = 0.0;
        
        std::vector<unsigned> fixed_nodes;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if( fabs(mesh.GetNode(i)->rGetLocation()[0])<1e-6)
            {
                fixed_nodes.push_back(i);
            }
        }
        
        NonlinearElasticityAssembler<3> assembler(&mesh, 
                                                  &law, 
                                                  body_force,
                                                  1.0,
                                                  "simple_nonlin_elas_3d",
                                                  fixed_nodes);
                                                  
        assembler.Solve();
        
        std::vector<c_vector<double,3> >& r_solution = assembler.rGetDeformedPosition();
        
        double xend = 1.23671;
        double yend = 0.00651; // the same as zend
        
        ////////////////////////////////////////////////////////////
        // compare the solution at the corners with the values 
        // obtained using the dealii finite elasticity assembler
        //
        // Results have been visually checked to see they agree 
        // (they do, pretty well - note a slightly more refined
        // mesh was in the dealii simulation).
        ////////////////////////////////////////////////////////////
        
        // node 0 should still be at (0,0,0)
        assert( fabs(mesh.GetNode(0)->rGetLocation()[0] - 0) < 1e-9 );
        assert( fabs(mesh.GetNode(0)->rGetLocation()[1] - 0) < 1e-9 );
        assert( fabs(mesh.GetNode(0)->rGetLocation()[2] - 0) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[0](0), 0.0, 1e-9 );
        TS_ASSERT_DELTA( r_solution[0](1), 0.0, 1e-9 );
        TS_ASSERT_DELTA( r_solution[0](2), 0.0, 1e-9 );
        
        // node 3 should still be at (0,1,0)
        assert( fabs(mesh.GetNode(3)->rGetLocation()[0] - 0) < 1e-9 );
        assert( fabs(mesh.GetNode(3)->rGetLocation()[1] - 1) < 1e-9 );
        assert( fabs(mesh.GetNode(3)->rGetLocation()[2] - 0) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[3](0), 0.0, 1e-9 );
        TS_ASSERT_DELTA( r_solution[3](1), 1.0, 1e-9 );
        TS_ASSERT_DELTA( r_solution[3](2), 0.0, 1e-9 );

        // DEALII value for X=(1,0,0) node is x=(1.23671,0.00651,0.00651)
        assert( fabs(mesh.GetNode(1)->rGetLocation()[0] - 1) < 1e-9 );
        assert( fabs(mesh.GetNode(1)->rGetLocation()[1] - 0) < 1e-9 );
        assert( fabs(mesh.GetNode(1)->rGetLocation()[2] - 0) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[1](0), xend, 1e-2 );
        TS_ASSERT_DELTA( r_solution[1](1), yend, 1e-2 );
        TS_ASSERT_DELTA( r_solution[1](2), yend, 1e-2 );

        // DEALII value for X=(1,1,0) node is x=(1.23671,0.00651,0.00651)
        assert( fabs(mesh.GetNode(2)->rGetLocation()[0] - 1) < 1e-9 );
        assert( fabs(mesh.GetNode(2)->rGetLocation()[1] - 1) < 1e-9 );
        assert( fabs(mesh.GetNode(2)->rGetLocation()[2] - 0) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[2](0),   xend, 1e-2 );
        TS_ASSERT_DELTA( r_solution[2](1), 1-yend, 1e-2 );
        TS_ASSERT_DELTA( r_solution[2](2),   yend, 1e-2 );

        // DEALII value for X=(1,0,1) node is x=(1.23671,0.00651,0.00651)
        assert( fabs(mesh.GetNode(5)->rGetLocation()[0] - 1) < 1e-9 );
        assert( fabs(mesh.GetNode(5)->rGetLocation()[1] - 0) < 1e-9 );
        assert( fabs(mesh.GetNode(5)->rGetLocation()[2] - 1) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[5](0),   xend, 1e-2 );
        TS_ASSERT_DELTA( r_solution[5](1),   yend, 1e-2 );
        TS_ASSERT_DELTA( r_solution[5](2), 1-yend, 1e-2 );

        // DEALII value for X=(1,1,1) node is x=(1.23671,0.00651,0.00651)
        assert( fabs(mesh.GetNode(6)->rGetLocation()[0] - 1) < 1e-9 );
        assert( fabs(mesh.GetNode(6)->rGetLocation()[1] - 1) < 1e-9 );
        assert( fabs(mesh.GetNode(6)->rGetLocation()[2] - 1) < 1e-9 );
        TS_ASSERT_DELTA( r_solution[6](0),   xend, 1e-2 );
        TS_ASSERT_DELTA( r_solution[6](1), 1-yend, 1e-2 );
        TS_ASSERT_DELTA( r_solution[6](2), 1-yend, 1e-2 );
    }
};


#endif /*TESTNONLINEARELASTICITYASSEMBLERLONG_HPP_*/
