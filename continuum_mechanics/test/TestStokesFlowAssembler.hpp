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

#ifndef TESTSTOKESFLOWASSEMBLER_HPP_
#define TESTSTOKESFLOWASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "StokesFlowAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "QuadraticMesh.hpp"
#include "TrianglesMeshReader.hpp"

class TestStokesFlowAssembler : public CxxTest::TestSuite
{
public:
    /*
     * Test that the matrix is calculated correctly on the cannonical triangle.
     * Tests against the analytical solution calculated by hand.
     */
    void TestAssembler()  throw(Exception)
    {
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/canonical_triangle_quadratic", 2, 2, false);
        mesh.ConstructFromMeshReader(mesh_reader);

        double mu = 2.0;
        c_vector<double,2> body_force = zero_vector<double>(2);

        StokesFlowProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetViscosity(mu);
        problem_defn.SetBodyForce(body_force);

        StokesFlowAssembler<2> assembler(&mesh, &problem_defn);

        Vec vec = PetscTools::CreateVec(15);
        Mat mat;
        PetscTools::SetupMat(mat, 15, 15, 15);

        assembler.SetVectorToAssemble(vec, true);
        assembler.SetMatrixToAssemble(mat, true);
        assembler.Assemble();
        PetscMatTools::Finalise(mat);

        double A[6][6] = {
                           {      1.0,  1.0/6.0,  1.0/6.0,      0.0, -2.0/3.0, -2.0/3.0},
                           {  1.0/6.0,  1.0/2.0,      0.0,      0.0,      0.0, -2.0/3.0},
                           {  1.0/6.0,      0.0,  1.0/2.0,      0.0, -2.0/3.0,      0.0},
                           {      0.0,      0.0,      0.0,  8.0/3.0, -4.0/3.0, -4.0/3.0},
                           { -2.0/3.0,      0.0, -2.0/3.0, -4.0/3.0,  8.0/3.0,      0.0},
                           { -2.0/3.0, -2.0/3.0,      0.0, -4.0/3.0,      0.0,  8.0/3.0}
                         };

        double Bx[6][3] = {
                            { -1.0/6.0,      0.0,      0.0},
                            {      0.0,  1.0/6.0,      0.0},
                            {      0.0,      0.0,      0.0},
                            {  1.0/6.0,  1.0/6.0,  1.0/3.0},
                            { -1.0/6.0, -1.0/6.0, -1.0/3.0},
                            {  1.0/6.0, -1.0/6.0,      0.0},
                         };

        double By[6][3] = {
                            { -1.0/6.0,      0.0,      0.0},
                            {      0.0,      0.0,      0.0},
                            {      0.0,      0.0,  1.0/6.0},
                            {  1.0/6.0,  1.0/3.0,  1.0/6.0},
                            {  1.0/6.0,      0.0, -1.0/6.0},
                            { -1.0/6.0, -1.0/3.0, -1.0/6.0},
                          };

        c_matrix<double,15,15> exact_ael = zero_matrix<double>(15);

        // The diagonal 6x6 blocks
        for (unsigned i=0; i<6; i++)
        {
            for (unsigned j=0; j<6; j++)
            {
                exact_ael(2*i,  2*j)   = mu*A[i][j];
                exact_ael(2*i+1,2*j+1) = mu*A[i][j];
            }
        }

        // The 6x3 Blocks
        for (unsigned i=0; i<6; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                exact_ael(2*i,12+j)   = -Bx[i][j];
                exact_ael(2*i+1,12+j) = -By[i][j];
                //- as -Div(U)=0
                exact_ael(12+j,2*i)   = -Bx[i][j];
                exact_ael(12+j,2*i+1) = -By[i][j];
            }
        }


        int lo, hi;
        MatGetOwnershipRange(mat, &lo, &hi);
        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            for (unsigned j=0; j<15; j++)
            {
                TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,i,j), exact_ael(i,j), 1e-9);
            }
            TS_ASSERT_DELTA(PetscVecTools::GetElement(vec,i), 0.0, 1e-9);
        }

        ReplicatableVector vec_repl(vec);
        for (unsigned i=0; i<15; i++)
        {
            TS_ASSERT_DELTA(vec_repl[i], 0.0, 1e-9);
        }

        VecDestroy(vec);
        MatDestroy(mat);
    }
};

#endif // TESTSTOKESFLOWASSEMBLER_HPP_
