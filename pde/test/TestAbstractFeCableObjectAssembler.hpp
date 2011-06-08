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
#ifndef TESTABSTRACTFECABLEOBJECTASSEMBLER_HPP_
#define TESTABSTRACTFECABLEOBJECTASSEMBLER_HPP_


#include <cxxtest/TestSuite.h>

#include "AbstractFeCableObjectAssembler.hpp"
#include "MixedDimensionMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscVecTools.hpp"
#include "ReplicatableVector.hpp"
#include "TrianglesMeshReader.hpp"

template<unsigned DIM>
class BasicCableVectorAssembler : public AbstractFeCableObjectAssembler<DIM,DIM,1,true,false,NORMAL>
{
private:
    double mCoefficient;
  
    c_vector<double,1*2> ComputeCableVectorTerm(
        c_vector<double, 2>& rPhi,
        c_matrix<double, DIM, 2>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double, 1, DIM>& rGradU,
        Element<1,DIM>* pElement)
    {
        return -mCoefficient*rPhi;
    }

public:
    BasicCableVectorAssembler(MixedDimensionMesh<DIM,DIM>* pMesh, double coefficient)
        : AbstractFeCableObjectAssembler<DIM,DIM,1,true,false,NORMAL>(pMesh),
          mCoefficient(coefficient)
    {
    }
};


class TestAbstractFeCableObjectAssembler : public CxxTest::TestSuite
{
public:
    void TestBasicVectorAssemblers() throw(Exception)
    {
        std::string mesh_base("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
        TrianglesMeshReader<2,2> reader(mesh_base);
        MixedDimensionMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader); 
              
        Vec vec = PetscTools::CreateVec(mesh.GetNumNodes());

        BasicCableVectorAssembler<2> assembler(&mesh, 2.0);

        assembler.SetVectorToAssemble(vec,true);
        assembler.Assemble();

        PetscVecTools::Assemble(vec);

        /*
         * Cables:
        0       55      56      1
        1       56      57      2
        2       57      58      3
        3       58      59      4
        4       59      60      5
        5       60      61      6
        6       61      62      7
        7       62      63      8
        8       63      64      9
        9       64      65      10
         *
         */

        ReplicatableVector vec_repl(vec);
        for (unsigned i = 0; i < 55; ++i)
        {
            TS_ASSERT_DELTA(vec_repl[i], 0.0, 1e-4); 
        }
        TS_ASSERT_DELTA(vec_repl[55], -0.01, 1e-4);
        for (unsigned i=56; i<65; i++)
        {
            TS_ASSERT_DELTA(vec_repl[i], -0.02, 1e-4);
        }
        TS_ASSERT_DELTA(vec_repl[65], -0.01, 1e-4);
        for (unsigned i = 66; i < mesh.GetNumNodes(); ++i)
        {
            TS_ASSERT_DELTA(vec_repl[i], 0.0, 1e-4); 
        }
        
        VecDestroy(vec);
    }
};




#endif /*TESTABSTRACTFECABLEOBJECTASSEMBLER_HPP_*/
