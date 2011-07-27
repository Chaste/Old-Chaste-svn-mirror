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

#ifndef TESTNATURALNEUMANNSURFACETERMASSEMBLER_HPP_
#define TESTNATURALNEUMANNSURFACETERMASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "NaturalNeumannSurfaceTermAssembler.hpp"
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestNaturalNeumannSurfaceTermAssembler : public CxxTest::TestSuite
{
public:
    void TestAssembler() throw (Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_mesh_5_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_neumann_boundary_condition = new ConstBoundaryCondition<1>(1.0);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorEnd();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, p_neumann_boundary_condition);

        NaturalNeumannSurfaceTermAssembler<1,1,1> assembler(&mesh, &bcc);

        Vec vec = PetscTools::CreateVec(mesh.GetNumNodes());
        assembler.SetVectorToAssemble(vec,true);

        assembler.Assemble();

    }
};

#endif // TESTNATURALNEUMANNSURFACETERMASSEMBLER_HPP_
