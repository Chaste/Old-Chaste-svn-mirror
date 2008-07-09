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
#ifndef TESTSOLVELAPLACIANWITHQUADRATICS_HPP_
#define TESTSOLVELAPLACIANWITHQUADRATICS_HPP_


#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.hpp"
#include <petsc.h>
#include <vector>
#include <cmath>
#include "BoundaryConditionsContainer.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "QuadraticBasisFunction.hpp"

class TestSolveLaplacianWithQuadratics : public CxxTest::TestSuite
{
public:
    void testSolveLaplacianWithQuadratics() throw (Exception)
    {
        // load mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // figure out number of unknowns
        unsigned num_vertices = mesh.GetNumNodes();
        unsigned num_internal_nodes = 0;
        
        for (ConformingTetrahedralMesh<2,2>::EdgeIterator edge_iterator=mesh.EdgesBegin();
             edge_iterator!=mesh.EdgesEnd();
             ++edge_iterator)
        { 
            num_internal_nodes++;
        }
        
        unsigned num_dofs = num_vertices + num_internal_nodes;
            
        TS_ASSERT_EQUALS(num_dofs, 289u);
                
        // initialise solution vector and linear system
        LinearSystem linear_system(num_dofs);
        linear_system.SetMatrixIsConstant(true);
        
        // compute and store the extra nodes somewhere??
        std::vector<std::vector<double> > lnods(mesh.GetNumElements());
        for(unsigned i=0; i<lnods.size(); i++)
        {
            lnods[i].resize(6);
        }
        


        // set up map from ith component in solution vector to node
                
        // loop over elements and call AssembleOnElement, then add contributions to linear system
        // (AssembleOnElement looks like it'll be the same as with linears, except the sizes of the
        // vector of basis functions and gradients, and therefore Aelem and belem, will be different.
        // + Number of gauss points should change. 
        
        // apply boundary conditions
        
        // solve linear system
        
        // collect solution values at vertices (ie nodes of original mesh)
    }
};
#endif /*TESTSOLVELAPLACIANWITHQUADRATICS_HPP_*/
