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

#ifndef TESTMONODOMAINPURKINJEASSEMBLERS_HPP_
#define TESTMONODOMAINPURKINJEASSEMBLERS_HPP_

#include <cxxtest/TestSuite.h>

#include "MonodomainPurkinjeVolumeAssembler.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscTools.hpp"
#include "LuoRudy1991.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "PetscMatTools.hpp"

class TestMonodomainPurkinjeAssemblers : public CxxTest::TestSuite
{
public:
    void TestVolumeIntegralPartOfMatrix() throw (Exception)
    {
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRegularSlabMesh(0.05, 0.1, 0.1);

    	PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory;
    	cell_factory.SetMesh(&mesh);
    	MonodomainTissue<2> tissue(&cell_factory);


    	Mat purkinje_mat;
    	PetscTools::SetupMat(purkinje_mat, 2*mesh.GetNumNodes(), 2*mesh.GetNumNodes(), 9);
    	Mat normal_mat;
    	PetscTools::SetupMat(normal_mat, mesh.GetNumNodes(), mesh.GetNumNodes(), 9);

    	MonodomainPurkinjeVolumeAssembler<2,2> purkinje_vol_assembler(&mesh, &tissue);
    	purkinje_vol_assembler.SetMatrixToAssemble(purkinje_mat, true);
    	purkinje_vol_assembler.Assemble();

    	MonodomainAssembler<2,2> normal_vol_assembler(&mesh, &tissue);
    	normal_vol_assembler.SetMatrixToAssemble(normal_mat, true);
    	normal_vol_assembler.Assemble();

    	PetscMatTools::Finalise(purkinje_mat);
    	PetscMatTools::Finalise(normal_mat);

    	PetscInt lo, hi;
    	PetscMatTools::GetOwnershipRange(purkinje_mat, lo, hi);

    	for(unsigned i=0; i<mesh.GetNumNodes(); i++)
    	{
    		if(lo<=(int)(2*i) && (int)(2*i)<hi)
    		{
				for(unsigned j=0; j<mesh.GetNumNodes(); j++)
				{
					TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i,2*j),   PetscMatTools::GetElement(normal_mat,i,j), 1e-8);
					TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i,2*j+1), 0.0, 1e-8);
				}
    		}

    		if(lo<=(int)(2*i+1) && (int)(2*i+1)<hi)
    		{
				for(unsigned j=0; j<mesh.GetNumNodes(); j++)
				{
					TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i+1,2*j),   0.0, 1e-8);
					TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i+1,2*j+1), 0.0, 1e-8);
				}
    		}
    	}

    	//PetscMatTools::Display(purkinje_mat);
    	//PetscMatTools::Display(normal_mat);

    	MatDestroy(purkinje_mat);
    	MatDestroy(normal_mat);

    }
};

#endif // TESTMONODOMAINPURKINJEASSEMBLERS_HPP_
