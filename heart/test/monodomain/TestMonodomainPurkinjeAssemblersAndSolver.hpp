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

#ifndef TESTMONODOMAINPURKINJEASSEMBLERSANDSOLVER_HPP_
#define TESTMONODOMAINPURKINJEASSEMBLERSANDSOLVER_HPP_

#include <cxxtest/TestSuite.h>

#include "MonodomainPurkinjeVolumeAssembler.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "TetrahedralMesh.hpp"
#include "MixedDimensionMesh.hpp"
#include "PetscTools.hpp"
#include "LuoRudy1991.hpp"
#include "DiFrancescoNoble1985.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "MonodomainPurkinjeCableAssembler.hpp"
#include "PetscMatTools.hpp"
#include "MonodomainPurkinjeSolver.hpp"


class PurkinjeCellFactory : public AbstractPurkinjeCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    PurkinjeCellFactory()
        : AbstractPurkinjeCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-6000.0, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        ChastePoint<2> location = GetMesh()->GetNode(node)->GetPoint();

        if (fabs(location[0])<1e-6)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }

    AbstractCardiacCell* CreatePurkinjeCellForTissueNode(unsigned node)
    {
        return new CellDiFrancescoNoble1985FromCellML(mpSolver, mpZeroStimulus);
    }
};



class TestMonodomainPurkinjeAssemblersAndSolvers : public CxxTest::TestSuite
{
public:
    void TestMonodomainPurkinjeVolumeAssembler() throw (Exception)
    {
    	PdeSimulationTime::SetPdeTimeStep(0.01);

    	DistributedTetrahedralMesh<2,2> mesh;
        mesh.ConstructRegularSlabMesh(0.05, 0.1, 0.1);

    	PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory;
    	cell_factory.SetMesh(&mesh);
    	MonodomainTissue<2> tissue(&cell_factory);

    	// Make sure that a 2Nx2N matrix is partitioned in the same place as an NxN matrix.
    	unsigned num_local_nodes = mesh.GetDistributedVectorFactory()->GetLocalOwnership();
    	Mat purkinje_mat;
    	PetscTools::SetupMat(purkinje_mat, 2*mesh.GetNumNodes(), 2*mesh.GetNumNodes(), 9, 2*num_local_nodes, 2*num_local_nodes);
    	Mat normal_mat;
    	PetscTools::SetupMat(normal_mat, mesh.GetNumNodes(), mesh.GetNumNodes(), 9, num_local_nodes, num_local_nodes);

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

    	//Check that the partitioning is exactly as expected
    	TS_ASSERT_EQUALS((unsigned)lo, 2*mesh.GetDistributedVectorFactory()->GetLow());
    	TS_ASSERT_EQUALS((unsigned)hi, 2*mesh.GetDistributedVectorFactory()->GetHigh());

        for (AbstractTetrahedralMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
        		node_iter != mesh.GetNodeIteratorEnd(); ++node_iter)
    	{
    		unsigned i = node_iter->GetIndex();
        	assert(lo<=(int)(2*i) && (int)(2*i)<hi);
			for(unsigned j=0; j<mesh.GetNumNodes(); j++)
			{
				TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i,2*j),   PetscMatTools::GetElement(normal_mat,i,j), 1e-8);
				TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i,2*j+1), 0.0, 1e-8);
			}

    		assert(lo<=(int)(2*i+1) && (int)(2*i+1)<hi);
			for(unsigned j=0; j<mesh.GetNumNodes(); j++)
			{
				TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i+1,2*j),   0.0, 1e-8);
				TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i+1,2*j+1), 0.0, 1e-8);
			}
    	}

    	MatDestroy(purkinje_mat);
    	MatDestroy(normal_mat);

    }

    void TestMonodomainPurkinjeCableAssembler() throw(Exception)
    {
    	PdeSimulationTime::SetPdeTimeStep(0.01);

        std::string mesh_base("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
        TrianglesMeshReader<2,2> reader(mesh_base);

        ///\todo There are named indices in the test, so we need predictable numbering...
        MixedDimensionMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(reader);

        Mat purkinje_mat;
    	unsigned num_local_nodes = mesh.GetDistributedVectorFactory()->GetLocalOwnership();
        PetscTools::SetupMat(purkinje_mat, 2*mesh.GetNumNodes(), 2*mesh.GetNumNodes(), 9, 2*num_local_nodes, 2*num_local_nodes);

		MonodomainPurkinjeCableAssembler<2,2> purkinje_cable_assembler(&mesh);

		purkinje_cable_assembler.SetMatrixToAssemble(purkinje_mat, true);
		purkinje_cable_assembler.Assemble();
		PetscMatTools::Finalise(purkinje_mat);

		PetscInt lo, hi;
		PetscMatTools::GetOwnershipRange(purkinje_mat, lo, hi);
		for (AbstractTetrahedralMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
				node_iter != mesh.GetNodeIteratorEnd(); ++node_iter)
		{
			unsigned i = node_iter->GetIndex();
			assert(lo<=(int)(2*i) && (int)(2*i)<hi);
			for(unsigned j=0; j<mesh.GetNumNodes(); j++)
			{
				TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i,2*j), 0 , 1e-8);
				TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i,2*j+1), 0.0, 1e-8);
			}

			assert(lo<=(int)(2*i+1) && (int)(2*i+1)<hi);
			for(unsigned j=0; j<mesh.GetNumNodes(); j++)
			{
				//Non-Purkinje are all zero
				TS_ASSERT_DELTA( PetscMatTools::GetElement(purkinje_mat,2*i+1,2*j),   0.0, 1e-8);

				//Make sure that columns associated with cable node have non-zero Purkinje entries
				if ( (i>55) && (i<65) && (j>=i-1) && (j<=i+1))
				{
					TS_ASSERT_DIFFERS( PetscMatTools::GetElement(purkinje_mat,2*i+1,2*j+1), 0.0);
				}
				else if ((i==55) && (j>=55) && (j<=56) )
				{
					TS_ASSERT_DIFFERS( PetscMatTools::GetElement(purkinje_mat,2*i+1,2*j+1), 0.0);
				}
				else if ((i==65) && (j>=64) && (j<=65) )
				{
					TS_ASSERT_DIFFERS( PetscMatTools::GetElement(purkinje_mat,2*i+1,2*j+1), 0.0);
				}
				else
				{
					//Other entries are zero
					TS_ASSERT_DELTA(PetscMatTools::GetElement(purkinje_mat,2*i+1,2*j+1), 0.0 ,1.0e-8);
				}
			}
		}

		MatDestroy(purkinje_mat);
	}

    void TestMonodomainPurkinjeSolver() throw(Exception)
    {
        PdeSimulationTime::SetPdeTimeStep(0.01);

        std::string mesh_base("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
        TrianglesMeshReader<2,2> reader(mesh_base);

        MixedDimensionMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        PurkinjeCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        MonodomainTissue<2> tissue( &cell_factory );

        // don't
        BoundaryConditionsContainer<2,2,2> bcc;

        MonodomainPurkinjeSolver<2,2> solver(&mesh,&tissue,&bcc);

        ///\todo #1898 set up initial conditions vector, then call Solve()...

    }
};

#endif // TESTMONODOMAINPURKINJEASSEMBLERSANDSOLVER_HPP_
