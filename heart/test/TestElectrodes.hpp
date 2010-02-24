/*

Copyright (C) University of Oxford, 2005-2010

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

#ifndef TESTELECTRODES_HPP_
#define TESTELECTRODES_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <vector>

#include "Electrodes.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "DistributedVector.hpp"
#include "ArchiveOpener.hpp"

class TestElectrodes : public CxxTest::TestSuite
{
public:
    void TestElectrodeGrounded2dAndSwitchOffSwitchOn() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        DistributedTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double magnitude = 543.324;
        double start_time = 1.0; //ms
        double duration = 2.0; //ms
        Electrodes<2> electrodes(mesh,true,0,0.0,10.0,magnitude,start_time,duration);

        boost::shared_ptr<BoundaryConditionsContainer<2,2,2> > p_bcc = electrodes.GetBoundaryConditionsContainer();

        for(DistributedTetrahedralMesh<2,2>::BoundaryElementIterator iter
           = mesh.GetBoundaryElementIteratorBegin();
           iter != mesh.GetBoundaryElementIteratorEnd();
           iter++)
        {
            if ( fabs((*iter)->CalculateCentroid()[0] - 0.0) < 1e-6 )
            {
                double value = p_bcc->GetNeumannBCValue(*iter,(*iter)->CalculateCentroid(),1);
                TS_ASSERT_DELTA(value,magnitude,1e-12);
            }
        }

        unsigned num_grounded_nodes = 0u;

        for (AbstractTetrahedralMesh<2,2>::NodeIterator iter=mesh.GetNodeIteratorBegin();
             iter != mesh.GetNodeIteratorEnd();
             ++iter)
        {
            Node<2>* p_node = &(*iter);
            if (p_bcc->HasDirichletBoundaryCondition(p_node, 1))
            {
                double x_val = p_node->rGetLocation()[0];
                TS_ASSERT_DELTA(x_val, 10.0, 1e-12);
                num_grounded_nodes++;
                TS_ASSERT_EQUALS(p_bcc->GetDirichletBCValue(p_node, 1), 0.0);
            }
        }

        unsigned num_grounded_nodes_reduced;
        int mpi_ret = MPI_Allreduce(&num_grounded_nodes, &num_grounded_nodes_reduced, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        TS_ASSERT_EQUALS(mpi_ret, MPI_SUCCESS);

        TS_ASSERT_EQUALS(num_grounded_nodes_reduced, 11u);

        TS_ASSERT_THROWS_THIS(Electrodes<2> bad_electrodes(mesh,true,0,5.0,10.0,magnitude,start_time,duration),
                "Minimum value of coordinate is not the value given");
        TS_ASSERT_THROWS_THIS(Electrodes<2> bad_electrodes(mesh,true,0,0.0,30.0,magnitude,start_time,duration),
                "Maximum value of coordinate is not the value given");

        // Nothing at the beginning of the simulation
        TS_ASSERT_EQUALS(electrodes.SwitchOff(0.0), false); // t<end time
        TS_ASSERT_EQUALS(electrodes.SwitchOn(0.0), false); // t<start time
               
        // Time to switch on
        TS_ASSERT_EQUALS(electrodes.SwitchOff(1.0), false); // false as t<end time
        TS_ASSERT_EQUALS(electrodes.SwitchOn(1.0), true); // true as t>start time

        // Implemented to switch off at times extrictly bigger than starting point + duration
        TS_ASSERT_EQUALS(electrodes.SwitchOff(3.0), false); // false as t<end time
        TS_ASSERT_EQUALS(electrodes.SwitchOn(3.0),false); // false as electrode already switched on
        
        // Time to switch off
        TS_ASSERT_EQUALS(electrodes.SwitchOff(3.0+1e-12), true); // true as t>end_time
        TS_ASSERT_EQUALS(electrodes.SwitchOn(3.0+1e-12), false); // false as electrode already switched on

        // Everything is over now...
        TS_ASSERT_EQUALS(electrodes.SwitchOff(4.0), false); // false as electrodes has been switched off
        TS_ASSERT_EQUALS(electrodes.SwitchOn(4.0), false); // false as electrodes has been switched on
    }


    void TestElectrodeUngrounded2d() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        DistributedTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double magnitude = 543.324;
        double duration = 2.0;
        Electrodes<2> electrodes(mesh,false,0,0,10,magnitude,0.0,duration);

        boost::shared_ptr<BoundaryConditionsContainer<2,2,2> >  p_bcc = electrodes.GetBoundaryConditionsContainer();

        for(DistributedTetrahedralMesh<2,2>::BoundaryElementIterator iter
                = mesh.GetBoundaryElementIteratorBegin();
           iter != mesh.GetBoundaryElementIteratorEnd();
           iter++)
        {
            if ( fabs((*iter)->CalculateCentroid()[0] - 0.0) < 1e-6 )
            {
                double value = p_bcc->GetNeumannBCValue(*iter,(*iter)->CalculateCentroid(),1);

                TS_ASSERT_DELTA(value,magnitude,1e-12);
            }


            if ( fabs((*iter)->CalculateCentroid()[0] - 10.0) < 1e-6 )
            {
                double value = p_bcc->GetNeumannBCValue(*iter,(*iter)->CalculateCentroid(),1);
                TS_ASSERT_DELTA(value,-magnitude,1e-12);
            }
        }
    }

    void TestElectrodeUngrounded2dDifferentAreas() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        DistributedTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        try
        {
            // Node 120 is the bottom right corner (10,10) in the original mesh.
            unsigned index = 120;
            
            // In parallel work out the new index for 120
            if (!PetscTools::IsSequential())
            {
                index =  mesh.rGetNodePermutation()[120];
            }
            
            // Move the node slightlty to the left so it's not considered to be located at the edge by 
            // the Electrodes constructor. Node might be halo as well...
            c_vector<double, 2>& corner_node_location = mesh.GetNode(index)->rGetModifiableLocation();
            corner_node_location[0] -= 0.5;
        }
        catch(Exception& e)
        {
            // Don't do anything if you don't own (or halo-own) the node
        } 

        // Arbitrary input flux
        double flux_in_magnitude = 543.324;
        
        // area of the left electrode is 10, area of the right electrode is 9
        // input flux * left area = output flux * right area
        double flux_out_magnitude = flux_in_magnitude * 10/9;
                
        double duration = 2.0;
        Electrodes<2> electrodes(mesh,false,0,0,10,flux_in_magnitude,0.0,duration);

        TS_ASSERT_THROWS_ANYTHING(electrodes.ComputeElectrodesAreasAndCheckEquality(0,0,10));
        
        boost::shared_ptr<BoundaryConditionsContainer<2,2,2> >  p_bcc = electrodes.GetBoundaryConditionsContainer();

        for(DistributedTetrahedralMesh<2,2>::BoundaryElementIterator iter
                = mesh.GetBoundaryElementIteratorBegin();
           iter != mesh.GetBoundaryElementIteratorEnd();
           iter++)
        {
            if ( fabs((*iter)->CalculateCentroid()[0] - 0.0) < 1e-6 )
            {
                double value = p_bcc->GetNeumannBCValue(*iter,(*iter)->CalculateCentroid(),1);

                TS_ASSERT_DELTA(value,flux_in_magnitude,1e-12);
            }


            if ( fabs((*iter)->CalculateCentroid()[0] - 10.0) < 1e-6 )
            {
                double value = p_bcc->GetNeumannBCValue(*iter,(*iter)->CalculateCentroid(),1);
                TS_ASSERT_DELTA(value,-flux_out_magnitude,1e-12);
            }
        }
    }

    void TestElectrodeGrounded3d() throw (Exception)
    {
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(10,10,10);

        double magnitude = 543.324;
        double duration = 2.0;
        Electrodes<3> electrodes(mesh,true,1,0,10,magnitude,0.0,duration);

        boost::shared_ptr<BoundaryConditionsContainer<3,3,2> >  p_bcc = electrodes.GetBoundaryConditionsContainer();

        for(DistributedTetrahedralMesh<3,3>::BoundaryElementIterator iter
                = mesh.GetBoundaryElementIteratorBegin();
           iter != mesh.GetBoundaryElementIteratorEnd();
           iter++)
        {
            if ( fabs((*iter)->CalculateCentroid()[1] - 0.0) < 1e-6 )
            {
                double value = p_bcc->GetNeumannBCValue(*iter,(*iter)->CalculateCentroid(),1);

                TS_ASSERT_DELTA(value,magnitude,1e-12);
            }
        }

        unsigned num_grounded_nodes = 0u;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            Node<3>* p_node = mesh.GetNode(i);
            if (p_bcc->HasDirichletBoundaryCondition(p_node, 1))
            {
                double y_val = p_node->rGetLocation()[1];
                TS_ASSERT_DELTA(y_val, 10.0, 1e-12);
                num_grounded_nodes++;
                TS_ASSERT_EQUALS(p_bcc->GetDirichletBCValue(p_node, 1), 0.0);
            }
        }
        TS_ASSERT_EQUALS(num_grounded_nodes, 121u);
    }

    void TestElectrodeUngrounded3d() throw (Exception)
    {
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(10,10,10);

        double magnitude = 543.324;
        double duration = 2.0;
        Electrodes<3> electrodes(mesh,false,1,0,10,magnitude,0.0,duration);

        boost::shared_ptr<BoundaryConditionsContainer<3,3,2> > p_bcc = electrodes.GetBoundaryConditionsContainer();

        for(DistributedTetrahedralMesh<3,3>::BoundaryElementIterator iter
                = mesh.GetBoundaryElementIteratorBegin();
           iter != mesh.GetBoundaryElementIteratorEnd();
           iter++)
        {
            if ( fabs((*iter)->CalculateCentroid()[1] - 0.0) < 1e-6 )
            {
                double value = p_bcc->GetNeumannBCValue(*iter,(*iter)->CalculateCentroid(),1);

                TS_ASSERT_DELTA(value,magnitude,1e-12);
            }


            if ( fabs((*iter)->CalculateCentroid()[1] - 10.0) < 1e-6 )
            {
                double value = p_bcc->GetNeumannBCValue(*iter,(*iter)->CalculateCentroid(),1);
                TS_ASSERT_DELTA(value,-magnitude,1e-12);
            }
        }
    }
    
    /**
     * Note - we have had to make the Electrodes class archive the mesh as well so it 
     * can properly reconstruct itself. Becasue boost is clever it should still only 
     * archive one copy of the mesh in an archive file.
     */
    void TestArchivingElectrodes() throw(Exception)
    {
        std::string archive_dir = "archive";
        std::string archive_file = "Electrodes.arch";
        
        // These values used for construction and testing later.
        double magnitude = 543.324;
        double duration = 2.0;
        
        {
            // Make a mesh
            TetrahedralMesh<3,3> mesh;
            mesh.ConstructCuboid(10,10,10);
            // Create Electrodes class - the const is required to prevent boost errors on compilation.
            Electrodes<3>* const p_electrodes = new Electrodes<3>(mesh,false,1,0,10,magnitude,0.0,duration);
            
            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Save the electrodes to the archive...   
            (*p_arch) << p_electrodes;
            delete p_electrodes;
        }
    
        {
            Electrodes<3>* p_electrodes;
    
            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
    
            // restore from the archive
            (*p_arch) >> p_electrodes;
    
            // Repeat the above test TestElectrodeUngrounded3d() on the loaded Electrodes object.    
            boost::shared_ptr<BoundaryConditionsContainer<3,3,2> > p_bcc = 
                p_electrodes->GetBoundaryConditionsContainer();
        
            // (We can do this p_electrodes->mpMesh only because we are a friend class,
            // this is the only place it is used so didn't want to make a GetMesh() method).
            for(DistributedTetrahedralMesh<3,3>::BoundaryElementIterator iter 
                     = p_electrodes->mpMesh->GetBoundaryElementIteratorBegin();
               iter != p_electrodes->mpMesh->GetBoundaryElementIteratorEnd();
               iter++)
            {
                if ( fabs((*iter)->CalculateCentroid()[1] - 0.0) < 1e-6 )
                {
                    double value = p_bcc->GetNeumannBCValue(*iter,(*iter)->CalculateCentroid(),1);    
                    TS_ASSERT_DELTA(value,magnitude,1e-12);
                }   
    
                if ( fabs((*iter)->CalculateCentroid()[1] - 10.0) < 1e-6 )
                {
                    double value = p_bcc->GetNeumannBCValue(*iter,(*iter)->CalculateCentroid(),1);
                    TS_ASSERT_DELTA(value,-magnitude,1e-12);
                }
            }    
            
            // Free memory
            delete p_electrodes->mpMesh;
            delete p_electrodes;
        }
    }
    
};



#endif /*TESTELECTRODES_HPP_*/
