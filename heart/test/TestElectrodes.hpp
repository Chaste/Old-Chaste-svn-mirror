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

#ifndef TESTELECTRODES_HPP_
#define TESTELECTRODES_HPP_


#include <cxxtest/TestSuite.h>
#include <vector>
#include "Electrodes.hpp"


class TestElectrodes : public CxxTest::TestSuite
{
public: 
    void TestElectrodeGrounded2d() throw (Exception)
    {
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(10,10);
        
        double magnitude = 543.324;
        double duration = 2.0;        
        Electrodes<2> electrodes(mesh,true,0,0.0,10.0,magnitude,duration);
        
        BoundaryConditionsContainer<2,2,2>* p_bcc = electrodes.GetBoundaryConditionsContainer();
        
        for(TetrahedralMesh<2,2>::BoundaryElementIterator iter 
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
        
        std::vector<unsigned> grounded_nodes = electrodes.GetGroundedNodes();
        TS_ASSERT_EQUALS(grounded_nodes.size(), 11u);
        for(unsigned i=0; i<grounded_nodes.size(); i++)
        {
            double x_val = mesh.GetNode(grounded_nodes[i])->rGetLocation()[0];
            TS_ASSERT_DELTA(x_val, 10.0, 1e-12);
        }
    }



    void TestElectrodeUngrounded2d() throw (Exception)
    {
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(10,10);
        
        double magnitude = 543.324;
        double duration = 2.0;
        Electrodes<2> electrodes(mesh,false,0,0,10,magnitude,duration);
        
        BoundaryConditionsContainer<2,2,2>* p_bcc = electrodes.GetBoundaryConditionsContainer();
        
        for(TetrahedralMesh<2,2>::BoundaryElementIterator iter 
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

    void TestElectrodeGrounded3d() throw (Exception)
    {
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(10,10,10);
        
        double magnitude = 543.324;
        double duration = 2.0;
        Electrodes<3> electrodes(mesh,true,1,0,10,magnitude,duration);
        
        BoundaryConditionsContainer<3,3,2>* p_bcc = electrodes.GetBoundaryConditionsContainer();
        
        for(TetrahedralMesh<3,3>::BoundaryElementIterator iter 
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

        std::vector<unsigned> grounded_nodes = electrodes.GetGroundedNodes();
        TS_ASSERT_EQUALS(grounded_nodes.size(), 121u);
        for(unsigned i=0; i<grounded_nodes.size(); i++)
        {
            double y_val = mesh.GetNode(grounded_nodes[i])->rGetLocation()[1];
            TS_ASSERT_DELTA(y_val, 10.0, 1e-12);
        }
    }
    
    void TestElectrodeUngrounded3d() throw (Exception)
    {
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(10,10,10);
        
        double magnitude = 543.324;
        double duration = 2.0;
        Electrodes<3> electrodes(mesh,false,1,0,10,magnitude,duration);
        
        BoundaryConditionsContainer<3,3,2>* p_bcc = electrodes.GetBoundaryConditionsContainer();
        
        for(TetrahedralMesh<3,3>::BoundaryElementIterator iter 
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
};
  


#endif /*TESTELECTRODES_HPP_*/
