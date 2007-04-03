#ifndef TESTSPRINGS3D_HPP_
#define TESTSPRINGS3D_HPP_


#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "TrianglesMeshWriter.cpp"
#include <cmath>

#include <vector>
#include "OutputFileHandler.hpp"

#include "CancerParameters.hpp"
#include "ColumnDataReader.hpp"
#include "CancerParameters.hpp"
#include "NodeMap.hpp"


class TestSprings3d : public CxxTest::TestSuite
{
    
   ConformingTetrahedralMesh<3,3> Make3dMesh()
   {
       ConformingTetrahedralMesh<3,3> mesh;
       unsigned width=2;
       unsigned height=2;
       unsigned depth=1;

//     unsigned num_boundary_nodes =   2*( (width+1)*(height+1) + (width+1)*(depth+1) + (depth+1)*(height+1) )
//                                     - 4*(width-1 + height-1 + depth-1)
//                                     - 16;

       mesh.ConstructCuboid(width,height,depth,true);

       TrianglesMeshWriter<3,3> mesh_writer("","3dSpringMesh");
       mesh_writer.WriteFilesUsingMesh(mesh);

       return mesh;

   }
   
   /**
    * Calculates the forces on each node
    *
    * @return drdt the x and y force components on each node
    */
    std::vector<std::vector<double> > CalculateForcesOnEachNode(ConformingTetrahedralMesh<3,3>& rMesh)
    {
        std::vector<std::vector<double> > drdt(rMesh.GetNumAllNodes());
        for (unsigned i=0; i<rMesh.GetNumAllNodes(); i++)
        {
            drdt[i].resize(3);
        }
        
        std::vector<std::vector<unsigned> > node_pairs_checked;
        //////////////////////////////////////////////////////////////////
        // loop over element and for each one loop over its three edges
        ////////////////////////////////////////////////////////////////////
        for (unsigned elem_index = 0; elem_index<rMesh.GetNumAllElements(); elem_index++)
        {
            Element<3,3>* p_element = rMesh.GetElement(elem_index);
            
            for (unsigned k=0; k<4; k++)
            {
                unsigned nodeA = k;
                
                for (unsigned l=k+1; l<k+4; l++)
                {
                    unsigned nodeB = l%4;
                    unsigned nodeA_global_index = p_element->GetNode(nodeA)->GetIndex();
                    unsigned nodeB_global_index = p_element->GetNode(nodeB)->GetIndex();
                
                
                    // check whether we have already worked out the force between these two...
                    bool is_force_already_calculated = false;
                    
                    for (unsigned i=0 ; i<node_pairs_checked.size() ; i++)
                    {
                       std::vector<unsigned> node_pair = node_pairs_checked[i];
                       if(node_pair[0]==nodeA_global_index || node_pair[1]==nodeA_global_index)
                       { // first node is in node_pair
                            if(node_pair[0]==nodeB_global_index || node_pair[1]==nodeB_global_index)
                            { // both are in node_pair
                                is_force_already_calculated = true;
                            } 
                       } 
                    }                
                
                    if (!is_force_already_calculated)
                    {
                        c_vector<double, 3> drdt_contribution = CalculateForceInThisSpring(p_element,nodeA,nodeB);
                        std::vector<unsigned> this_pair;
                        this_pair.push_back(nodeA_global_index);
                        this_pair.push_back(nodeB_global_index);
                        node_pairs_checked.push_back(this_pair);
                        
                        drdt[ nodeB_global_index ][0] -= drdt_contribution(0);
                        drdt[ nodeB_global_index ][1] -= drdt_contribution(1);
                        drdt[ nodeB_global_index ][2] -= drdt_contribution(2);
                                
                        drdt[ nodeA_global_index ][0] += drdt_contribution(0);
                        drdt[ nodeA_global_index ][1] += drdt_contribution(1);
                        drdt[ nodeA_global_index ][2] += drdt_contribution(2);
                    }
                }
            }
        }
    
        return drdt;
    }

    /**
     * @return the x and y forces in this spring
     */
    c_vector<double, 3> CalculateForceInThisSpring(Element<3,3>*& rPElement,const unsigned& rNodeA,const unsigned& rNodeB)
    {
        c_vector<double, 3> drdt_contribution;
        c_vector<double, 3> unit_difference;
        unit_difference(0)=rPElement->GetNodeLocation(rNodeB,0)-rPElement->GetNodeLocation(rNodeA,0);
        unit_difference(1)=rPElement->GetNodeLocation(rNodeB,1)-rPElement->GetNodeLocation(rNodeA,1);
        unit_difference(2)=rPElement->GetNodeLocation(rNodeB,2)-rPElement->GetNodeLocation(rNodeA,2);
        double distance_between_nodes=sqrt(unit_difference(0)*unit_difference(0)+unit_difference(1)*unit_difference(1)+unit_difference(2)*unit_difference(2));
        unit_difference=unit_difference/distance_between_nodes;
        
        double rest_length = 1.0;
        
        CancerParameters *p_params = CancerParameters::Instance();
        
        return drdt_contribution = p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal() * unit_difference * (distance_between_nodes - rest_length);
    }
    
    /**
     * Moves each node to a new position for this timestep
     *
     * @param rDrDt the x, y and z force components on each node.
     */
    void UpdateNodePositions(const std::vector< std::vector<double> >& rDrDt, ConformingTetrahedralMesh<3,3>& rMesh)
    {
        for (unsigned index = 0; index<rMesh.GetNumAllNodes(); index++)
        {
            Point<3> new_point = GetNewNodeLocation(index,rDrDt,rMesh);
            rMesh.SetNode(index, new_point, false);      
        }
        
        
    }
    
    Point<3> GetNewNodeLocation(const unsigned& rOldNodeIndex, const std::vector< std::vector<double> >& rDrDt, ConformingTetrahedralMesh<3,3>& rMesh)
    {
        double dt = 0.01;
        
        c_vector<double,3> old_point = rMesh.GetNode(rOldNodeIndex)->rGetLocation();
        Point<3> new_point;
        
        // Euler style update to node position
        new_point.rGetLocation()[0] = old_point[0] + dt*rDrDt[rOldNodeIndex][0];
        new_point.rGetLocation()[1] = old_point[1] + dt*rDrDt[rOldNodeIndex][1];
        new_point.rGetLocation()[2] = old_point[2] + dt*rDrDt[rOldNodeIndex][2];
        
        return new_point;
    }
   

public:

   void Test3dSpringMesh() throw (Exception)
   {

        ConformingTetrahedralMesh<3,3> mesh = Make3dMesh();
        
//        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
//        ConformingTetrahedralMesh<3,3> mesh;
//        mesh.ConstructFromMeshReader(mesh_reader);
//        
        
        {
            TrianglesMeshWriter<3,3> mesh_writer("","3dSpringMeshStart");
            mesh_writer.WriteFilesUsingMesh(mesh);

        }
        
        
        
        for (unsigned i=0 ; i<2 ; i++)
        {
//            std::cout << "step = " << i << std::endl;
            std::vector<std::vector<double> > forces = CalculateForcesOnEachNode(mesh);
            UpdateNodePositions(forces,mesh);
            
//            NodeMap map(mesh.GetNumAllNodes());
//            mesh.ReMesh(map);

//            std::ostringstream time_stamp;
//            time_stamp << i;
//            TrianglesMeshWriter<3,3> mesh_writer("","3dSpringMesh"+time_stamp.str());
//            mesh_writer.WriteFilesUsingMesh(mesh);
        }
        
        
        
        TrianglesMeshWriter<3,3> mesh_writer("","3dSpringMeshEnd");
        mesh_writer.WriteFilesUsingMesh(mesh);
        
        /*
         * Use the command:
         * 
         * tetview 3dSpringMesh
         * 
         * to visualize final output.
         */
        
   }
void TestOne3dElement() throw (Exception)
   {
        
//        ConformingTetrahedralMesh<3,3> mesh;
        
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        
        TrianglesMeshWriter<3,3> mesh_writer1("","3dSpringTetrahedronMeshStart");
        mesh_writer1.WriteFilesUsingMesh(mesh);
        
               
        // Test forces on springs
        unsigned nodeA = 0, nodeB = 1 ;
        Element<3,3>* p_element = mesh.GetElement(0);
        c_vector<double, 3> drdt_contribution = CalculateForceInThisSpring(p_element,nodeA,nodeB);
        
        for(unsigned i=0; i < 3;i++)
        {
            TS_ASSERT_DELTA(drdt_contribution[i],0.0,1e-6);
        }
        
        //  Test forces on nodes        
        for (unsigned i=0 ; i<1 ; i++)
        {
//            std::cout << "step = " << i << std::endl;
            std::vector<std::vector<double> > forces = CalculateForcesOnEachNode(mesh);
            UpdateNodePositions(forces,mesh);
            
            for (unsigned j=0; j<4; j++)
            {
                for(unsigned k=0;k<3;k++)
                {
                    TS_ASSERT_DELTA(forces[j][k],0.0,1e-6);
                }
            }
            
                     
            
//            NodeMap map(mesh.GetNumAllNodes());
//            mesh.ReMesh(map);

//            std::ostringstream time_stamp;
//            time_stamp << i;
//            TrianglesMeshWriter<3,3> mesh_writer("","3dSpringMesh"+time_stamp.str());
//            mesh_writer.WriteFilesUsingMesh(mesh);
        }
        TrianglesMeshWriter<3,3> mesh_writer("","3dSpringTetrahedronMeshEnd");
        mesh_writer.WriteFilesUsingMesh(mesh);    
        
        // Scale entire mesh and check that forces are correctly calculated
        CancerParameters *p_params = CancerParameters::Instance();
        double scale_factor = 1.5;
        
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            c_vector<double,3> old_point = mesh.GetNode(i)->rGetLocation();
            Point<3> new_point;
            new_point.rGetLocation()[0] = scale_factor*old_point[0];
            new_point.rGetLocation()[1] = scale_factor*old_point[1];
            new_point.rGetLocation()[2] = scale_factor*old_point[2];
            mesh.SetNode(i, new_point, false);            
        }
               
        std::vector<std::vector<double> > new_forces = CalculateForcesOnEachNode(mesh);
        UpdateNodePositions(new_forces,mesh);
        for (unsigned j=0; j<4; j++)
        {
            for(unsigned k=0;k<3;k++)
            {
                TS_ASSERT_DELTA(fabs(new_forces[j][k]),p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal()*(scale_factor-1)*sqrt(2),1e-6);
            }
        }
        
        // Move one node and check that forces are correctly calculated
        ConformingTetrahedralMesh<3,3> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader);
        
        c_vector<double,3> old_point = mesh2.GetNode(0)->rGetLocation();
        Point<3> new_point;
        new_point.rGetLocation()[0] = 0.0;
        new_point.rGetLocation()[1] = 0.0;
        new_point.rGetLocation()[2] = 0.0;
        mesh2.SetNode(0, new_point, false);   
          
        
        
         // Test forces on springs
        unsigned nodeA2 = 0, nodeB2 = 1 ;
        Element<3,3>* p_element2 = mesh2.GetElement(0);
        c_vector<double, 3> drdt_contribution2 = CalculateForceInThisSpring(p_element2,nodeA2,nodeB2);
        
        for(unsigned i=0; i < 3;i++)
        {
            TS_ASSERT_DELTA(fabs(drdt_contribution2[i]),p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal()*(1 - sqrt(3)/(2*sqrt(2)))/sqrt(3.0),1e-6);
        }
        
        std::vector<std::vector<double> > new_forces2 = CalculateForcesOnEachNode(mesh2);
        
        for(unsigned k=0;k<3;k++)
        {
            TS_ASSERT_DELTA(new_forces2[0][k],p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal()*(1 - sqrt(3)/(2*sqrt(2)))/sqrt(3.0),1e-6);
        }
        UpdateNodePositions(new_forces2,mesh2);
        TrianglesMeshWriter<3,3> mesh_writer2("","3dSpringTetrahedronMeshEnd2");
        mesh_writer2.WriteFilesUsingMesh(mesh);  
        
   }
};

#endif /*TESTSPRINGS3D_HPP_*/

