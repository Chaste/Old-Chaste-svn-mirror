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
#include "SpheroidSimulation3D.hpp"
#include "NodeMap.hpp"

#include "MeinekeCryptCell.hpp"
#include "FixedCellCycleModel.hpp"
//#include "StochasticCellCycleModel.hpp"
//#include "WntCellCycleModel.hpp"
//#include "WntGradient.hpp"
//#include "WntCellCycleOdeSystem.hpp"
//#include "TysonNovakCellCycleModel.hpp"
#include "ColumnDataReader.hpp"
#include "SimulationTime.hpp"

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
   
   
    

   

public:

//   void Test3dSpringMesh() throw (Exception)
//   {
//
//        ConformingTetrahedralMesh<3,3> mesh = Make3dMesh();
//        
////        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
////        ConformingTetrahedralMesh<3,3> mesh;
////        mesh.ConstructFromMeshReader(mesh_reader);
////        
//        
//        {
//            TrianglesMeshWriter<3,3> mesh_writer("","3dSpringMeshStart");
//            mesh_writer.WriteFilesUsingMesh(mesh);
//
//        }
//        
//        
//        
//        for (unsigned i=0 ; i<2 ; i++)
//        {
////            std::cout << "step = " << i << std::endl;
//            std::vector<std::vector<double> > forces = CalculateForcesOnEachNode(mesh);
//            UpdateNodePositions(forces,mesh);
//            
////            NodeMap map(mesh.GetNumAllNodes());
////            mesh.ReMesh(map);
//
////            std::ostringstream time_stamp;
////            time_stamp << i;
////            TrianglesMeshWriter<3,3> mesh_writer("","3dSpringMesh"+time_stamp.str());
////            mesh_writer.WriteFilesUsingMesh(mesh);
//        }
//        
//        
//        
//        TrianglesMeshWriter<3,3> mesh_writer("","3dSpringMeshEnd");
//        mesh_writer.WriteFilesUsingMesh(mesh);
//        
//        /*
//         * Use the command:
//         * 
//         * tetview 3dSpringMesh
//         * 
//         * to visualize final output.
//         */
//        
//   }
void TestOne3dElement() throw (Exception)
   {
        
//        ConformingTetrahedralMesh<3,3> mesh;
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element");
        
        ConformingTetrahedralMesh<3,3> mesh;
        
        
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TrianglesMeshWriter<3,3> mesh_writer1("","3dSpringTetrahedronMeshStart");
        mesh_writer1.WriteFilesUsingMesh(mesh);
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        std::vector<MeinekeCryptCell> cells;
        CryptCellType cell_type;
        cell_type = STEM;
        
        MeinekeCryptCell cell(cell_type, HEALTHY, 0u, new FixedCellCycleModel());
        for(unsigned i=0; i< mesh.GetNumNodes(); i++)
        {
            cell.SetNodeIndex(i);
            cell.SetBirthTime(-50.0);
            cells.push_back(cell);
        }
        
        
        
        SpheroidSimulation3D simulator = SpheroidSimulation3D(mesh,cells);
        simulator.SetMaxCells(40);
        simulator.SetMaxElements(40);
        //simulator.SetPeriodicSides(false);
        //simulator.SetOutputDirectory("TestPrivateMemberDirectory");
        // Test forces on springs
        unsigned nodeA = 0, nodeB = 1 ;
        Element<3,3>* p_element = mesh.GetElement(0);
        c_vector<double, 3> drdt_contribution = simulator.CalculateForceInThisSpring(p_element,nodeA,nodeB);
        for(unsigned i=0; i < 3;i++)
        {
            TS_ASSERT_DELTA(drdt_contribution[i],0.0,1e-6);
        }
        
        //  Test forces on nodes        
        for (unsigned i=0 ; i<1 ; i++)
        {
//            std::cout << "step = " << i << std::endl;
            std::vector<std::vector<double> > forces = simulator.CalculateVelocitiesOfEachNode();
            simulator.UpdateNodePositions(forces);
            
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
               
        std::vector<std::vector<double> > new_forces = simulator.CalculateVelocitiesOfEachNode();
        simulator.UpdateNodePositions(new_forces);
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
        SpheroidSimulation3D simulator2 = SpheroidSimulation3D(mesh2,cells);
        simulator2.SetMaxCells(40);
        simulator2.SetMaxElements(40);
        
        c_vector<double,3> old_point = mesh2.GetNode(0)->rGetLocation();
        Point<3> new_point;
        new_point.rGetLocation()[0] = 0.0;
        new_point.rGetLocation()[1] = 0.0;
        new_point.rGetLocation()[2] = 0.0;
        mesh2.SetNode(0, new_point, false);   
          
        
        
         // Test forces on springs
        unsigned nodeA2 = 0, nodeB2 = 1 ;
        Element<3,3>* p_element2 = mesh2.GetElement(0);
        c_vector<double, 3> drdt_contribution2 = simulator2.CalculateForceInThisSpring(p_element2,nodeA2,nodeB2);
        
        for(unsigned i=0; i < 3;i++)
        {
            TS_ASSERT_DELTA(fabs(drdt_contribution2[i]),p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal()*(1 - sqrt(3)/(2*sqrt(2)))/sqrt(3.0),1e-6);
        }
        
        std::vector<std::vector<double> > new_forces2 = simulator2.CalculateVelocitiesOfEachNode();
        
        for(unsigned k=0;k<3;k++)
        {
            TS_ASSERT_DELTA(new_forces2[0][k],p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal()*(1 - sqrt(3)/(2*sqrt(2)))/sqrt(3.0),1e-6);
        }
        simulator2.UpdateNodePositions(new_forces2);
        TrianglesMeshWriter<3,3> mesh_writer2("","3dSpringTetrahedronMeshEnd2");
        mesh_writer2.WriteFilesUsingMesh(mesh2);  
        
   }

};

#endif /*TESTSPRINGS3D_HPP_*/

