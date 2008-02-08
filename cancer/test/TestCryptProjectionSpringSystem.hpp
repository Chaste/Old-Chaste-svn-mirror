#ifndef TESTCRYPTPROJECTIONSPRINGSYSTEM_HPP_
#define TESTCRYPTPROJECTIONSPRINGSYSTEM_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <cmath>
#include <vector>

#include "CryptProjectionSpringSystem.hpp"
#include "Meineke2001SpringSystem.hpp"
#include "FixedCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCancerTestSuite.hpp"

/**
 * Note that all these tests call setUp() and tearDown() before running,
 * so if you copy them into a new test suite be sure to copy these methods
 * too.
 */
class TestCryptProjectionSpringSystem : public AbstractCancerTestSuite
{
public:   
    
    void TestSpringSystemMethods() throw (Exception)
    {      
        CancerParameters* p_params = CancerParameters::Instance();

        // Create a mesh 
        unsigned num_cells_width = 10;
        unsigned num_cells_depth = 10;
        unsigned thickness_of_ghost_layer = 0;
        
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, thickness_of_ghost_layer, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
               
        // Centre the mesh at (0,0)
        c_vector<double,2> width_extremes = p_mesh->GetWidthExtremes(0u);
        c_vector<double,2> height_extremes = p_mesh->GetWidthExtremes(1u);  
              
        double width_of_mesh = (num_cells_width/(num_cells_width + 2.0*thickness_of_ghost_layer))*(width_extremes[1] - width_extremes[0]);
        double height_of_mesh = (num_cells_depth/(num_cells_depth + 2.0*thickness_of_ghost_layer))*(height_extremes[1] - height_extremes[0]);
                
        p_mesh->Translate(-width_of_mesh/2, -height_of_mesh/2);
                
        // Create some cells
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            
            if (i==4 || i==5)
            {
                cell.SetBirthTime(-0.5);
            }
            else
            {
                cell.SetBirthTime(-10.0);
            }
            cells.push_back(cell);
        }

        // Create a tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        tissue.MarkSpring(tissue.rGetCellAtNodeIndex(4), tissue.rGetCellAtNodeIndex(5));
        tissue.SetGhostNodes(ghost_node_indices);
        
        // Test it is possible to construct a spring system
        TS_ASSERT_THROWS_NOTHING(CryptProjectionSpringSystem spring_system(tissue));
                
        // Create a spring system with crypt surface z = 2*r
        p_params->SetCryptProjectionParameterA(2.0);
        p_params->SetCryptProjectionParameterB(1.0);
        CryptProjectionSpringSystem spring_system(tissue);
        
        TS_ASSERT(!spring_system.NeedsVoronoiTessellation()) // for coverage        
                
        // Test get methods
        TS_ASSERT_DELTA(spring_system.GetA(), 2.0, 1e-12);
        TS_ASSERT_DELTA(spring_system.GetB(), 1.0, 1e-12);      
        
        // Test crypt height and gradient calculations
        c_vector<double, 2> node_location_2d = p_mesh->GetNode(0)->rGetLocation();  
        TS_ASSERT_DELTA(spring_system.CalculateCryptSurfaceHeightAtPoint(node_location_2d), 2.0*pow(norm_2(node_location_2d),1.0), 1e-12);
        TS_ASSERT_DELTA(spring_system.CalculateCryptSurfaceDerivativeAtPoint(node_location_2d), 2.0, 1e-12);
                        
        // Test updating of mNode3dLocationMap         
        spring_system.UpdateNode3dLocationMap();
                
        // Move a node slightly
        ChastePoint<2> new_point;
        new_point.rGetLocation()[0] = node_location_2d[0]+0.05;
        new_point.rGetLocation()[1] = node_location_2d[1];  
        p_mesh->SetNode(0, new_point, false);
        
        
        // Test UpdateNode3dLocationMap()       
              
        c_vector<double, 2> new_node_location_2d;
        new_node_location_2d[0] = new_point.rGetLocation()[0];
        new_node_location_2d[1] = new_point.rGetLocation()[1];  
              
        spring_system.UpdateNode3dLocationMap();         
        
        // Check the map updates correctly (note that we have used no ghost nodes, so the map does contain 0)        
        c_vector<double, 3> calculated_new_node_location_3d = spring_system.mNode3dLocationMap[0];        
        c_vector<double, 3> correct_new_node_location_3d;
        
        correct_new_node_location_3d[0] = new_node_location_2d[0];
        correct_new_node_location_3d[1] = new_node_location_2d[1];
        correct_new_node_location_3d[2] = spring_system.CalculateCryptSurfaceHeightAtPoint(new_node_location_2d);
         
        TS_ASSERT_DELTA(calculated_new_node_location_3d[0], correct_new_node_location_3d[0], 1e-12);
        TS_ASSERT_DELTA(calculated_new_node_location_3d[1], correct_new_node_location_3d[1], 1e-12);
        TS_ASSERT_DELTA(calculated_new_node_location_3d[2], correct_new_node_location_3d[2], 1e-12);
        
        
        // Test force calculation on a normal spring
                      
        c_vector<double,2> force_on_spring; // between nodes 0 and 1
        
        // Find one of the elements that nodes 0 and 1 live on
        ChastePoint<2> new_point2;
        new_point2.rGetLocation()[0] = new_point[0] + 0.01;
        new_point2.rGetLocation()[1] = new_point[1] + 0.01 ;
        
        unsigned elem_index = p_mesh->GetContainingElementIndex(new_point2, false);
        Element<2,2>* p_element = p_mesh->GetElement(elem_index);
        
        force_on_spring = spring_system.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),p_element->GetNodeGlobalIndex(0));
              
        TS_ASSERT_DELTA(force_on_spring[0], -5.7594, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0230 , 1e-4);   
            

        // Test force calculation with a cutoff
        
        double dist = norm_2( p_mesh->GetVectorFromAtoB(p_element->GetNode(0)->rGetLocation(), p_element->GetNode(1)->rGetLocation()) );   
        spring_system.UseCutoffPoint(dist-0.1);

        force_on_spring = spring_system.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),p_element->GetNodeGlobalIndex(0));
        TS_ASSERT_DELTA(force_on_spring[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);
        
        
        // Test force calculation for a pair of newly born neighbouring cells 
        force_on_spring = spring_system.CalculateForceBetweenNodes(4,5);               
        TS_ASSERT_DELTA(force_on_spring[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0 , 1e-4);   

        tissue.UnmarkSpring(tissue.rGetCellAtNodeIndex(4), tissue.rGetCellAtNodeIndex(5));

        // Test velocity calculation for a particular node
        std::vector<c_vector<double, 2> >& velocities_on_each_node = spring_system.rCalculateVelocitiesOfEachNode();
                
        TS_ASSERT_DELTA(velocities_on_each_node[0][0], 0.0, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[0][1], 0.0, 1e-4);
                       
                
        // Test that in the case of a flat crypt surface (mA=mB=0), the results are the same as for Meineke2001SpringSystem
        p_params->SetCryptProjectionParameterA(0.001);
        p_params->SetCryptProjectionParameterB(0.001);
        CryptProjectionSpringSystem flat_crypt_spring_system(tissue);
        Meineke2001SpringSystem<2> meineke_spring_system(tissue);
        
        // Normally this would be set up at the start of rCalculateVelocitiesOfEachNode
        flat_crypt_spring_system.UpdateNode3dLocationMap();
                
        for(MeshBasedTissue<2>::SpringIterator spring_iterator = tissue.SpringsBegin();
            spring_iterator != tissue.SpringsEnd();
            ++spring_iterator)
        {        
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            
            c_vector<double, 2> force_flat = flat_crypt_spring_system.CalculateForceBetweenNodes(nodeA_global_index,nodeB_global_index);
            c_vector<double, 2> force_meineke = meineke_spring_system.CalculateForceBetweenNodes(nodeA_global_index,nodeB_global_index);
            
            TS_ASSERT_DELTA( force_flat[0], force_meineke[0], 1e-3);
            TS_ASSERT_DELTA( force_flat[1], force_meineke[1], 1e-3);
        }
    }
    
    
    void xTestWntChemotaxis() throw (Exception)
    {      
        CancerParameters* p_params = CancerParameters::Instance();
        
        // Create a mesh 
        unsigned num_cells_width = 10;
        unsigned num_cells_depth = 10;
        unsigned thickness_of_ghost_layer = 0;
        
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, thickness_of_ghost_layer, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
               
        // Centre the mesh at (0,0)
        c_vector<double,2> width_extremes = p_mesh->GetWidthExtremes(0u);
        c_vector<double,2> height_extremes = p_mesh->GetWidthExtremes(1u);  
              
        double width_of_mesh = (num_cells_width/(num_cells_width + 2.0*thickness_of_ghost_layer))*(width_extremes[1] - width_extremes[0]);
        double height_of_mesh = (num_cells_depth/(num_cells_depth + 2.0*thickness_of_ghost_layer))*(height_extremes[1] - height_extremes[0]);
                
        p_mesh->Translate(-width_of_mesh/2, -height_of_mesh/2);
                
        // Create some cells
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(-10.0);
            cells.push_back(cell);
        }

        // Create a tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        tissue.MarkSpring(tissue.rGetCellAtNodeIndex(4), tissue.rGetCellAtNodeIndex(5));
        tissue.SetGhostNodes(ghost_node_indices);
        
        WntConcentration::Instance()->SetType(RADIAL); 
        WntConcentration::Instance()->SetTissue(tissue);  
                
        // Create a spring system with crypt surface z = 2*r
        p_params->SetCryptProjectionParameterA(2.0);
        p_params->SetCryptProjectionParameterB(1.0);
        CryptProjectionSpringSystem spring_system(tissue);
        
        spring_system.SetWntChemotaxis(true);
        
        // Test velocity calculation for a particular node
        std::vector<c_vector<double, 2> >& velocities_on_each_node = spring_system.rCalculateVelocitiesOfEachNode();
                
        TS_ASSERT_DELTA(velocities_on_each_node[0][0], 0.0, 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[0][1], 0.0, 1e-4);
       
    }
    
    
    
    void TestArchiving() throw (Exception)
    {   
        CancerParameters* p_params = CancerParameters::Instance();
        
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "crypt_projection_spring_system.arch";

        unsigned num_nodes;
        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        
            ConformingTetrahedralMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            num_nodes = mesh.GetNumNodes();

            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

            std::vector<TissueCell> cells;
            TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
            for(unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                cell.SetNodeIndex(i);
                cell.SetBirthTime(-50.0);
                cells.push_back(cell);
            }
        
            MeshBasedTissue<2> crypt(mesh,cells);
            p_params->SetCryptProjectionParameterA(1.0);
            p_params->SetCryptProjectionParameterB(2.0);
            CryptProjectionSpringSystem spring_system(crypt);
         
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            CryptProjectionSpringSystem * const p_spring_system = &spring_system;  
            
            p_spring_system->UseCutoffPoint(1.1);

            output_arch << p_spring_system;
        }
       
        {
            MeshBasedTissue<2>::meshPathname = "mesh/test/data/square_2_elements";
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            CryptProjectionSpringSystem* p_spring_system;
            
            // Restore from the archive
            input_arch >> p_spring_system;
            
            // Test the member data
            TS_ASSERT_EQUALS(p_spring_system->mUseCutoffPoint,true);
            TS_ASSERT_DELTA(p_spring_system->mCutoffPoint,1.1,1e-12);  
            TS_ASSERT_DELTA(p_spring_system->GetA(),1.0,1e-12); 
            TS_ASSERT_DELTA(p_spring_system->GetB(),2.0,1e-12);
            
            delete p_spring_system->mpTissue;
            delete p_spring_system;
        }
    } 

};


#endif /*TESTCRYPTPROJECTIONSPRINGSYSTEM_HPP_*/
