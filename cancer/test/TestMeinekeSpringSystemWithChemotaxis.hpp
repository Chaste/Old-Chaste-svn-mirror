#ifndef TESTMEINEKESPRINGSYSTEMWITHCHEMOTAXIS_HPP_
#define TESTMEINEKESPRINGSYSTEMWITHCHEMOTAXIS_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <cmath>
#include <vector>

#include "MeinekeSpringSystemWithChemotaxis.hpp"
#include "MeshBasedTissueWithGhostNodes.cpp"
#include "TrianglesMeshReader.hpp"
#include "CellsGenerator.hpp"
#include "AbstractCancerTestSuite.hpp"


class TestMeinekeSpringSystemWithChemotaxis : public AbstractCancerTestSuite
{    
public:

    void TestBasicChemotaxis() throw (Exception)
    {
        unsigned cells_across = 7;
        unsigned cells_up = 5;
        unsigned thickness_of_ghost_layer = 0;
        
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(-10);
            cells.push_back(cell);
        }

        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, ghost_node_indices);
        
        // Set up cellwisedata and associate it with the tissue
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(p_mesh->GetNumNodes(), 1);
        p_data->SetTissue(tissue);
        
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            double x = p_mesh->GetNode(i)->rGetLocation()[0];
            p_data->SetValue(x/50.0, p_mesh->GetNode(i));
        }
        
        MeinekeSpringSystemWithChemotaxis<2> meineke_spring_system_with_chemotaxis(tissue);
        std::vector<c_vector<double, 2> >& velocities_on_each_node = meineke_spring_system_with_chemotaxis.rCalculateVelocitiesOfEachNode();
 
        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            bool is_a_ghost_node = tissue.rGetGhostNodes()[i];

            if (!is_a_ghost_node)
            {
                double x = p_mesh->GetNode(i)->rGetLocation()[0];
                double c = x/50;
                double norm_grad_c = 1.0/50.0;
                double force_magnitude = meineke_spring_system_with_chemotaxis.ChemotacticForceMagnitude(c,norm_grad_c);
                double damping = CancerParameters::Instance()->GetDampingConstantNormal();
                
                // Fc = force_magnitude*(1,0), Fspring=0 => velocity = damping*force_magnitude*(1,0)
                TS_ASSERT_DELTA(velocities_on_each_node[i][0], force_magnitude/damping, 1e-4);
                TS_ASSERT_DELTA(velocities_on_each_node[i][1], 0.0, 1e-4);
            }
        }
    }
    
    void TestArchiving() throw (Exception)
    {   
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "chemotaxis_spring_system.arch";

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
        
            MeshBasedTissue<2> tissue(mesh,cells);
                    
            MeinekeSpringSystemWithChemotaxis<2> chemotaxis_spring_system(tissue);
         
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            MeinekeSpringSystemWithChemotaxis<2> * const p_chemotaxis_spring_system = &chemotaxis_spring_system;  
            
            p_chemotaxis_spring_system->UseCutoffPoint(1.1);
            p_chemotaxis_spring_system->SetAreaBasedViscosity(true);
            p_chemotaxis_spring_system->SetMutantSprings(true,0.2,0.3);
            p_chemotaxis_spring_system->SetBCatSprings(true);
            p_chemotaxis_spring_system->SetNecroticSprings(true);
            
            output_arch << p_chemotaxis_spring_system;
        }
       
        {
            MeshBasedTissue<2>::meshPathname = "mesh/test/data/square_2_elements";
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            MeinekeSpringSystemWithChemotaxis<2>* p_chemotaxis_spring_system;
            
            // Restore from the archive
            input_arch >> p_chemotaxis_spring_system;
            
            // Test the member data
            TS_ASSERT_EQUALS(p_chemotaxis_spring_system->mUseCutoffPoint,true);
            TS_ASSERT_DELTA(p_chemotaxis_spring_system->mCutoffPoint,1.1,1e-12);            
            TS_ASSERT_EQUALS(p_chemotaxis_spring_system->mUseEdgeBasedSpringConstant, false);
            TS_ASSERT_EQUALS(p_chemotaxis_spring_system->mUseAreaBasedViscosity, true);
            TS_ASSERT_EQUALS(p_chemotaxis_spring_system->mUseMutantSprings, true);
            TS_ASSERT_DELTA(p_chemotaxis_spring_system->mMutantMutantMultiplier, 0.2, 1e-12);
            TS_ASSERT_DELTA(p_chemotaxis_spring_system->mNormalMutantMultiplier, 0.3, 1e-12);
            TS_ASSERT_EQUALS(p_chemotaxis_spring_system->mUseBCatSprings, true);
            TS_ASSERT_EQUALS(p_chemotaxis_spring_system->mUseNecroticSprings, true);
            
            delete p_chemotaxis_spring_system->mpTissue;
            delete p_chemotaxis_spring_system;
        }
    }
};


#endif /*TESTMEINEKESPRINGSYSTEMWITHCHEMOTAXIS_HPP_*/
