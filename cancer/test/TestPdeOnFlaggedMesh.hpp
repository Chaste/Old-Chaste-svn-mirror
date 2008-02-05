#ifndef TESTPDEONFLAGGEDMESH_HPP_
#define TESTPDEONFLAGGEDMESH_HPP_

#include <cxxtest/TestSuite.h>

#include "HoneycombMeshGenerator.hpp"
#include "ParabolicFlaggedMeshAssembler.hpp"
#include "EllipticFlaggedMeshAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "FlaggedMeshBoundaryConditionsContainer.hpp"
#include "SimpleDataWriter.hpp"
#include "PetscTools.hpp"



/*
 *  A simple oxygen diffusion pde, with parameters chosen for the scaling in tissue simulations
 *  
 *  oxygen diffusion: u_t = D (u_xx+u_yy) - A
 * 
 *  Non-dim =>  (1/T) u_t = (D/L^2) (u_xx+u_yy) - A
 *          =>        u_t = (DT/L^2)(u_xx+u_yy) - AT
 * 
 *  Then choose A such AT~(DT/L^2), and divide by (DT/L^2), so
 *  SourceTerm=-1, DiffusionTerm=I, DuDtCoefficient=L^2/DT = (1e-3)^2/1e-5.1e2 = 1e-3
 */ 
class SimpleDiffusionPde : public AbstractLinearParabolicPde<2>
{

public:
    double ComputeLinearSourceTerm(const ChastePoint<2>& )
    {
        return -1;
    }
    
    c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& )
    {
        return identity_matrix<double>(2);
    }
    
    double ComputeDuDtCoefficientFunction(const ChastePoint<2>& )
    {
        return 1;
    }
    
};

class SimpleEllipticPde : public AbstractLinearEllipticPde<2>
{

public:
    double ComputeConstantInUSourceTerm(const ChastePoint<2>& )
    {
        return -1;
    }
    
    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>& )
    {
        return 0.0;
    }    

    c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& )
    {
        return identity_matrix<double>(2);
    }
        
};


class TestPdeOnFlaggedMesh: public CxxTest::TestSuite
{

public:    
    void TestPdeOnHoneycombMesh() throw (Exception)
    {
        if (!PetscTools::IsSequential())
        {
            TS_TRACE("This test does not pass in parallel yet. HoneycombMeshGenerator writes and reads files which produce a race-condition");
            return;
        }
        int num_cells_depth = 11;
        int num_cells_width = 6;
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();

        //for testing on a fine mesh
        //ConformingTetrahedralMesh<2,2>* p_mesh=new ConformingTetrahedralMesh<2,2>(); 
        //p_mesh->ConstructRectangularMesh(50,50);
        //p_mesh->Scale(1.0/50, 1.0/50); 
        //std::set<unsigned> ghost_node_indices;
        
        p_mesh->FlagElementsNotContainingNodes(ghost_node_indices);
        
        unsigned counter=0; 
        ConformingTetrahedralMesh<2,2>::ElementIterator
            iter = p_mesh->GetElementIteratorBegin();
       
        while (iter != p_mesh->GetElementIteratorEnd())
        {
            Element<2,2>& element = **iter;
            if(element.IsFlagged())
            {
                counter++;
            }
            ++iter;
        }
        TS_ASSERT_EQUALS(counter,100u);

        // Set up boundary conditions
        FlaggedMeshBoundaryConditionsContainer<2,1> flagged_bcc(*p_mesh, 1.0);
               
        
        ////////////////////////////////////////////////////////
        // solve with elliptic equation
        ////////////////////////////////////////////////////////

        // Instantiate PDE object
        SimpleEllipticPde elliptic_pde;
        
        // Assembler for fine mesh flagged region. use bbc from before
        EllipticFlaggedMeshAssembler<2> elliptic_assembler(p_mesh, &elliptic_pde, &flagged_bcc);
        
        Vec result_elliptic_restricted = elliptic_assembler.Solve();
        ReplicatableVector result_elliptic_repl(result_elliptic_restricted);

        std::map<unsigned, unsigned>& map2 = p_mesh->rGetSmasrmMap();
        std::map<unsigned, unsigned>::iterator map_iter2 = map2.begin();

        std::vector<double> x2;
        std::vector<double> y2;
        std::vector<double> u2;

        while (map_iter2!=map2.end())
        {
            unsigned node_index = map_iter2->first;
            unsigned smasrm_index = map_iter2->second;

            // check the node isn't a ghost node
            TS_ASSERT_EQUALS( ghost_node_indices.find(node_index), ghost_node_indices.end() );
            
            x2.push_back(p_mesh->GetNode(node_index)->rGetLocation()[0]);
            y2.push_back(p_mesh->GetNode(node_index)->rGetLocation()[1]);
            u2.push_back(result_elliptic_repl[smasrm_index]);

            map_iter2++;
        }
        

        std::vector<std::vector<double> > data2;
        data2.push_back(x2);
        data2.push_back(y2);
        data2.push_back(u2);
        
        //SimpleDataWriter writer2("temp", "ellip", data2, false);

        //TS_ASSERT_EQUALS(result_elliptic_repl.size(),result_diffusion_repl.size());
        //for (unsigned i=0; i< result_diffusion_repl.size(); i++)
        //{
        //    TS_ASSERT_DELTA(result_elliptic_repl[i],result_diffusion_repl[i],1e-2);
        //}
        
        VecDestroy(result_elliptic_restricted);
    }
};

#endif /*TESTPDEONFLAGGEDMESH_HPP_*/
