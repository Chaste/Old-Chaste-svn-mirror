#ifndef TESTCELLWISEDATAGRADIENT_HPP_
#define TESTCELLWISEDATAGRADIENT_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include <cmath>
#include <vector>
#include "Tissue.cpp"
#include "CellwiseData.cpp"
#include "CellwiseDataGradient.hpp"
#include "CellsGenerator.hpp"
#include "CommonCancerTestSetup.hpp"


class TestCellwiseDataGradient : public AbstractCancerTestSuite
{    

public:
    void TestCellwiseDataGradientVerySmallMesh() throw(Exception)
    {        
        // create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // create a tissue
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateBasic(cells, mesh);
        Tissue<2> tissue(mesh,cells);

        // set up data: C(x,y) = x^2
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(mesh.GetNumNodes(), 1);
        p_data->SetTissue(tissue);     
    
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            p_data->SetValue(x*x, mesh.GetNode(i));
        }

        CellwiseDataGradient<2> gradient;
        gradient.SetupGradients();
        
        // With the algorithm being used, the numerical gradient is (1,0)
        // for each of the nodes 
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA( gradient.rGetGradient(i)(0), 1.0, 1e-9);
            TS_ASSERT_DELTA( gradient.rGetGradient(i)(1), 0.0, 1e-9);
        }
        
        CellwiseData<2>::Destroy();
    }


    void TestCellwiseDataGradientFineMesh() throw(Exception)
    {        
        // create a mesh: [0,2]x[0,2]
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4096_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // create a tissue
        std::vector<TissueCell> cells;
        CellsGenerator<2>::GenerateBasic(cells, mesh);
        Tissue<2> tissue(mesh,cells);

        //////////////////////////////////
        // C(x,y) = const
        //////////////////////////////////
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(mesh.GetNumNodes(), 1);
        p_data->SetTissue(tissue);     
    
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, mesh.GetNode(i));
        }

        CellwiseDataGradient<2> gradient;
        gradient.SetupGradients();
        
        // check gradient 
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA( gradient.rGetGradient(i)(0), 0.0, 1e-9);
            TS_ASSERT_DELTA( gradient.rGetGradient(i)(1), 0.0, 1e-9);
        }


        //////////////////////////////////
        // C(x,y) = x-y
        //////////////////////////////////
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            p_data->SetValue(x-y, mesh.GetNode(i));
        }

        // check gradient 
        gradient.SetupGradients();
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA( gradient.rGetGradient(i)(0),  1.0, 1e-9);
            TS_ASSERT_DELTA( gradient.rGetGradient(i)(1), -1.0, 1e-9);
        }

        //////////////////////////////////
        // C(x,y) = x^2 - y^2
        //////////////////////////////////
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            p_data->SetValue(x*x - y*y, mesh.GetNode(i));
        }

        // check gradient - here there is some numerical error
        gradient.SetupGradients();
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];

            double tol = 0.3;
            if(x==0 || x==2 || y==0 || y==2) //ie on boundary
            {
                tol = 0.6;
            }
                
            TS_ASSERT_DELTA( gradient.rGetGradient(i)(0),  2*x, tol);
            TS_ASSERT_DELTA( gradient.rGetGradient(i)(1), -2*y, tol);
        }

        CellwiseData<2>::Destroy();
    }
    
//    void TestCellwiseDataGradientWithGhostNodes() throw(Exception)
//    {        
//        // create a mesh: [0,2]x[0,2]
//        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4096_elements");
//        ConformingTetrahedralMesh<2,2> mesh;
//        mesh.ConstructFromMeshReader(mesh_reader);
//        
//        // create a tissue
//        std::vector<TissueCell> cells;
//        CellsGenerator<2>::GenerateBasic(cells, mesh);
//        
//        Tissue<2> tissue(mesh,cells);
//        // set boundary nodes to be ghosts
//        std::set< unsigned > ghost_node_indices;
//        for (Tissue<2>::Iterator iter=tissue.Begin();
//             iter != tissue.End();
//             ++iter)
//        {
//            if (iter.GetNode()->IsBoundaryNode())
//            {
//                ghost_node_indices.insert( iter.GetNode()->GetIndex() );
//                std::cout << iter.GetNode()->GetIndex() << " ";
//            }
//        }
//        tissue.SetGhostNodes(ghost_node_indices);
//        
//        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
//        p_data->SetNumNodesAndVars(mesh.GetNumNodes(), 1);
//        p_data->SetTissue(tissue);     
//    
//
//        //////////////////////////////////
//        // C(x,y) = x^2 - y^2
//        //////////////////////////////////
//        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
//        {
//            double x = mesh.GetNode(i)->rGetLocation()[0];
//            double y = mesh.GetNode(i)->rGetLocation()[1];
//            if (mesh.GetNode(i)->IsBoundaryNode())
//            {
//                p_data->SetValue(DBL_MAX, mesh.GetNode(i));
//            }
//            else
//            {
//                p_data->SetValue(x*x - y*y, mesh.GetNode(i));
//            }
//        }
//
//        // check gradient - here there is some numerical error
//        CellwiseDataGradient<2> gradient;
//        gradient.SetupGradients();
//        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
//        {
//            double x = mesh.GetNode(i)->rGetLocation()[0];
//            double y = mesh.GetNode(i)->rGetLocation()[1];
//
//            double tol = 0.6;
//            if ( !mesh.GetNode(i)->IsBoundaryNode() )//ie not ghost
//            {   
//                TS_ASSERT_DELTA( gradient.rGetGradient(i)(0),  2*x, tol);
//                TS_ASSERT_DELTA( gradient.rGetGradient(i)(1), -2*y, tol);
//            }
//        }
//
//        CellwiseData<2>::Destroy();
//    }

};        

#endif /*TESTCELLWISEDATAGRADIENT_HPP_*/
