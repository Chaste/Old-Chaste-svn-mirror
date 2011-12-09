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


#ifndef TESTCARDIACELECTROMECHANICSPROBLEMLONG_HPP_
#define TESTCARDIACELECTROMECHANICSPROBLEMLONG_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechProbRegularGeom.hpp"
#include "LuoRudy1991.hpp"
#include "NonlinearElasticityTools.hpp"
#include "NobleVargheseKohlNoble1998WithSac.hpp"

class TestCardiacElectroMechanicsProblemLong : public CxxTest::TestSuite
{
public:
    void Test2dHardcodedResult() throw(Exception)
    {
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-1000*1000);

        // run to 125 ms - about where the width is at its minimum (see figures
        // in "A numerical method for cardiac mechano–electric simulations", Annals of Biomed. Imaging

        HeartConfig::Instance()->SetSimulationDuration(125.0);

        CardiacElectroMechProbRegularGeom<2> problem(INCOMPRESSIBLE,
                                                     1.0,  /* width */
                                                     5,    /* mech mesh size */
                                                     60,   /* elec elem each dir */
                                                     &cell_factory,
                                                     NHS,
                                                     1.0,  /* mechanics solve timestep */
                                                     1.0,  /* contraction model ode dt */
                                                     "TestCardiacEmNhs2dLong");

        problem.SetNoElectricsOutput();
        problem.Solve();

        // test by checking the length of the tissue against hardcoded value
        std::vector<c_vector<double,2> >& r_deformed_position = problem.rGetDeformedPosition();
        TS_ASSERT_DELTA(r_deformed_position[5](0), 0.8257, 1e-3);

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }

    void Test2dVariableFibres() throw(Exception)
    {
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-1000*1000);

        TetrahedralMesh<2,2> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(1.0/96.0/*stepsize*/, 1.0/*length*/, 1.0/*width*/, 1.0/*depth*/);

        QuadraticMesh<2> mechanics_mesh;
        mechanics_mesh.ConstructRegularSlabMesh(0.2, 1.0, 1.0, 1.0 /*as above with a different stepsize*/);

        std::vector<unsigned> fixed_nodes
            = NonlinearElasticityTools<2>::GetNodesByComponentValue(mechanics_mesh, 0, 0.0); // all the X=0.0 nodes

        ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(NHS,1.0);
        problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetMechanicsSolveTimestep(1.0);
        problem_defn.SetVariableFibreSheetDirectionsFile("heart/test/data/fibre_tests/5by5mesh_curving_fibres.ortho", false);

        HeartConfig::Instance()->SetSimulationDuration(125.0);

        CardiacElectroMechanicsProblem<2> problem(INCOMPRESSIBLE,
                                                  &electrics_mesh,
                                                  &mechanics_mesh,
                                                  &cell_factory,
                                                  &problem_defn,
                                                  "TestCardiacEmVaryingFibres");


        // fibres going from (1,0) at X=0 to (1,1)-direction at X=1
        /* the fibres file was created with the code (inside a class that owns a mesh)
        for(unsigned elem_index=0; elem_index<this->mpQuadMesh->GetNumElements(); elem_index++)
        {
            double X = this->mpQuadMesh->GetElement(elem_index)->CalculateCentroid()[0];
            double theta = M_PI*X/4;
            std::cout << cos(theta) << " " << sin(theta) << " " << -sin(theta) << " " << cos(theta) << "\n" << std::flush;
        }
        assert(0);
        */

        // problem.SetNoElectricsOutput();
        problem.Solve();

        // test by checking the length of the tissue against hardcoded value
        std::vector<c_vector<double,2> >& r_deformed_position = problem.rGetDeformedPosition();
        // visualised, looks good - contracts in X-direction near the fixed surface,
        // but on the other side the fibres are in the (1,1) direction, so contraction
        // pulls the tissue downward a bit
        TS_ASSERT_DELTA(r_deformed_position[5](0), 0.9019, 2e-3); 
        //IntelProduction differs by about 1.6e-3...

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }



    void Test3d() throw(Exception)
    {
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory(-1000*1000);

        // set up two meshes of 1mm by 1mm by 1mm
        TetrahedralMesh<3,3> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.01, 0.1, 0.1, 0.1);

        QuadraticMesh<3> mechanics_mesh(0.1, 0.1, 0.1, 0.1);

        // fix the nodes on x=0
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<3>::GetNodesByComponentValue(mechanics_mesh,0,0);


        HeartConfig::Instance()->SetSimulationDuration(50.0);

        ElectroMechanicsProblemDefinition<3> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(NHS,1.0);
        problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetMechanicsSolveTimestep(1.0);

        CardiacElectroMechanicsProblem<3> problem(INCOMPRESSIBLE,
                                                  &electrics_mesh,
                                                  &mechanics_mesh,
                                                  &cell_factory,
                                                  &problem_defn,
                                                  "TestCardiacElectroMech3d");

        problem.SetNoElectricsOutput();
        problem.Solve();

        // test by checking the length of the tissue against hardcoded value
        std::vector<c_vector<double,3> >& r_deformed_position = problem.rGetDeformedPosition();
        TS_ASSERT_DELTA(r_deformed_position[1](0), 0.0879, 1e-3);

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }
    
    void TestTwistingCube() throw(Exception)
    {
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory(-1000*1000);

        // set up two meshes of 1mm by 1mm by 1mm
        TetrahedralMesh<3,3> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.01, 0.1, 0.1, 0.1);
    
        QuadraticMesh<3> mechanics_mesh(0.02, 0.1, 0.1, 0.1);

        // fix the nodes on Z=0
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<3>::GetNodesByComponentValue(mechanics_mesh,2,0.0);

        HeartConfig::Instance()->SetSimulationDuration(50.0);

        ElectroMechanicsProblemDefinition<3> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(KERCHOFFS2003,1.0);
        problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetMechanicsSolveTimestep(1.0);
        problem_defn.SetVariableFibreSheetDirectionsFile("heart/test/data/fibre_tests/5by5by5_fibres_by_quadpt.orthoquad", true);

        CardiacElectroMechanicsProblem<3> problem(INCOMPRESSIBLE,
                                                  &electrics_mesh,
                                                  &mechanics_mesh,
                                                  &cell_factory,
                                                  &problem_defn,
                                                  "TestCardiacElectroMech3dTwistingCube");


/////// Use the following to set up the fibres file        
////        GaussianQuadratureRule<3> quad_rule(3);
////        QuadraturePointsGroup<3> quad_points(mechanics_mesh, quad_rule);
////        std::cout << quad_points.Size() << "\n"; 
////        for(unsigned i=0; i<quad_points.Size(); i++)
////        {
////            ////std::cout << quad_points.Get(i)(0) << " " << quad_points.Get(i)(1) << " " << quad_points.Get(i)(2) << " ";
////            double x = quad_points.Get(i)(0);
////            double theta = M_PI/3 - 10*x*2*M_PI/3; // 60 degrees when x=0, -60 when x=0.1; 
////            std::cout <<  "0 " << cos(theta)  << " " << sin(theta) 
////                      << " 0 " << -sin(theta) << " " << cos(theta)
////                      << " 1 0 0\n";  
////        }
////        return;


        problem.Solve();

        // verified that it twists by visualising, some hardcoded values here..

        std::vector<c_vector<double,3> >& r_deformed_position = problem.rGetDeformedPosition();
        TS_ASSERT_DELTA(r_deformed_position[6*6*5](0),  0.0116, 1e-3);
        TS_ASSERT_DELTA(r_deformed_position[6*6*5](1), -0.0141, 1e-3);
        TS_ASSERT_DELTA(r_deformed_position[6*6*5](2),  0.1007, 1e-3);

        TS_ASSERT_DELTA(r_deformed_position[6*6*6-1](0), 0.0872, 1e-3);
        TS_ASSERT_DELTA(r_deformed_position[6*6*6-1](1), 0.1138, 1e-3);
        TS_ASSERT_DELTA(r_deformed_position[6*6*6-1](2), 0.1015, 1e-3);

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();

    }


    // longer running, finer-mesh version of TestWithCompressibleApproach() in TestCardiacElectroMechanicsProblem.hpp
    void TestWithCompressibleApproachLong() throw(Exception)
    {
        EXIT_IF_PARALLEL; // #1913 currently, the compressible preconditioner is ICC, which is only supported in sequential

        HeartEventHandler::Disable();

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-1000*1000);

        HeartConfig::Instance()->SetSimulationDuration(50.0);

        CardiacElectroMechProbRegularGeom<2> problem(COMPRESSIBLE,
                                                     0.05, /* width (cm) */
                                                     20,   /* mech mesh size*/
                                                     5,    /* elec elem each dir */
                                                     &cell_factory,
                                                     KERCHOFFS2003,
                                                     1.0,   /* mechanics solve timestep */
                                                     0.01,  /* Kerchoffs ode timestep */
                                                     "TestCompressibleWithKerchoffsLong");

        problem.Solve();

        // Mainly just testing no errors when Solve was called.
        // The results of this test can be visually compared with the results of the
        // equivalent incompressible simulation in TestWithKerchoffs.

        TS_ASSERT_DELTA(problem.rGetDeformedPosition()[20](0), 0.0429, 0.0002);
        TS_ASSERT_DELTA(problem.rGetDeformedPosition()[20](1),-0.0038, 0.0002);
    }

    // runs 5 different 3d tests with fibres read to be in various directions, and
    // checks contraction occurs in the right directions (and bulging occurs in the
    // other directions)
    void TestFibreRead() throw(Exception)
    {
        EXIT_IF_PARALLEL; ///\todo #1913

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML,3> cell_factory(-5000*1000);

        double tissue_initial_size = 0.05;
        TetrahedralMesh<3,3> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.01/*stepsize*/, tissue_initial_size/*length*/, tissue_initial_size/*width*/, tissue_initial_size);

        QuadraticMesh<3> mechanics_mesh;
        mechanics_mesh.ConstructRegularSlabMesh(0.05, tissue_initial_size, tissue_initial_size, tissue_initial_size /*as above with a different stepsize*/);

        std::vector<unsigned> fixed_nodes
            = NonlinearElasticityTools<3>::GetNodesByComponentValue(mechanics_mesh, 2, 0.0);


        HeartConfig::Instance()->SetSimulationDuration(10.0);

        ElectroMechanicsProblemDefinition<3> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(KERCHOFFS2003,1.0);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetMechanicsSolveTimestep(1.0);
        problem_defn.SetUseDefaultCardiacMaterialLaw(COMPRESSIBLE);

        // This is how the fibres files that are used in the simulations are created.
        // alongX.ortho is:  fibres in X axis, sheet in XY
        //    -- same as default directions: should give shortening in X direction
        //       and identical results to default fibre setting
        // alongY1.ortho is: fibres in Y axis, sheet in YX
        // alongY2.ortho is: fibres in Y axis, sheet in YZ
        //    -- as material law is transversely isotropic, these two should
        //       give identical results, shortening in Y-direction
        // alongZ.ortho is:  fibres in Z axis, (sheet in XZ)
        //    -- should shorten in Z-direction
        OutputFileHandler handler("FibreFiles");
        out_stream p_X_file  = handler.OpenOutputFile("alongX.ortho");
        out_stream p_Y1_file = handler.OpenOutputFile("alongY1.ortho");
        out_stream p_Y2_file = handler.OpenOutputFile("alongY2.ortho");
        out_stream p_Z_file  = handler.OpenOutputFile("alongZ.ortho");
        *p_X_file  << mechanics_mesh.GetNumElements() << "\n";
        *p_Y1_file << mechanics_mesh.GetNumElements() << "\n";
        *p_Y2_file << mechanics_mesh.GetNumElements() << "\n";
        *p_Z_file  << mechanics_mesh.GetNumElements() << "\n";
        for(unsigned i=0; i<mechanics_mesh.GetNumElements(); i++)
        {
            //double X = mechanics_mesh.GetElement(i)->CalculateCentroid()(0);
            *p_X_file  << "1 0 0 0 1 0 0 0 1\n";
            *p_Y1_file << "0 1 0 1 0 0 0 0 1\n";
            *p_Y2_file << "0 1 0 0 0 1 1 0 0\n";
            *p_Z_file  << "0 0 1 1 0 0 0 1 0\n";
        }
        p_X_file->close();
        p_Y1_file->close();
        p_Y2_file->close();
        p_Z_file->close();


        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();

        //////////////////////////////////////////////////////////////////
        // Solve with no fibres read.
        //////////////////////////////////////////////////////////////////
        std::vector<c_vector<double,3> > r_deformed_position_no_fibres;
        {
            HeartConfig::Instance()->SetSimulationDuration(20.0);

            CardiacElectroMechanicsProblem<3> problem(COMPRESSIBLE,
                                                      &electrics_mesh,
                                                      &mechanics_mesh,
                                                      &cell_factory,
                                                      &problem_defn,
                                                      "TestCardiacEmFibreRead");


            problem.Solve();
            r_deformed_position_no_fibres = problem.rGetDeformedPosition();
        }

        //////////////////////////////////////////////////////////////////
        // Solve with fibres read: fibres in X-direction
        //////////////////////////////////////////////////////////////////
        std::vector<c_vector<double,3> > r_deformed_position_fibres_alongX;
        {
            HeartConfig::Instance()->SetSimulationDuration(20.0);
            problem_defn.SetVariableFibreSheetDirectionsFile("heart/test/data/fibre_tests/alongX.ortho", false);

            CardiacElectroMechanicsProblem<3> problem(COMPRESSIBLE,
                                                      &electrics_mesh,
                                                      &mechanics_mesh,
                                                      &cell_factory,
                                                      &problem_defn,
                                                      "TestCardiacEmFibreRead");


            problem.Solve();
            r_deformed_position_fibres_alongX = problem.rGetDeformedPosition();
        }

        // test the two results are identical
        for(unsigned i=0; i<r_deformed_position_no_fibres.size(); i++)
        {
            for(unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(r_deformed_position_no_fibres[i](j), r_deformed_position_fibres_alongX[i](j), 1e-8);
            }
        }

        // check node 7 starts is the far corner
        assert(fabs(mechanics_mesh.GetNode(7)->rGetLocation()[0] - tissue_initial_size)<1e-8);
        assert(fabs(mechanics_mesh.GetNode(7)->rGetLocation()[1] - tissue_initial_size)<1e-8);
        assert(fabs(mechanics_mesh.GetNode(7)->rGetLocation()[2] - tissue_initial_size)<1e-8);

        // Test that contraction occurred in the X-direction
        TS_ASSERT_LESS_THAN(r_deformed_position_fibres_alongX[7](0), tissue_initial_size);
        TS_ASSERT_LESS_THAN(tissue_initial_size, r_deformed_position_fibres_alongX[7](1));
        TS_ASSERT_LESS_THAN(tissue_initial_size, r_deformed_position_fibres_alongX[7](2));

        // hardcoded test to check nothing changes
        TS_ASSERT_DELTA(r_deformed_position_fibres_alongX[7](0), 0.0484, 1e-4);
        TS_ASSERT_DELTA(r_deformed_position_fibres_alongX[7](1), 0.0508, 1e-4);

        ////////////////////////////////////////////////////////////////////
        // Solve with fibres read: fibres in Y-direction, sheet in YX plane
        ////////////////////////////////////////////////////////////////////
        std::vector<c_vector<double,3> > r_deformed_position_fibres_alongY1;
        {
            HeartConfig::Instance()->SetSimulationDuration(20.0);
            problem_defn.SetVariableFibreSheetDirectionsFile("heart/test/data/fibre_tests/alongY1.ortho", false);

            CardiacElectroMechanicsProblem<3> problem(COMPRESSIBLE,
                                                      &electrics_mesh,
                                                      &mechanics_mesh,
                                                      &cell_factory,
                                                      &problem_defn,
                                                      "TestCardiacEmFibreRead");


            problem.Solve();
            r_deformed_position_fibres_alongY1 = problem.rGetDeformedPosition();
        }


        ////////////////////////////////////////////////////////////////////
        // Solve with fibres read: fibres in Y-direction, sheet in YZ plane
        ////////////////////////////////////////////////////////////////////
        std::vector<c_vector<double,3> > r_deformed_position_fibres_alongY2;
        {
            HeartConfig::Instance()->SetSimulationDuration(20.0);
            problem_defn.SetVariableFibreSheetDirectionsFile("heart/test/data/fibre_tests/alongY2.ortho", false);

            CardiacElectroMechanicsProblem<3> problem(COMPRESSIBLE,
                                                      &electrics_mesh,
                                                      &mechanics_mesh,
                                                      &cell_factory,
                                                      &problem_defn,
                                                      "TestCardiacEmFibreRead");


            problem.Solve();
            r_deformed_position_fibres_alongY2 = problem.rGetDeformedPosition();
        }


        // test the two results are identical
        for(unsigned i=0; i<r_deformed_position_no_fibres.size(); i++)
        {
            for(unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(r_deformed_position_fibres_alongY1[i](j), r_deformed_position_fibres_alongY2[i](j), 1e-8);
            }
        }

        // Test that contraction occurred in the Y-direction
        TS_ASSERT_LESS_THAN(tissue_initial_size, r_deformed_position_fibres_alongY1[7](0));
        TS_ASSERT_LESS_THAN(r_deformed_position_fibres_alongY1[7](1), tissue_initial_size);
        TS_ASSERT_LESS_THAN(tissue_initial_size, r_deformed_position_fibres_alongY1[7](2));

        // hardcoded test to check nothing changes
        TS_ASSERT_DELTA(r_deformed_position_fibres_alongY1[7](1), 0.0486, 1e-4);
        TS_ASSERT_DELTA(r_deformed_position_fibres_alongY1[7](0), 0.0505, 1e-4);


        //////////////////////////////////////////////////////////////////
        // Solve with fibres read: fibres in Z-direction
        //////////////////////////////////////////////////////////////////
        std::vector<c_vector<double,3> > r_deformed_position_fibres_alongZ;
        {
            HeartConfig::Instance()->SetSimulationDuration(20.0);
            problem_defn.SetVariableFibreSheetDirectionsFile("heart/test/data/fibre_tests/alongZ.ortho", false);

            CardiacElectroMechanicsProblem<3> problem(COMPRESSIBLE,
                                                      &electrics_mesh,
                                                      &mechanics_mesh,
                                                      &cell_factory,
                                                      &problem_defn,
                                                      "TestCardiacEmFibreRead");


            problem.Solve();
            r_deformed_position_fibres_alongZ = problem.rGetDeformedPosition();
        }

        // Test that contraction occurred in the X-direction
        TS_ASSERT_LESS_THAN(tissue_initial_size, r_deformed_position_fibres_alongZ[7](0));
        TS_ASSERT_LESS_THAN(tissue_initial_size, r_deformed_position_fibres_alongZ[7](1));
        TS_ASSERT_LESS_THAN(r_deformed_position_fibres_alongZ[7](2), tissue_initial_size);

        // hardcoded test to check nothing changes
        TS_ASSERT_DELTA(r_deformed_position_fibres_alongZ[7](2), 0.0456, 1e-4);
        TS_ASSERT_DELTA(r_deformed_position_fibres_alongZ[7](0), 0.0503, 1e-4);
    }


//    void Test3dWithNoble98SacAndImpact() throw(Exception)
//    {
//        // zero stimuli
//        PlaneStimulusCellFactory<CML_noble_varghese_kohl_noble_1998_basic_with_sac, 3> cell_factory(0);
//
//        // set up two meshes of 1cm by 1cm by 1cm
//        TetrahedralMesh<3,3> electrics_mesh;
//        unsigned num_elem = 10;
//        electrics_mesh.ConstructCuboid(num_elem,num_elem,num_elem);
//        electrics_mesh.Scale(1.0/num_elem, 1.0/num_elem, 1.0/num_elem);
//
//        QuadraticMesh<3> mechanics_mesh(0.2, 1.0, 1.0, 1.0);
//
//        // fix the nodes on x=0
//        std::vector<unsigned> fixed_nodes
//          = NonlinearElasticityTools<3>::GetNodesByComponentValue(mechanics_mesh,0,0);
//
//        std::vector<BoundaryElement<2,3>*> impact_region;        
//        for (TetrahedralMesh<3,3>::BoundaryElementIterator iter
//              = mechanics_mesh.GetBoundaryElementIteratorBegin();
//             iter != mechanics_mesh.GetBoundaryElementIteratorEnd();
//             ++iter)
//        {
//            c_vector<double,3> centroid =(*iter)->CalculateCentroid(); 
//            if (    (fabs(centroid[1])<1e-4) 
//                 && (centroid[0] < 0.05)
//                 && (centroid[2] < 0.05) )
//            {
//                BoundaryElement<2,3>* p_element = *iter;
//                impact_region.push_back(p_element);
//            }
//        }
//        assert(impact_region.size() > 0);
//
//        CardiacElectroMechanicsProblem<3> problem(INCOMPRESSIBLE,
//                                                  KERCHOFFS2003,
//                                                  &electrics_mesh,
//                                                  &mechanics_mesh,
//                                                  fixed_nodes,
//                                                  &cell_factory,
//                                                  10,   /* end time */
//                                                  0.01, /* electrics timestep (ms) */
//                                                  1.0,  /* mechanics solve timestep */
//                                                  1.0,  /* contraction model ode dt */
//                                                  "TestCardiacElectroMech3dImpact");
//
//        problem.SetImpactRegion(impact_region);
//
//        problem.Solve();
//
//        // test by checking the length of the tissue against hardcoded value
//        std::vector<c_vector<double,3> >& r_deformed_position = problem.rGetDeformedPosition();
//        TS_ASSERT_DELTA(r_deformed_position[1](0), 0.0879, 1e-3);
//
//        MechanicsEventHandler::Headings();
//        MechanicsEventHandler::Report();
//    }
};
#endif /*TESTCARDIACELECTROMECHANICSPROBLEMLONG_HPP_*/
