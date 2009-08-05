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

#ifndef TESTCARDIACDOMAINFASTSLOWPROBLEM_HPP_
#define TESTCARDIACDOMAINFASTSLOWPROBLEM_HPP_

#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "SimpleStimulus.hpp"
#include "FastSlowBackwardEulerNoble98.hpp"
#include "BackwardEulerNobleVargheseKohlNoble1998.hpp"
#include "MonodomainFastSlowProblem.hpp"
#include "BidomainFastSlowProblem.hpp"
#include "NobleVargheseKohlNoble1998.hpp"
#include "NobleVargheseKohlNoble1998Optimised.hpp"
#include "TrianglesMeshReader.hpp"

// simple cell factory that creates fast-slow cells.
template <class CELL>
class MyCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
    bool mForSquare; // whether this factory is for the square or disk

public:

    MyCellFactory(bool forSquare=true) :
        AbstractCardiacCellFactory<2>(),
        mpStimulus(new SimpleStimulus(-10000*250.0, 0.2)),        //Where did this come from?
        //Will it work with LR91 and N98?
        mForSquare(forSquare)
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        if(mForSquare)
        {
            double x = GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];

            // stimulus is zero unless x=0
            boost::shared_ptr<AbstractStimulusFunction> p_stim;
            if(fabs(x)<1e-6)
            {
                p_stim = mpStimulus;
            }
            else
            {
                p_stim = mpZeroStimulus;
            }

            return new CELL(mpSolver, p_stim);
        }
        else
        {
            double y = GetMesh()->GetNode(nodeIndex)->rGetLocation()[1];

            // stimulus is zero unless at edge vertex
            boost::shared_ptr<AbstractStimulusFunction> p_stim;
            if(y < -0.09)
            {
                p_stim = mpStimulus;
            }
            else
            {
                p_stim = mpZeroStimulus;
            }

            return new CELL(mpSolver, p_stim);
        }
    }
};


class TestCardiacFastSlowProblem : public CxxTest::TestSuite
{
private:
    //Make sure that Mono/Bi are solving with the same parameters
    double mSlowOdeTimestep;
    MixedTetrahedralMesh<2,2> mMixedMesh;
    TetrahedralMesh<2,2> *mpFineMesh;
    unsigned mProbeIndex;

    void SetHeartConfigForTest()
    {
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");

        //Found this figure of absolute tolerance by running a standard Monodomain
        //simulation with default KSP parameters (reltol etc) and then finding an
        //absolute tolerance with the same CPU time (same number of KSP iterations)
        HeartConfig::Instance()->SetUseAbsoluteTolerance(2e-5);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.05, 0.05, 0.05);
        HeartConfig::Instance()->SetSimulationDuration(0.2); //ms
        mSlowOdeTimestep = 0.1;
    }


    void CheckProblem(std::vector<double> &voltageFastSlow,
                      std::vector<double> &voltageNormal)
    {
        //////////////////////////////////////////////////
        // compare
        //////////////////////////////////////////////////
        TS_ASSERT_EQUALS(voltageFastSlow.size(), voltageNormal.size() );

        bool some_voltage_greater_than_zero = false;
        for(unsigned i=0; i<voltageFastSlow.size(); i++)
        {
            TS_ASSERT_DELTA(voltageFastSlow[i], voltageNormal[i], 5.0);
            if(voltageFastSlow[i] > 0.0)
            {
                some_voltage_greater_than_zero = true;
            }
        }
        TS_ASSERT(!some_voltage_greater_than_zero);
    }

    void RunAndCheckMonodomain(std::vector<double> &voltageFastSlow,
                                AbstractCardiacCellFactory<2>* pCellFactory,
                                std::string outputDirectory
                                )
    {

        MonodomainProblem<2> monodomain_prob(pCellFactory);
        HeartConfig::Instance()->SetOutputDirectory(outputDirectory);
        monodomain_prob.SetMesh(mpFineMesh);
        monodomain_prob.ConvertOutputToMeshalyzerFormat(true);
        monodomain_prob.Initialise();
        monodomain_prob.Solve();
        HeartEventHandler::Headings();
        HeartEventHandler::Report();
        HeartEventHandler::Reset();
        /* Read the vector at a typical downstream node */
        Hdf5DataReader results_reader=monodomain_prob.GetDataReader();
        std::vector<double> voltage_normal=results_reader.GetVariableOverTime("V", mProbeIndex);
        CheckProblem(voltageFastSlow, voltage_normal);
    }

    void RunAndCheckBidomain(std::vector<double> &voltageFastSlow,
                                AbstractCardiacCellFactory<2>* pCellFactory,
                                std::string outputDirectory)
    {
        std::vector<double> voltage_normal;
        HeartConfig::Instance()->SetOutputDirectory(outputDirectory);

        BidomainProblem<2> bidomain_prob( pCellFactory);
        bidomain_prob.SetMesh(mpFineMesh);

        bidomain_prob.ConvertOutputToMeshalyzerFormat(true);
        bidomain_prob.Initialise();
        bidomain_prob.Solve();
        HeartEventHandler::Headings();
        HeartEventHandler::Report();
        HeartEventHandler::Reset();
        /* Read the vector at a typical downstream node */
        Hdf5DataReader results_reader=bidomain_prob.GetDataReader();
        voltage_normal=results_reader.GetVariableOverTime("V", mProbeIndex);
        CheckProblem(voltageFastSlow, voltage_normal);
    }

public:
     void TestConstructAndThrowAwayMixedMesh3D()
    {
        double width = 1.; //cm
        double height = 1.; //cm
        double depth = 1.; //cm
        double coarse_ds = 1.; //cm
        double fine_ds = 0.5;

        unsigned num_coarse_elem_each_dir = unsigned(width/coarse_ds + 0.5);
        unsigned num_fine_elem_each_dir = unsigned(width/fine_ds + 0.5);

        TS_ASSERT_EQUALS(num_coarse_elem_each_dir, 1U);
        TS_ASSERT_EQUALS(num_fine_elem_each_dir, 2U);
        MixedTetrahedralMesh<3,3> mixed_mesh;
        mixed_mesh.ConstructCuboidMeshes(width, height, depth, num_coarse_elem_each_dir, num_fine_elem_each_dir);
        TetrahedralMesh<3,3> *p_fine_mesh=mixed_mesh.GetFineMesh();

        TS_ASSERT_EQUALS(mixed_mesh.GetNumNodes(), 8U);
        TS_ASSERT_EQUALS(p_fine_mesh->GetNumNodes(), 27U);

        unsigned node_ratio = p_fine_mesh->GetNumNodes() / mixed_mesh.GetNumNodes(); // Note -- rounded
        TS_ASSERT_EQUALS(node_ratio, 3U);//Rounded down from 3^3/2^3=3 3/8


     }
    void TestConstructAndKeepMixedMesh2D()
    {
        //This test is needed for later tests to operate
        double width = 0.5; //cm
        double height = 0.5; //cm
        double coarse_ds = 0.1; //cm
        double fine_ds = 0.01;

        unsigned num_coarse_elem_each_dir = unsigned(width/coarse_ds + 0.5);
        unsigned num_fine_elem_each_dir = unsigned(width/fine_ds + 0.5);

        TS_ASSERT_EQUALS(num_coarse_elem_each_dir, 5U);
        TS_ASSERT_EQUALS(num_fine_elem_each_dir, 50U);
        mMixedMesh.ConstructRectangularMeshes(width, height, num_coarse_elem_each_dir, num_fine_elem_each_dir);
        mpFineMesh=mMixedMesh.GetFineMesh();

        TS_ASSERT_EQUALS(mpFineMesh->GetNumNodes(), 2601U);
        TS_ASSERT_EQUALS(mMixedMesh.GetNumNodes(), 36U);


        //Find a place to probe on the middle of the right-hand side
        mProbeIndex = UINT_MAX;
        for (unsigned i=0; i<mpFineMesh->GetNumNodes(); i++){
            double x=mpFineMesh->GetNode(i)->rGetLocation()[0];
            if ( fabs (x - width) < 1e-6)
            {
                double y=mpFineMesh->GetNode(i)->rGetLocation()[1];
                if ( fabs(y - (height/2.0)) < 1e-6)
                {
                    mProbeIndex = i;
                    break;
                }
            }
        }
        TS_ASSERT(mProbeIndex != UINT_MAX);
    }

    // Run the Monodomain fast slow problem and compare solution with a normal problem
    //
    // This test only checks for consistency in both solutions.
    //
    void TestMonodomainFastSlowProblemAgainstNormal() throw (Exception)
    {
        SetHeartConfigForTest();

        //////////////////////////////////////////////////
        // solve a mixed mesh, fast/slow problem
        //////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputDirectory("MonodomainFastSlowTest");
        HeartConfig::Instance()->SetOutputFilenamePrefix("res");

        std::vector<double> voltage_fast_slow;
        {
            MyCellFactory<FastSlowBackwardEulerNoble98> cell_factory;

            MonodomainFastSlowProblem<2> monodomain_fast_slow_prob( &cell_factory, mMixedMesh, mSlowOdeTimestep );
            monodomain_fast_slow_prob.ConvertOutputToMeshalyzerFormat(true);
            monodomain_fast_slow_prob.Initialise();
            monodomain_fast_slow_prob.Solve();

            HeartEventHandler::Headings();
            HeartEventHandler::Report();
            HeartEventHandler::Reset();

            /* Read the vector at a typical downstream node */
            Hdf5DataReader results_reader=monodomain_fast_slow_prob.GetDataReader();

            voltage_fast_slow=results_reader.GetVariableOverTime("V", mProbeIndex);

            unsigned expected_size = (unsigned)(HeartConfig::Instance()->GetSimulationDuration()/HeartConfig::Instance()->GetPrintingTimeStep() ) + 1 ;
            TS_ASSERT_EQUALS(voltage_fast_slow.size(),  expected_size);
        }

        /////////////////////////////////////////////////////////
        // solve using normal monodomain problem - backward model
        /////////////////////////////////////////////////////////
        MyCellFactory<BackwardEulerNobleVargheseKohlNoble1998> cell_factory_backward;
        RunAndCheckMonodomain(voltage_fast_slow, &cell_factory_backward, "MonodomainBackwardToCompareWithFastSlowTest");

        ////In long tests we require a different time step for forward models
        //HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.05);

        //MyCellFactory<CML_noble_varghese_kohl_noble_1998_basic_pe_lut> cell_factory_lu;
        //RunAndCheckMonodomain(voltage_fast_slow, &cell_factory_lu, "MonodomainLookupToCompareWithFastSlowTest");

        //MyCellFactory<CML_noble_varghese_kohl_noble_1998_basic> cell_factory_fe;
        //RunAndCheckMonodomain(voltage_fast_slow, &cell_factory_fe, "MonodomainForwardEulerToCompareWithFastSlowTest");
    }



    // Run the Bidomain fast slow problem and compare solution with a normal problem
    //
    // This test only checks for consistency in both solutions. Bigger meshes and simulation times needed
    // to study performance
    //
    void TestBidomainFastSlowProblemAgainstNormal() throw (Exception)
    {
        SetHeartConfigForTest();

        //////////////////////////////////////////////////
        // solve a mixed mesh, fast/slow problem
        //////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputDirectory("BidomainFastSlowTest");
        HeartConfig::Instance()->SetOutputFilenamePrefix("res");

        std::vector<double> voltage_fast_slow;
        {
            MyCellFactory<FastSlowBackwardEulerNoble98> cell_factory;

            BidomainFastSlowProblem<2> bidomain_fast_slow_prob( &cell_factory, mMixedMesh, mSlowOdeTimestep );

            bidomain_fast_slow_prob.ConvertOutputToMeshalyzerFormat(true);
            bidomain_fast_slow_prob.Initialise();
            bidomain_fast_slow_prob.Solve();

            HeartEventHandler::Headings();
            HeartEventHandler::Report();
            HeartEventHandler::Reset();

           /* Read the vector at a typical downstream node */
            Hdf5DataReader results_reader=bidomain_fast_slow_prob.GetDataReader();
            voltage_fast_slow=results_reader.GetVariableOverTime("V", mProbeIndex);

            unsigned expected_size = (unsigned)(HeartConfig::Instance()->GetSimulationDuration()/HeartConfig::Instance()->GetPrintingTimeStep() ) + 1 ;
            TS_ASSERT_EQUALS(voltage_fast_slow.size(),  expected_size);
         }

        //////////////////////////////////////////////////////////
        // solve using normal bidomain problem - backward model
        /////////////////////////////////////////////////////////

        MyCellFactory<BackwardEulerNobleVargheseKohlNoble1998> cell_factory_normal;
        RunAndCheckBidomain(voltage_fast_slow, &cell_factory_normal, "BidomainBackwardToCompareWithFastSlowTest");

        ////In long tests we require a different time step for forward models
        //HeartConfig::Instance()->SetOdeTimeStep(0.01);
        //HeartConfig::Instance()->SetPdeTimeStep(0.01);

        //MyCellFactory<CML_noble_varghese_kohl_noble_1998_basic_pe_lut> cell_factory_lu;
        //RunAndCheckBidomain(voltage_fast_slow, &cell_factory_lu, "BidomainLookupToCompareWithFastSlowTest");

        //MyCellFactory<CML_noble_varghese_kohl_noble_1998_basic> cell_factory_fe;
        //RunAndCheckBidomain(voltage_fast_slow, &cell_factory_fe, "BidomainForwardEulerToCompareWithFastSlowTest");
    }


    // A test where the coarse mesh does not completely cover the fine mesh,
    // to check extrapolation (and resetting if extrapolated values go out
    // of range works ok and doesn't lead to big errors). We use a circle as
    // the fine mesh and a two-element square with corners at north,south,east
    // and west as the coarse mesh
    void TestExtrapolationOnDisk() throw(Exception)
    {
        HeartConfig::Instance()->SetOutputDirectory("MonodomainFastSlowOnDisk");
        HeartConfig::Instance()->SetOutputFilenamePrefix("res");
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
        HeartConfig::Instance()->SetSimulationDuration(5.0); //ms


        MixedTetrahedralMesh<2,2> mixed_mesh;

        // set up the mesh with node (-1/10,0),(0,-1/10),(1/10,0),(0,1/10)
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_2_elements"); // [0,0 - 1,1]
        mixed_mesh.ConstructFromMeshReader(reader);
        mixed_mesh.Scale(sqrt(2)/10,sqrt(2)/10);
        mixed_mesh.Rotate(M_PI/4.0);
        mixed_mesh.Translate(-1.0/10, 0);

        // Disk is radius 1, centred at origin
        TrianglesMeshReader<2,2> fine_mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructFromMeshReader(fine_mesh_reader);
        fine_mesh.Scale(1.0/10, 1.0/10); // now radius 0.1

        mixed_mesh.SetFineMesh(&fine_mesh);

        ReplicatableVector res_fastslow;
        ReplicatableVector res_normal;

        {
            MyCellFactory<FastSlowBackwardEulerNoble98> cell_factory(false);
            MonodomainFastSlowProblem<2> monodomain_fast_slow_prob( &cell_factory, mixed_mesh, 0.1 );
            monodomain_fast_slow_prob.ConvertOutputToMeshalyzerFormat(true);
            monodomain_fast_slow_prob.Initialise();
            monodomain_fast_slow_prob.Solve();
            res_fastslow.ReplicatePetscVector( monodomain_fast_slow_prob.GetSolution() );
        }

        {
            MyCellFactory<BackwardEulerNobleVargheseKohlNoble1998> cell_factory(false);

            MonodomainProblem<2> monodomain_prob( &cell_factory);
            monodomain_prob.ConvertOutputToMeshalyzerFormat(true);
            monodomain_prob.SetMesh(&fine_mesh);
            monodomain_prob.Initialise();
            monodomain_prob.Solve();
            res_normal.ReplicatePetscVector( monodomain_prob.GetSolution() );
        }
        assert(res_fastslow.size()==res_normal.size());

        for(unsigned i=0; i<res_fastslow.size(); i++)
        {
            TS_ASSERT_DELTA( res_fastslow[i], res_normal[i], 2); // was 5e-2 with no extrapolation
        }
    }
};

#endif /*TESTCARDIACDOMAINFASTSLOWPROBLEM_HPP_*/
