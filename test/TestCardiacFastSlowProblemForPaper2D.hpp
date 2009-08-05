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

#ifndef TESTCARDIACDOMAINFASTSLOWPROBLEMFORPAPER2D_HPP_
#define TESTCARDIACDOMAINFASTSLOWPROBLEMFORPAPER2D_HPP_

#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "SimpleStimulus.hpp"
#include "FastSlowLuoRudyIModel1991.hpp"
#include "FastSlowBackwardEulerNoble98.hpp"
#include "BackwardEulerNobleVargheseKohlNoble1998.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "MonodomainFastSlowProblem.hpp"
#include "BidomainFastSlowProblem.hpp"
#include "NobleVargheseKohlNoble1998.hpp"
#include "NobleVargheseKohlNoble1998Optimised.hpp"
#include "CellProperties.hpp"
#include "PropagationPropertiesCalculator.hpp"

// simple cell factory that creates fast-slow cells.
template <class CELL>
class MyCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:

    MyCellFactory() :
        AbstractCardiacCellFactory<2>(),
        mpStimulus(new SimpleStimulus(-10000*250.0, 0.2))        //Where did this come from?
        //Will it work with LR91 and N98?
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
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
};


class TestCardiacFastSlowProblemForPaper2D : public CxxTest::TestSuite
{
private:
    static const bool mRunForwardTests=false;
    //Make sure that Mono/Bi are solving with the same parameters
    double mSlowOdeTimestep;
    MixedTetrahedralMesh<2,2> mMixedMesh;
    TetrahedralMesh<2,2> *mpFineMesh;
    unsigned mMiddleIndex, mRhsIndex;
    std::vector<unsigned> mOutputNodes;

    /** Physiological comparison variables */
    double mConductionVelocityFastSlow;
    double mApd90FastSlow;
    double mMaxUpstrokeVelocityFastSlow;
    double mProbeDistance;

    void SetHeartConfigForTest()
    {
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");

        //Found this figure of absolute tolerance by running a standard Monodomain
        //simulation with default KSP parameters (reltol etc) and then finding an
        //absolute tolerance with the same CPU time (same number of KSP iterations)
        HeartConfig::Instance()->SetUseAbsoluteTolerance(2e-4);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1, 0.1, 0.1);
        HeartConfig::Instance()->SetSimulationDuration(350); //ms
        mSlowOdeTimestep = 1.0;
    }

    void DumpVector(std::vector<double> &voltage, std::string name)
    {
        if (PetscTools::AmMaster())
        {
            OutputFileHandler handler(HeartConfig::Instance()->GetOutputDirectory()+"/output", false);
            out_stream file= handler.OpenOutputFile(name);
            for (unsigned i=0;i<voltage.size();i++)
            {
                (*file)<<i*HeartConfig::Instance()->GetPrintingTimeStep()<<"\t";
                (*file)<<voltage[i]<<"\n";
            }
            file->close();
        }
        MPI_Barrier(PETSC_COMM_WORLD);
    }



    void CheckPropertiesAndDump(Hdf5DataReader &rResultsReader)
    {
        std::vector<double> voltage_normal=rResultsReader.GetVariableOverTime("V", mMiddleIndex);
        DumpVector(voltage_normal, "middle_node.dat");

        std::vector<double> voltage_rhs=rResultsReader.GetVariableOverTime("V", mRhsIndex);
        DumpVector(voltage_rhs, "rhs_node.dat");

        double max_upstroke_velocity = 0.0;
        double conduction_velocity = 0.0;
        double apd_90 = 0.0;

        PropagationPropertiesCalculator ppc(&rResultsReader);
        conduction_velocity=ppc.CalculateConductionVelocity(mMiddleIndex, mRhsIndex, mProbeDistance);
        max_upstroke_velocity = ppc.CalculateMaximumUpstrokeVelocity(mMiddleIndex);
        apd_90 = ppc.CalculateActionPotentialDuration(90, mMiddleIndex);

        if (PetscTools::AmMaster())
        {
            std::cout<<"Normal test\n";
            std::cout<<"Conduction Velocity = "<<conduction_velocity<<"\n";
            std::cout<<"Max Upstroke Velocity = "<<max_upstroke_velocity<<"\n";
            std::cout<<"APD 90 = "<<apd_90<<"\n";
        }
        TS_ASSERT_DELTA(mConductionVelocityFastSlow, conduction_velocity, 0.05);
        TS_ASSERT_DELTA(mMaxUpstrokeVelocityFastSlow, max_upstroke_velocity, 35.0);
        TS_ASSERT_DELTA(mApd90FastSlow, apd_90, 2.0);
    }

    void RunAndCheckMonodomain(AbstractCardiacCellFactory<2>* pCellFactory,
                               std::string outputDirectory)
    {

        MonodomainProblem<2> monodomain_prob(pCellFactory);
        HeartConfig::Instance()->SetOutputDirectory(outputDirectory);
        monodomain_prob.SetMesh(mpFineMesh);
        monodomain_prob.SetOutputNodes(mOutputNodes);
        monodomain_prob.ConvertOutputToMeshalyzerFormat(true);
        monodomain_prob.Initialise();
        monodomain_prob.Solve();
        HeartEventHandler::Headings();
        HeartEventHandler::Report();
        HeartEventHandler::Reset();
        /* Read the vector at a typical downstream node */
        Hdf5DataReader results_reader=monodomain_prob.GetDataReader();
        CheckPropertiesAndDump(results_reader);
    }

    void RunAndCheckBidomain(AbstractCardiacCellFactory<2>* pCellFactory,
                             std::string outputDirectory)
    {
        std::vector<double> voltage_normal;
        HeartConfig::Instance()->SetOutputDirectory(outputDirectory);

        BidomainProblem<2> bidomain_prob( pCellFactory);
        bidomain_prob.SetMesh(mpFineMesh);
        bidomain_prob.SetOutputNodes(mOutputNodes);
        bidomain_prob.ConvertOutputToMeshalyzerFormat(true);
        bidomain_prob.Initialise();
        bidomain_prob.Solve();
        HeartEventHandler::Headings();
        HeartEventHandler::Report();
        HeartEventHandler::Reset();
        /* Read the vector at a typical downstream node */
        Hdf5DataReader results_reader=bidomain_prob.GetDataReader();
        CheckPropertiesAndDump(results_reader);
    }

public:
    void TestConstructMixedMesh()
    {
        double width = 2; //cm
        double height = 2; //cm
        double coarse_ds = 0.1; //cm
        double fine_ds = 0.01;

        unsigned num_coarse_elem_each_dir = unsigned(width/coarse_ds + 0.5);
        unsigned num_fine_elem_each_dir = unsigned(width/fine_ds + 0.5);

        TS_ASSERT_EQUALS(num_coarse_elem_each_dir, 20U);
        TS_ASSERT_EQUALS(num_fine_elem_each_dir, 200U);
        mMixedMesh.ConstructRectangularMeshes(width, height, num_coarse_elem_each_dir, num_fine_elem_each_dir);
        mpFineMesh=mMixedMesh.GetFineMesh();

        TS_ASSERT_EQUALS(mpFineMesh->GetNumNodes(), 40401U);
        TS_ASSERT_EQUALS(mMixedMesh.GetNumNodes(), 441U);

        unsigned node_ratio = mpFineMesh->GetNumNodes() / mMixedMesh.GetNumNodes(); // Note -- rounded
        TS_ASSERT_EQUALS(node_ratio, 91U);//Rounded down


        //Find a place to probe in the middle of the mesh and on the middle of the right-hand side
        mMiddleIndex = UINT_MAX;
        mRhsIndex  = UINT_MAX;
        for (unsigned i=0; i<mpFineMesh->GetNumNodes(); i++){
            double y=mpFineMesh->GetNode(i)->rGetLocation()[1];
            if ( fabs(y - (height/2.0)) < 1e-6)
            {
                //Half way up
                double x=mpFineMesh->GetNode(i)->rGetLocation()[0];
                if ( fabs (x - width) < 1e-6)
                {
                    //On righthand side
                    mRhsIndex = i;
                }
                if ( fabs (x - width/2.0) < 1e-6)
                {
                    //Bang in the middle
                    mMiddleIndex = i;
                }
            }
        }
        TS_ASSERT(mMiddleIndex != UINT_MAX);
        TS_ASSERT(mRhsIndex != UINT_MAX);
        TS_ASSERT_EQUALS(mMiddleIndex, 20200U);
        TS_ASSERT_EQUALS(mRhsIndex, 20300U);

        mProbeDistance = mpFineMesh->GetDistanceBetweenNodes(mMiddleIndex, mRhsIndex);
        TS_ASSERT_DELTA(mProbeDistance, 1.0, 1e-10);

        mOutputNodes.push_back(mMiddleIndex);
        mOutputNodes.push_back(mRhsIndex);

        SetHeartConfigForTest();

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
        HeartConfig::Instance()->SetOutputDirectory("BidomainFastSlow");
        HeartConfig::Instance()->SetOutputFilenamePrefix("res");

        {
            MyCellFactory<FastSlowBackwardEulerNoble98> cell_factory;

            BidomainFastSlowProblem<2> bidomain_fast_slow_prob( &cell_factory, mMixedMesh, mSlowOdeTimestep );
            bidomain_fast_slow_prob.SetOutputNodes(mOutputNodes);
            bidomain_fast_slow_prob.ConvertOutputToMeshalyzerFormat(true);
            bidomain_fast_slow_prob.Initialise();
            bidomain_fast_slow_prob.Solve();

            HeartEventHandler::Headings();
            HeartEventHandler::Report();
            HeartEventHandler::Reset();

           /* Read the vector at a typical downstream node */
            Hdf5DataReader results_reader=bidomain_fast_slow_prob.GetDataReader();
            std::vector<double> voltage_rhs=results_reader.GetVariableOverTime("V", mRhsIndex);
            DumpVector(voltage_rhs, "rhs_node.dat");

            std::vector<double> voltage_fast_slow=results_reader.GetVariableOverTime("V", mMiddleIndex);
            DumpVector(voltage_fast_slow, "middle_node.dat");

            unsigned expected_size = (unsigned)(HeartConfig::Instance()->GetSimulationDuration()/HeartConfig::Instance()->GetPrintingTimeStep() ) + 1 ;
            TS_ASSERT_EQUALS(voltage_fast_slow.size(),  expected_size);

            PropagationPropertiesCalculator ppc(&results_reader);
            mConductionVelocityFastSlow=ppc.CalculateConductionVelocity(mMiddleIndex, mRhsIndex, mProbeDistance);
            mMaxUpstrokeVelocityFastSlow = ppc.CalculateMaximumUpstrokeVelocity(mMiddleIndex);
            mApd90FastSlow = ppc.CalculateActionPotentialDuration(90, mMiddleIndex);

            if (PetscTools::AmMaster())
            {
                std::cout<<"Fast slow test\n";
                std::cout<<"Conduction Velocity = "<<mConductionVelocityFastSlow<<"\n";
                std::cout<<"Max Upstroke Velocity = "<<mMaxUpstrokeVelocityFastSlow<<"\n";
                std::cout<<"APD 90 = "<<mApd90FastSlow<<"\n";
            }
        }

        //////////////////////////////////////////////////////////
        // solve using normal bidomain problem - backward model
        /////////////////////////////////////////////////////////

        MyCellFactory<BackwardEulerNobleVargheseKohlNoble1998> cell_factory_normal;
        RunAndCheckBidomain(&cell_factory_normal, "BidomainBackwardToCompareWithFastSlow");


        if (mRunForwardTests)
        {
            ////In long tests we require a different time step for forward models
            HeartConfig::Instance()->SetOdeTimeStep(0.01);
            HeartConfig::Instance()->SetPdeTimeStep(0.01);

            MyCellFactory<CML_noble_varghese_kohl_noble_1998_basic_pe_lut> cell_factory_lu;
            RunAndCheckBidomain(&cell_factory_lu, "BidomainLookupToCompareWithFastSlow");

            MyCellFactory<CML_noble_varghese_kohl_noble_1998_basic> cell_factory_fe;
            RunAndCheckBidomain(&cell_factory_fe, "BidomainForwardEulerToCompareWithFastSlow");
        }
        else
        {
            if (PetscTools::AmMaster())
            {
                std::cout<<"Warning: Foward integration tests are disabled\n";
            }
        }
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
        HeartConfig::Instance()->SetOutputDirectory("MonodomainFastSlow");
        HeartConfig::Instance()->SetOutputFilenamePrefix("res");

         {
            MyCellFactory<FastSlowBackwardEulerNoble98> cell_factory;

            MonodomainFastSlowProblem<2> monodomain_fast_slow_prob( &cell_factory, mMixedMesh, mSlowOdeTimestep );
            monodomain_fast_slow_prob.SetOutputNodes(mOutputNodes);
            monodomain_fast_slow_prob.ConvertOutputToMeshalyzerFormat(true);
            monodomain_fast_slow_prob.Initialise();
            monodomain_fast_slow_prob.Solve();

            HeartEventHandler::Headings();
            HeartEventHandler::Report();
            HeartEventHandler::Reset();

            /* Read the vector at a typical downstream node */
            Hdf5DataReader results_reader=monodomain_fast_slow_prob.GetDataReader();
            std::vector <double> voltage_rhs=results_reader.GetVariableOverTime("V", mRhsIndex);
            DumpVector(voltage_rhs, "rhs_node.dat");

            std::vector<double> voltage_fast_slow=results_reader.GetVariableOverTime("V", mMiddleIndex);
            DumpVector(voltage_fast_slow, "middle_node.dat");

            unsigned expected_size = (unsigned)(HeartConfig::Instance()->GetSimulationDuration()/HeartConfig::Instance()->GetPrintingTimeStep() ) + 1 ;
            TS_ASSERT_EQUALS(voltage_fast_slow.size(),  expected_size);

            PropagationPropertiesCalculator ppc(&results_reader);
            mConductionVelocityFastSlow = ppc.CalculateConductionVelocity(mMiddleIndex, mRhsIndex, mProbeDistance);
            mMaxUpstrokeVelocityFastSlow = ppc.CalculateMaximumUpstrokeVelocity(mMiddleIndex);
            mApd90FastSlow = ppc.CalculateActionPotentialDuration(90, mMiddleIndex);

            if (PetscTools::AmMaster())
            {
                std::cout<<"Fast slow test\n";
                std::cout<<"Conduction Velocity = "<<mConductionVelocityFastSlow<<"\n";
                std::cout<<"Max Upstroke Velocity = "<<mMaxUpstrokeVelocityFastSlow<<"\n";
                std::cout<<"APD 90 = "<<mApd90FastSlow<<"\n";
            }

        }

        /////////////////////////////////////////////////////////
        // solve using normal monodomain problem - backward model
        /////////////////////////////////////////////////////////
        MyCellFactory<BackwardEulerNobleVargheseKohlNoble1998> cell_factory_backward;
        RunAndCheckMonodomain(&cell_factory_backward, "MonodomainBackwardToCompareWithFastSlow");

        if (mRunForwardTests)
        {
            ////In long tests we require a different time step for forward models
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);

            MyCellFactory<CML_noble_varghese_kohl_noble_1998_basic_pe_lut> cell_factory_lu;
            RunAndCheckMonodomain(&cell_factory_lu, "MonodomainLookupToCompareWithFastSlow");

            MyCellFactory<CML_noble_varghese_kohl_noble_1998_basic> cell_factory_fe;
            RunAndCheckMonodomain(&cell_factory_fe, "MonodomainForwardEulerToCompareWithFastSlow");
        }
        else
        {
            if (PetscTools::AmMaster())
            {
                std::cout<<"Warning: Foward integration tests are disabled\n";
            }
        }
    }
};

#endif /*TESTCARDIACDOMAINFASTSLOWPROBLEMFORPAPER2D_HPP_*/
