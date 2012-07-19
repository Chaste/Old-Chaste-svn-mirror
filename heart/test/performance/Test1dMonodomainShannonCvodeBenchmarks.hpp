#ifndef TEST1DMONODOMAINSHANNONCVODEBENCHMARKS_HPP_
#define TEST1DMONODOMAINSHANNONCVODEBENCHMARKS_HPP_

#include <cxxtest/TestSuite.h>
#include <sstream>
#include <iostream>

#include "TetrahedralMesh.hpp"
#include "MonodomainProblem.hpp"
#include "RegularStimulus.hpp"
#include "Shannon2004.hpp"
//#include "Shannon2004BackwardEuler.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "HeartConfig.hpp"
#include "CvodeAdaptor.hpp"
#include "EulerIvpOdeSolver.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * Cell Factory defining cells for 1d chain.
 *
 * This gives a shared default solver to all the cells, as standard.
 */
template <class CELL_MODEL>
class ShannonCardiacCellFactory : public AbstractCardiacCellFactory<1>
{
private:

	boost::shared_ptr<RegularStimulus> mpStimulus;

public:
	ShannonCardiacCellFactory()
	: AbstractCardiacCellFactory<1>(),
	  mpStimulus(new RegularStimulus(-250000,5,2000,1))
	  {
	  }

	AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
	{
		CELL_MODEL *p_cell;

		if (this->GetMesh()->GetNode(node)->GetPoint()[0] == 0.0)
		{
			p_cell = new CELL_MODEL(this->mpSolver,
					                mpStimulus);
		}
		else
		{
			p_cell = new CELL_MODEL(this->mpSolver,
					                this->mpZeroStimulus);
		}

		return p_cell;
	}

};

/**
 *
 * Cell Factory defining Shannon cells for 1d chain.
 *
 * In this version each node has its own Cvode solver
 */
class ShannonCvodeCellFactory : public AbstractCardiacCellFactory<1>
{
private:

    boost::shared_ptr<RegularStimulus> mpStimulus;

public:
    ShannonCvodeCellFactory()
    : AbstractCardiacCellFactory<1>(),
      mpStimulus(new RegularStimulus(-250000,5,2000,1))
      {
      }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        AbstractCardiacCell* p_cell;

        /*
         * HOW_TO_TAG Cardiac/Problem definition
         * Use a CVODE solver in a tissue simulation
         */

        // Here we add a new Cvode Adaptor solver for every node.
        boost::shared_ptr<CvodeAdaptor> cvode_solver(new CvodeAdaptor());
        cvode_solver->SetMinimalReset(true);

        if (this->GetMesh()->GetNode(node)->GetPoint()[0] == 0.0)
        {
            p_cell = new CellShannon2004FromCellML(cvode_solver,
                                                   mpStimulus);
        }
        else
        {
            p_cell = new CellShannon2004FromCellML(cvode_solver,
                                                   this->mpZeroStimulus);
        }

        return p_cell;
    }

};

/**
 * This class tests that a CVODE tissue simulation gets comparable results to those
 * of Forward and Backward Euler.
 *
 * At present FE and BE are more different than FE and CVODE.
 */
class Test1dMonodomainShannonCvodeBenchmarks : public CxxTest::TestSuite
{
private:
    bool CompareBenchmarkResults(const std::vector<double>& rReferenceTrace, const std::vector<double>& rTestTrace, double tol)
    {
        assert(rReferenceTrace.size()==rTestTrace.size());

        // See what the maximum difference in voltages recorded at final time at this node is
        unsigned last_time = rReferenceTrace.size() - 1u;
        double difference = fabs(rReferenceTrace[last_time] - rTestTrace[last_time]);

        // I was looping over all of the time trace here, but conduction velocity must
        // be slightly affected by the different solvers, giving rise to massive
        // (height of upstroke ~135mV) differences, so test at 'plateau' instead.

        // Return a false if it is too big.
        bool result = true;
        if (difference > tol)
        {
            std::cout << "\nDifference in recorded voltages = " << difference << std::endl << std::flush;
            result = false;
        }

        // Also return a false if the cell didn't get stimulated (assuming test only runs to plateau).
        if (rTestTrace[last_time] < 0)
        {
            std::cout << "Cell is not depolarized\n" << std::endl << std::flush;
            result = false;
        }

        return result;
    }

public:

    void TestWithDifferentCellsAndSolvers() throw(Exception)
    {
        double duration = 25;
        HeartConfig::Instance()->SetSimulationDuration(duration); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        double pde_time_step = 0.1; //ms
        double printing_time_step = 0.1; //ms

        std::vector<double> fe_node_0;
        std::vector<double> fe_node_100;
        std::vector<double> times;

        {
            HeartConfig::Instance()->SetOutputDirectory("ShannonBenchmark/forward_euler");
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0025,pde_time_step,printing_time_step);
            ShannonCardiacCellFactory<CellShannon2004FromCellML> cell_factory;
            MonodomainProblem<1> monodomain_problem( &cell_factory );

            monodomain_problem.Initialise();
            double start_time = std::clock();
            monodomain_problem.Solve();
            double end_time = std::clock();
            double elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
            std::cout << "1. Forward Euler elapsed time = " << elapsed_time << " secs for " << duration << " ms\n";

            Hdf5DataReader data_reader = monodomain_problem.GetDataReader();
            times = data_reader.GetUnlimitedDimensionValues();
            fe_node_0 = data_reader.GetVariableOverTime("V", 0);
            fe_node_100 = data_reader.GetVariableOverTime("V", 100);
        }

        /*
         * To test backward Euler you need to generate an extra Chaste variant of the model from the CellML.
         *
         * i.e. go to heart/src/odes/cellml/Shannon2004-conf.xml and add
         * <arg>--backward-euler</arg>
         *
         * then uncomment the //#include "Shannon2004BackwardEuler.hpp" and the below block of code.
         *
         * Same idea for checking optimised models.
         */

//        {
//            HeartConfig::Instance()->SetOutputDirectory("ShannonBenchmark/backward_euler");
//            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(pde_time_step,pde_time_step,printing_time_step);
//            ShannonCardiacCellFactory<CellShannon2004FromCellMLBackwardEuler> cell_factory;
//            MonodomainProblem<1> monodomain_problem( &cell_factory );
//
//            monodomain_problem.Initialise();
//            double start_time = std::clock();
//            monodomain_problem.Solve();
//            double end_time = std::clock();
//            double elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
//            std::cout << "1b. Backward Euler elapsed time = " << elapsed_time << " secs for " << duration << " ms\n";
//
//            Hdf5DataReader data_reader = monodomain_problem.GetDataReader();
//            std::vector<double> be_node_0 = data_reader.GetVariableOverTime("V", 0);
//            std::vector<double> be_node_100 = data_reader.GetVariableOverTime("V", 100);
//
//            TS_ASSERT(CompareBenchmarkResults(fe_node_0, be_node_0, 1.1));     // More than a milliVolt of difference with B.E.
//            TS_ASSERT(CompareBenchmarkResults(fe_node_100, be_node_100, 1.1)); // More than a milliVolt of difference with B.E.
//        }

        {
            HeartConfig::Instance()->SetOutputDirectory("ShannonBenchmark/cvode");
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(pde_time_step,pde_time_step,printing_time_step);
            ShannonCvodeCellFactory cell_factory;
            MonodomainProblem<1> monodomain_problem( &cell_factory );

            monodomain_problem.Initialise();
            double start_time = std::clock();
            monodomain_problem.Solve();
            double end_time = std::clock();
            double elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
            std::cout << "2. CVODE adaptor elapsed time = " << elapsed_time << " secs for " << duration << " ms\n";

            Hdf5DataReader data_reader = monodomain_problem.GetDataReader();
            std::vector<double> cvode_node_0 = data_reader.GetVariableOverTime("V", 0);
            std::vector<double> cvode_node_100 = data_reader.GetVariableOverTime("V", 100);

            TS_ASSERT(CompareBenchmarkResults(fe_node_0, cvode_node_0, 0.4));     // Only 0.2 -- 0.4mV difference with CVODE
            TS_ASSERT(CompareBenchmarkResults(fe_node_100, cvode_node_100, 0.4)); // Only 0.2 -- 0.4mV difference with CVODE
        }

        {
            /// \todo #2116 Native CVODE timings.
        }

    }

};

#endif /*TEST1DMONODOMAINSHANNONCVODEBENCHMARKS_HPP_*/
