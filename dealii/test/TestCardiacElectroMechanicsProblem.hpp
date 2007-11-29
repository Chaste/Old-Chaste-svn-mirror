#ifndef TESTCARDIACELECTROMECHANICSPROBLEM_HPP_
#define TESTCARDIACELECTROMECHANICSPROBLEM_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechanicsProblem1d.hpp"
#include "CardiacElectroMechanicsProblem.hpp"

class TestCardiacElectroMechanicsProblem : public CxxTest::TestSuite
{
public:
    // test the interface works and does what it should do.
    // We only test the implicit solver as the explicit is expected to work for very long
    void Test2dImplicit() throw(Exception)
    {
        PlaneStimulusCellFactory<2> cell_factory(0.01, -1000*1000);

        CardiacElectroMechanicsProblem<2> implicit_problem(&cell_factory, 
                                                           10, /* end time */
                                                           5, /*mech mesh size*/ 
                                                           false, /* implicit */
                                                           100, /* 100*0.01ms mech dt */
                                                           "TestCardiacElectroMechImplicit");
        implicit_problem.SetNoElectricsOutput();
        implicit_problem.Solve();

        // Test by looking at the results manually to see if they look ok and 
        // checking nothing has changed by comparing the log files.
        // 
        // note we have to get rid of the first line in the log (which has the date
        // of the simulation) before doing the comparison.
        OutputFileHandler handler("TestCardiacElectroMechImplicit",false);
        std::string results_dir = handler.GetOutputDirectoryFullPath();         
        std::string command = "sed \"2d\" " + results_dir + "log.txt > " + results_dir + "log2.txt";
        system(command.c_str());
        TS_ASSERT_EQUALS(system(("diff -bB " + results_dir + "log2.txt dealii/test/data/TestCardiacElectroMechImplicit/log.txt").c_str()), 0);     
    }
};
#endif /*TESTCARDIACELECTROMECHANICSPROBLEM_HPP_*/
