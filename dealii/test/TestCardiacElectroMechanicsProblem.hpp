#ifndef TESTCARDIACELECTROMECHANICSPROBLEM_HPP_
#define TESTCARDIACELECTROMECHANICSPROBLEM_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechanicsProblem.hpp"

class TestCardiacElectroMechanicsProblem : public CxxTest::TestSuite
{
public:
    void TestSimple() throw(Exception)
    {
        double time_step = 0.01;
        PlaneStimulusCellFactory<1> cell_factory(time_step, -1000*1000);
        CardiacElectroMechanicsProblem<1> problem(&cell_factory, 1, time_step, "CardiacElectroMech");
        problem.Solve();
    }
};
#endif /*TESTCARDIACELECTROMECHANICSPROBLEM_HPP_*/
