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
    // We only test the implicit solver as the explicit is not expected to work for very long
    void Test2dImplicit() throw(Exception)
    {
        PlaneStimulusCellFactory<2> cell_factory(0.01, -1000*1000);

        CardiacElectroMechanicsProblem<2> implicit_problem(&cell_factory, 
                                                           10, /* end time */
                                                           5, /*mech mesh size*/ 
                                                           false, /* implicit */
                                                           100, /* 100*0.01ms mech dt */
                                                           0.01,
                                                           "TestCardiacElectroMechImplicit");
        implicit_problem.SetNoElectricsOutput();
        implicit_problem.Solve();

        // test by checking the length of the tissue against hardcoded value
        AbstractElasticityAssembler<2>* assembler = dynamic_cast<AbstractElasticityAssembler<2>*>(implicit_problem.mpCardiacMechAssembler);
        std::vector<Vector<double> >& deformed_position = assembler->rGetDeformedPosition();
        TS_ASSERT_DELTA(deformed_position[0](5), 0.998313, 1e-4);
    }
};
#endif /*TESTCARDIACELECTROMECHANICSPROBLEM_HPP_*/
