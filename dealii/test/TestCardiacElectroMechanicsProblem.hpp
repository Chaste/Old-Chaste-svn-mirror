#ifndef TESTCARDIACELECTROMECHANICSPROBLEM_HPP_
#define TESTCARDIACELECTROMECHANICSPROBLEM_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechanicsProblem1d.hpp"
#include "CardiacElectroMechanicsProblem.hpp"
//#include "NodewiseData.hpp"

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


//    void TestCinverseDataStructure() throw(Exception)
//    {
//        PlaneStimulusCellFactory<2> cell_factory(0.01, -1000*1000);
//        CardiacElectroMechanicsProblem<2> implicit_problem(&cell_factory, 
//                                                           0.05, /* end time */
//                                                           5, /*mech mesh size*/ 
//                                                           false, /* implicit */
//                                                           1, /* 0.01ms mech dt */
//                                                           0.01,
//                                                           "TestCardiacElectroMechImplicitCinverse");
//        implicit_problem.SetNoElectricsOutput();
//        implicit_problem.Solve();
//        
//        NodewiseData<2>* p_nodewise_data = NodewiseData<2>::Instance();
//        unsigned num_nodes = 97*97; //hardcoded
//        TS_ASSERT_EQUALS(p_nodewise_data->rGetData().size(), num_nodes);
//
//        for(unsigned i=0; i<p_nodewise_data->rGetData().size(); i++)
//        {
//            TS_ASSERT_EQUALS(p_nodewise_data->rGetData()[i].size(), 3u);
//            TS_ASSERT_DELTA(p_nodewise_data->rGetData()[i][0], 1.0, 1e-6);
//            TS_ASSERT_DELTA(p_nodewise_data->rGetData()[i][1], 0.0, 1e-6);
//            TS_ASSERT_DELTA(p_nodewise_data->rGetData()[i][2], 1.0, 1e-6);
//        }
//    }
};
#endif /*TESTCARDIACELECTROMECHANICSPROBLEM_HPP_*/
