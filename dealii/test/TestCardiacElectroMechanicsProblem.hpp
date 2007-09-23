#ifndef TESTCARDIACELECTROMECHANICSPROBLEM_HPP_
#define TESTCARDIACELECTROMECHANICSPROBLEM_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"
#include "OneDimCardiacElectroMechanicsProblem.hpp"

class TestCardiacElectroMechanicsProblem : public CxxTest::TestSuite
{
public:
    void TestCompareExplicitVsImplicit() throw(Exception)
    {
        double time_step = 0.01;
        PlaneStimulusCellFactory<1> cell_factory(time_step, -1000*1000);

        // instabilities appear at about 6.8
        OneDimCardiacElectroMechanicsProblem<1> explicit_problem(&cell_factory, 5, time_step, true,  "ExplicitCardiacElectroMech");
        explicit_problem.Solve();

        OneDimCardiacElectroMechanicsProblem<1> implicit_problem(&cell_factory, 5, time_step, false, "ImplicitCardiacElectroMech");
        implicit_problem.Solve();
        
        // temporary test - needs rewriting when output format is finalised
        for(unsigned i=0; i<500; i++)
        {
            OutputFileHandler handler("ExplicitCardiacElectroMech",false);
            std::string full_path1 = handler.GetTestOutputDirectory();
            
            std::stringstream file1;
            file1 << full_path1 << "/results_" << i << ".dat";
            
            std::ifstream ifs1(file1.str().c_str());
            double unused, length_of_fibre1;
            ifs1 >> unused;
            ifs1 >> length_of_fibre1;                 // the second entry is the length

            OutputFileHandler handler2("ImplicitCardiacElectroMech",false);
            std::string full_path2 = handler2.GetTestOutputDirectory();
            
            std::stringstream file2;
            file2 << full_path2 << "/results_" << i << ".dat";

            std::ifstream ifs2(file2.str().c_str());
            double length_of_fibre2;
            ifs2 >> unused;
            ifs2 >> length_of_fibre2;                // the second entry is the length
            
            TS_ASSERT_DELTA(length_of_fibre1, length_of_fibre2, fabs(length_of_fibre1*1e-5));
            
            // hardcoded test
            if(i==450)
            {
                std::cout << "LENGTH = " << length_of_fibre2 << "\n";
                TS_ASSERT_DELTA(length_of_fibre2, 0.999378, 1e-5);
            }
        }
    }
};
#endif /*TESTCARDIACELECTROMECHANICSPROBLEM_HPP_*/
