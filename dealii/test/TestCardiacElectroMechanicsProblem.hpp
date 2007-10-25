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
    void Test1dCompareExplicitVsImplicit() throw(Exception)
    {
        double time_step = 0.01;
        PlaneStimulusCellFactory<1> cell_factory(time_step, -1000*1000);

        // put these inside braces so that there destructor gets called and all output
        // files get closed
        {
            // instabilities appear at about 6.8
            CardiacElectroMechanicsProblem1d explicit_problem(&cell_factory, 0.5, time_step, true,  "ExplicitCardiacElectroMech");
            explicit_problem.Solve();

            CardiacElectroMechanicsProblem1d implicit_problem(&cell_factory, 0.5, time_step, false, "ImplicitCardiacElectroMech");
            implicit_problem.Solve();
        }

        std::string file1 = OutputFileHandler::GetChasteTestOutputDirectory() + 
                            "ExplicitCardiacElectroMech/length.txt";

        std::string file2 = OutputFileHandler::GetChasteTestOutputDirectory() + 
                            "ImplicitCardiacElectroMech/length.txt";
        
        std::ifstream ifs1(file1.c_str());
        std::ifstream ifs2(file2.c_str());
        
        double length1;
        double length2;
        
        ifs1 >> length1;
        ifs2 >> length2;

        double counter = 0;
        
        TS_ASSERT(!ifs1.eof());
        TS_ASSERT(!ifs2.eof());
        
        while(!ifs1.eof())
        {
            std::cout << length1 << " " << length2 << "\n";
            TS_ASSERT_DELTA(length1, length2,  fabs(length1*1e-5));
            ifs1 >> length1;
            ifs2 >> length2;
            
            // hardcoded test
            if(counter==450)
            {
                TS_ASSERT_DELTA(length2, 0.999378, 1e-5);
            }
            counter++;
        }
    }
    
    void Test2dCompareExplicitVsImplicit() throw(Exception)
    {
        double time_step = 0.01;
        PlaneStimulusCellFactory<2> cell_factory(time_step, -1000*1000);

        CardiacElectroMechanicsProblem<2> explicit_problem(&cell_factory, 100, time_step, true, 40, 16, "CardiacElectroMech2dExplicit");
        explicit_problem.Solve();

//        CardiacElectroMechanicsProblem<2> implicit_problem(&cell_factory, 100, time_step, false, 40, 16, "CardiacElectroMech2dImplicit");
  //      implicit_problem.Solve();
    }
};
#endif /*TESTCARDIACELECTROMECHANICSPROBLEM_HPP_*/
