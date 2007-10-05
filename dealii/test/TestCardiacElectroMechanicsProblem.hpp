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
//    void TestFindingGuassPoints()
//    {
//    }

    void Test_1D_CompareExplicitVsImplicit() throw(Exception)
    {
        double time_step = 0.01;
        PlaneStimulusCellFactory<1> cell_factory(time_step, -1000*1000);

        // instabilities appear at about 6.8
        CardiacElectroMechanicsProblem1d explicit_problem(&cell_factory, 5, time_step, true,  "ExplicitCardiacElectroMech");
        explicit_problem.Solve();

        CardiacElectroMechanicsProblem1d implicit_problem(&cell_factory, 5, time_step, false, "ImplicitCardiacElectroMech");
        implicit_problem.Solve();
          
        // Get the length of the fibre in both simulations and compare 
        for(unsigned i=0; i<500; i++)
        {
            std::string full_path1 = OutputFileHandler::GetChasteTestOutputDirectory() + 
                                     "ExplicitCardiacElectroMech/deformation/";
            
            // a bit nasty: we want to ready the second column of the last row (turns
            // out that this is where the x value for the 2nd node (ie the X=1) node is,
            // so we copy the last row (using tail) to a temp file, and read that in.
            // There's probably a much nicer way. 
            std::stringstream file1;
            file1 << full_path1 << "/solution_" << i << ".nodes";
            std::string temp_file1 = full_path1 + "/temp.txt";
            system(("tail -1 " + file1.str() + " > " + temp_file1).c_str());

            std::ifstream ifs1(temp_file1.c_str());
            double unused, length_of_fibre1;
            ifs1 >> unused;
            ifs1 >> length_of_fibre1;                 // the second entry is the length
            system(("rm -f " + temp_file1).c_str());

            std::string full_path2 = OutputFileHandler::GetChasteTestOutputDirectory() +
                                     "ImplicitCardiacElectroMech/deformation/";
            
            std::stringstream file2;
            file2 << full_path2 << "/solution_" << i << ".nodes";

            std::string temp_file2 = full_path2 + "/temp.txt";
            system(("tail -1 " + file2.str() + " > " + temp_file2).c_str());

            std::ifstream ifs2(temp_file2.c_str());
            double length_of_fibre2;
            ifs2 >> unused;
            ifs2 >> length_of_fibre2;                // the second entry is the length
            system(("rm -f " + temp_file2).c_str());
            
            if(fabs(length_of_fibre1-1.0)>0.01)
            {
                // must be an error in the read, fibre should be near 1.0 in length
                std::cout << length_of_fibre1;
                TS_FAIL("Probable error in length of fibre");
            }
            
            TS_ASSERT_DELTA(length_of_fibre1, length_of_fibre2, fabs(length_of_fibre1*1e-5));
            
            // hardcoded test
            if(i==450)
            {
                std::cout << "LENGTH = " << length_of_fibre2 << "\n";
                TS_ASSERT_DELTA(length_of_fibre2, 0.999378, 1e-5);
            }
        }
    }
    
    void Test_2D() throw(Exception)
    {
        double time_step = 0.01;
        PlaneStimulusCellFactory<2> cell_factory(time_step, -600*1000);

        CardiacElectroMechanicsProblem<2> problem(&cell_factory, 5, time_step, true,  "CardiacElectroMech2d");
        problem.Solve();
    }
};
#endif /*TESTCARDIACELECTROMECHANICSPROBLEM_HPP_*/
