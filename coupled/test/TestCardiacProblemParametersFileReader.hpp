#ifndef TESTCARDIACPROBLEMPARAMETERSFILEREADER_HPP_
#define TESTCARDIACPROBLEMPARAMETERSFILEREADER_HPP_


// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"
#include <cxxtest/TestSuite.h>
#include <petscvec.h>
//#include <iostream>
#include "CardiacProblemParametersFileReader.hpp"


class TestCardiacProblemParametersFileReader : public CxxTest::TestSuite 
{
public :
    void testCardiacProblemParametersFileReader() throw(Exception)
    {
        CardiacProblemParametersFileReader<3> reader("coupled/test/data/sample_input_file.txt");
 
        TS_ASSERT_EQUALS(reader.IsMonodomainProblem(), false);
 
        TS_ASSERT_DELTA(reader.GetSurfaceAreaToVolumeRatio(), 1500, 1e-10);
        TS_ASSERT_DELTA(reader.GetCapacitance(), 1.1, 1e-10);
    
        TS_ASSERT_DELTA(reader.GetPdeTimeStep(), 0.01, 1e-10);
        TS_ASSERT_DELTA(reader.GetOdeTimeStep(), 0.001, 1e-10);
        TS_ASSERT_DELTA(reader.GetPrintingTimeStep(), 0.1, 1e-10 );
        TS_ASSERT_DELTA(reader.GetEndTime(), 1, 1e-10);

        TS_ASSERT_EQUALS(reader.GetMeshFilename(),         "some_mesh");
        TS_ASSERT_EQUALS(reader.GetOutputDirectory(),      "some_dir");
        TS_ASSERT_EQUALS(reader.GetOutputFilenamePrefix(), "some_prefix");
        
        TS_ASSERT_EQUALS(reader.ThereAreFixedNodes(), false);
        
        TS_ASSERT_EQUALS(reader.GetStimulusStartTime(), 0.1);
        TS_ASSERT_EQUALS(reader.GetStimulusDuration(),  0.4);
        TS_ASSERT_EQUALS(reader.GetStimulusMagnitude(), -1e6);
        
        std::vector<unsigned> stimulated_nodes = reader.GetStimulatedNodes();
        TS_ASSERT_EQUALS(stimulated_nodes.size(), 1);
        TS_ASSERT_EQUALS(stimulated_nodes[0],  0);

            
        c_matrix<double,3,3> sigma_i = reader.GetIntracellularConductivity();
        c_matrix<double,3,3> sigma_e = reader.GetExtracellularConductivity();

        TS_ASSERT_DELTA( sigma_i(0,0), 1.1, 1e-10 );
        TS_ASSERT_DELTA( sigma_i(0,1), 0.1, 1e-10 );
        TS_ASSERT_DELTA( sigma_i(0,2), 0.2, 1e-10 );
        TS_ASSERT_DELTA( sigma_i(1,0), 0.1, 1e-10 );
        TS_ASSERT_DELTA( sigma_i(1,1), 1.2, 1e-10 );
        TS_ASSERT_DELTA( sigma_i(1,2), 0.3, 1e-10 );
        TS_ASSERT_DELTA( sigma_i(2,0), 0.2, 1e-10 );
        TS_ASSERT_DELTA( sigma_i(2,1), 0.3, 1e-10 );
        TS_ASSERT_DELTA( sigma_i(2,2), 1.3, 1e-10 );

        TS_ASSERT_DELTA( sigma_e(0,0), 2.1, 1e-10 );
        TS_ASSERT_DELTA( sigma_e(0,1), 0.4, 1e-10 );
        TS_ASSERT_DELTA( sigma_e(0,2), 0.5, 1e-10 );
        TS_ASSERT_DELTA( sigma_e(1,0), 0.4, 1e-10 );
        TS_ASSERT_DELTA( sigma_e(1,1), 2.2, 1e-10 );
        TS_ASSERT_DELTA( sigma_e(1,2), 0.6, 1e-10 );
        TS_ASSERT_DELTA( sigma_e(2,0), 0.5, 1e-10 );
        TS_ASSERT_DELTA( sigma_e(2,1), 0.6, 1e-10 );
        TS_ASSERT_DELTA( sigma_e(2,2), 2.3, 1e-10 );
        
        
        // test the conductivities for 2d and 1d (using same datafile - in 2d the reader will
        // read the first three conductivity values and ignore the last two, in 1d only the first
        // will be read
        CardiacProblemParametersFileReader<2> reader2d("coupled/test/data/sample_input_file.txt");
        CardiacProblemParametersFileReader<1> reader1d("coupled/test/data/sample_input_file.txt");

        c_matrix<double,2,2> sigma_i_2d = reader2d.GetIntracellularConductivity();
        c_matrix<double,2,2> sigma_e_2d = reader2d.GetExtracellularConductivity();

        c_matrix<double,1,1> sigma_i_1d = reader1d.GetIntracellularConductivity();
        c_matrix<double,1,1> sigma_e_1d = reader1d.GetExtracellularConductivity();

        TS_ASSERT_DELTA( sigma_i_2d(0,0), 1.1, 1e-10 );
        TS_ASSERT_DELTA( sigma_i_2d(0,1), 0.1, 1e-10 );
        TS_ASSERT_DELTA( sigma_i_2d(1,0), 0.1, 1e-10 );
        TS_ASSERT_DELTA( sigma_i_2d(1,1), 0.2, 1e-10 );  

        TS_ASSERT_DELTA( sigma_e_2d(0,0), 2.1, 1e-10 );
        TS_ASSERT_DELTA( sigma_e_2d(0,1), 0.4, 1e-10 );
        TS_ASSERT_DELTA( sigma_e_2d(1,0), 0.4, 1e-10 );
        TS_ASSERT_DELTA( sigma_e_2d(1,1), 0.5, 1e-10 );
        
        TS_ASSERT_DELTA( sigma_i_1d(0,0), 1.1, 1e-10 );

        TS_ASSERT_DELTA( sigma_e_1d(0,0), 2.1, 1e-10 );
        
        
        // previous input file had "NONE" for fixed nodes, also test input file which
        // does specify a fixed node file. This one says the fixed nodes are in fixed_nodes.txt
        CardiacProblemParametersFileReader<3> reader2("coupled/test/data/another_input_file.txt");
        TS_ASSERT_EQUALS(reader2.ThereAreFixedNodes(), true);
        std::vector<unsigned> fixed_nodes = reader2.GetFixedNodes();
        TS_ASSERT_EQUALS(fixed_nodes.size(), 4);
        TS_ASSERT_EQUALS(fixed_nodes[0], 10);
        TS_ASSERT_EQUALS(fixed_nodes[1], 11);
        TS_ASSERT_EQUALS(fixed_nodes[2], 12);
        TS_ASSERT_EQUALS(fixed_nodes[3], 13);
        
        // this reader should not have read the optional physiological constants 
        TS_ASSERT_EQUALS(reader2.IntracellularConductivityWasSet(), false);
        TS_ASSERT_EQUALS(reader2.ExtracellularConductivityWasSet(), false);
        TS_ASSERT_EQUALS(reader2.SurfaceAreaToVolumeRatioWasSet(), false);
        TS_ASSERT_EQUALS(reader2.CapacitanceWasSet(), false);
        
    }
};
#endif /*TESTCARDIACPROBLEMPARAMETERSFILEREADER_HPP_*/
