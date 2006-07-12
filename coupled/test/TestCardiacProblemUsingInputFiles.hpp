#ifndef TESTCARDIACPROBLEMUSINGINPUTFILES_HPP_
#define TESTCARDIACPROBLEMUSINGINPUTFILES_HPP_


#include "Element.hpp"
#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include <vector>
//#include <iostream>

#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "BidomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "ColumnDataReader.hpp"
#include "CardiacProblemParametersFileReader.hpp"


class CellFactoryUsingFileReader : public AbstractCardiacCellFactory<3>
{
private:
    CardiacProblemParametersFileReader<3>* mpReader;
    InitialStimulus *mpStimulus;
    
public:
    CellFactoryUsingFileReader(CardiacProblemParametersFileReader<3>* pReader) 
            : AbstractCardiacCellFactory<3>(0.01) // mTimetep is set to 0.01 but will be overwritten
    {
        mpReader   = pReader;
        mTimeStep  = mpReader->GetOdeTimeStep();

        mpStimulus = new InitialStimulus( mpReader->GetStimulusMagnitude(),
                                          mpReader->GetStimulusDuration(),
                                          mpReader->GetStimulusStartTime()
                                         );
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus, mpZeroStimulus);
    }
    
    void FinaliseCellCreation(std::vector<AbstractCardiacCell* >* pCellsDistributed, unsigned lo, unsigned hi)
    {
        for(unsigned i=0; i<mpReader->GetStimulatedNodes().size(); i++)
        {
            unsigned global_index = mpReader->GetStimulatedNodes()[i];
            if((global_index>=lo) && (global_index<hi))
            {
                int local_index = global_index - lo;
                (*pCellsDistributed)[ local_index ]->SetStimulusFunction(mpStimulus);
            }
        }
    }
    
    ~CellFactoryUsingFileReader(void)
    {
        delete mpStimulus;
    }
};



class TestCardiacProblemUsingInputFiles : public CxxTest::TestSuite 
{   
 
public:
    void testCardiacProblemUsingInputFiles()  throw(Exception)
    {
        CardiacProblemParametersFileReader<3> reader("coupled/test/data/bidomain_input_file.txt");
        TS_ASSERT_DELTA(reader.GetEndTime(), 1.0, 1e-6); //check the input file hasn't been changed

        CellFactoryUsingFileReader cell_factory(&reader);

        if( reader.IsMonodomainProblem() )
        {
            MonodomainProblem<3> monodomain_problem(&cell_factory);

            monodomain_problem.SetMeshFilename( reader.GetMeshFilename() );
            monodomain_problem.SetOutputDirectory( reader.GetOutputDirectory() );
            monodomain_problem.SetOutputFilenamePrefix( reader.GetOutputFilenamePrefix() );
   
            monodomain_problem.SetEndTime( reader.GetEndTime() );  
            monodomain_problem.SetPdeTimeStep( reader.GetPdeTimeStep() );
            monodomain_problem.SetPrintingTimeStep( reader.GetPrintingTimeStep() );

            monodomain_problem.Initialise();        
            monodomain_problem.Solve();
        }
        else
        {
            BidomainProblem<3> bidomain_problem(&cell_factory);

            bidomain_problem.SetMeshFilename( reader.GetMeshFilename() );
            bidomain_problem.SetOutputDirectory( reader.GetOutputDirectory() );
            bidomain_problem.SetOutputFilenamePrefix( reader.GetOutputFilenamePrefix() );
   
            bidomain_problem.SetEndTime( reader.GetEndTime() );  
            bidomain_problem.SetPdeTimeStep( reader.GetPdeTimeStep() );
            bidomain_problem.SetPrintingTimeStep( reader.GetPrintingTimeStep() );

            if( reader.ThereAreFixedNodes() )
            {
                bidomain_problem.SetFixedExtracellularPotentialNodes( reader.GetFixedNodes() );
            }

            bidomain_problem.Initialise();        
            bidomain_problem.Solve();
        }
    }
};

#endif /*TESTCARDIACPROBLEMUSINGINPUTFILES_HPP_*/
