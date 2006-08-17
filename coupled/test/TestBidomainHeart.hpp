#ifndef _TESTBIDOMAINHEART_HPP_
#define _TESTBIDOMAINHEART_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include <petscvec.h>
#include <vector>

#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"



class PointStimulusHeartCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    InitialStimulus *mpStimulus;
public:
    PointStimulusHeartCellFactory(double timeStep) : AbstractCardiacCellFactory<3>(timeStep)
    {
        mpStimulus = new InitialStimulus(-3000.0*1000, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus, mpZeroStimulus);
    }
    
    void FinaliseCellCreation(std::vector<AbstractCardiacCell* >* pCellsDistributed, unsigned lo, unsigned hi)
    {
        int stimulated_cells[] = {  37484-1,
                                    37499-1,
                                    37777-1,
                                    37779-1,
                                    38008-1,
                                    38332-1,
                                    38587-1,
                                    38588-1,
                                    39312-1,
                                    39314-1,
                                    39643-1,
                                    40588-1,
                                    40590-1,
                                    63885-1
                                 };
                                 
        for (int i=0; i<14; i++)
        {
            int global_index = stimulated_cells[i];
            if ((global_index>=(int)lo) && (global_index<(int)hi))
            {
                int local_index = global_index - lo;
                (*pCellsDistributed)[ local_index ]->SetStimulusFunction(mpStimulus);
            }
        }
    }
    
    ~PointStimulusHeartCellFactory(void)
    {
        delete mpStimulus;
    }
};

class TestBidomainHeart : public CxxTest::TestSuite
{

public:
    void TestBidomainDg0Heart()
    {
        PointStimulusHeartCellFactory cell_factory(0.01);
        BidomainProblem<3> bidomain_problem(&cell_factory);
        
        bidomain_problem.SetMeshFilename("mesh/test/data/heart_fifth");
        bidomain_problem.SetEndTime(100);   // 100 ms
        bidomain_problem.SetOutputDirectory("BiDg0_FifthHeart");
        bidomain_problem.SetOutputFilenamePrefix("BidomainLR91_FifthHeart");
        bidomain_problem.SetPdeTimeStep(0.01);
        bidomain_problem.Initialise();
        
        bidomain_problem.Solve();
    }
};

#endif //_TESTBIDOMAINHEART_HPP_
