#ifndef _TESTMONODOMAINHEART_HPP_
#define _TESTMONODOMAINHEART_HPP_

#include "Element.hpp"
#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include <vector>
//#include <iostream>

#include "ConformingTetrahedralMesh.cpp"
#include "PetscSetupAndFinalize.hpp"
#include "MonodomainPde.hpp"
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"

#include "LuoRudyIModel1991OdeSystem.hpp"

class PointStimulusHeartCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    InitialStimulus *mpStimulus;
public:
    PointStimulusHeartCellFactory(double timeStep) : AbstractCardiacCellFactory<3>(timeStep)
    {
        mpStimulus = new InitialStimulus(-1000.0*1000, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus);
    }
    
    void FinaliseCellCreation(std::vector<AbstractCardiacCell* >* pCellsDistributed, int lo, int hi)
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

        for(int i=0; i<14; i++)
        {
            int global_index = stimulated_cells[i];
            if((global_index>=lo) && (global_index<hi))
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

class TestMonodomainHeart : public CxxTest::TestSuite 
{   
 
public:
    void TestMonodomainDg0Heart()
    {
        PointStimulusHeartCellFactory cell_factory(0.01);
        MonodomainProblem<3> monodomain_problem(&cell_factory);

        monodomain_problem.SetMeshFilename("mesh/test/data/heart_fifth");
        monodomain_problem.SetEndTime(100);   // 100 ms
        monodomain_problem.SetOutputDirectory("/tmp/testoutput/MonoDg0Heart");
        monodomain_problem.SetOutputFilenamePrefix("MonodomainLR91_Heart");
        monodomain_problem.SetPdeTimeStep(0.01);
        monodomain_problem.Initialise();        

        monodomain_problem.Solve();
    }
};

#endif //_TESTMONODOMAINHEART_HPP_
