#ifndef _TESTMONODOMAINHEART_HPP_
#define _TESTMONODOMAINHEART_HPP_

#include "Element.hpp"

#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include <vector>
//#include <iostream>

#include "ConformingTetrahedralMesh.cpp"
#include "PetscSetupAndFinalize.hpp"
#include "MonodomainPdeIteration7.hpp"
#include "MonodomainProblemIteration7.hpp"
#include "AbstractCardiacCellFactory.hpp"


class PointStimulusHeartCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    InitialStimulus *mpStimulus;
public:
    PointStimulusHeartCellFactory(double timeStep) : AbstractCardiacCellFactory<3>(timeStep)
    {
        mpStimulus = new InitialStimulus(-300.0, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(int node)
    {
        return new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus, mTimeStep);
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
            if((stimulated_cells[i]>=lo) && (stimulated_cells[i]<hi))
            {
                (*pCellsDistributed)[ stimulated_cells[i] - lo ]->SetStimulusFunction(mpStimulus);
            }
        }
    }
    
    ~PointStimulusHeartCellFactory(void)
    {
        delete mpStimulus;
    }
};

/*
class PointStimulusHeart: public AbstractMonodomainProblemStimulus<3>
{
public:
    virtual void Apply(MonodomainPde<3> *pPde, 
                       ConformingTetrahedralMesh<3,3> *pMesh)
    {
        // Nodes to apply stimuus to on apex of heart found from 
        // Tulane data
        static InitialStimulus stimulus(-300.0, 0.5);
        pPde->SetStimulusFunctionAtNode(37484-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(37499-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(37777-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(37779-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(38008-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(38332-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(38587-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(38588-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(39312-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(39314-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(39643-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(40588-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(40590-1, &stimulus); 
        pPde->SetStimulusFunctionAtNode(63885-1, &stimulus); 

    }
};
*/
class TestMonodomainHeart : public CxxTest::TestSuite 
{   

    
public:
    void TestMonodomainDg0Heart()
    {
        PointStimulusHeartCellFactory cell_factory(0.005);
        MonodomainProblemIteration7<3> monodomain_problem(&cell_factory);

        monodomain_problem.SetMeshFilename("mesh/test/data/heart");
        monodomain_problem.SetEndTime(100);   // 100 ms
        monodomain_problem.SetOutputDirectory("testoutput/MonoDg0Heart");
        monodomain_problem.SetOutputFilenamePrefix("MonodomainLR91_Heart");
        monodomain_problem.SetPdeTimeStep(0.01);
        monodomain_problem.Initialise();        

        monodomain_problem.Solve();
        
        double* voltage_array;
    
        VecRestoreArray(monodomain_problem.mCurrentVoltage, &voltage_array);      
        VecAssemblyBegin(monodomain_problem.mCurrentVoltage);
        VecAssemblyEnd(monodomain_problem.mCurrentVoltage);
        VecDestroy(monodomain_problem.mCurrentVoltage);
    }
};

#endif //_TESTMONODOMAINHEART_HPP_
