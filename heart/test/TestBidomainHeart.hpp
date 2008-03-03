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
        mpStimulus = new InitialStimulus(-1000.0*1000, 0.5);
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


class PointStimulusHeartCellFactoryMetis : public AbstractCardiacCellFactory<3>
{
private:
    InitialStimulus *mpStimulus;
public:
    PointStimulusHeartCellFactoryMetis(double timeStep) : AbstractCardiacCellFactory<3>(timeStep)
    {
        mpStimulus = new InitialStimulus(-1000.0*1000, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus, mpZeroStimulus);
    }
    
    void FinaliseCellCreation(std::vector<AbstractCardiacCell* >* pCellsDistributed, unsigned lo, unsigned hi)
    {
/* 
 * Here's the list of stimulated cells from the heart mesh file with tetgen numbering:
37483   1.95075 0.02458 0.007709
37498   1.974   0.0669055       0.0212167
37776   1.92132 -0.0185282      0.0264612
37778   1.90362 -0.0457586      0.0653502
38007   1.93232 -0.0006106      0.0718023
38331   1.9199  -0.0110326      -0.0303119
38586   1.90586 -0.0489975      -0.0275832
38587   1.90054 -0.0704444      0.0187846
39311   1.95105 0.0306952       -0.0127931
39313   1.97209 0.072277        -0.0302588
39642   1.93571 0.0272909       -0.0672191
40587   1.95069 0.0286633       -0.0049338
40589   1.97168 0.0738751       -0.0122153
63884   1.93084 0       0

 * Here's the list of stimulated cells from the halfheart mesh file with tetgen numbering:
37483   0.975375        0.01229 0.0038545
37498   0.987   0.0334528       0.0106084
37776   0.96066 -0.0092641      0.0132306
37778   0.95181 -0.0228793      0.0326751
38007   0.96616 -0.0003053      0.0359011
38331   0.95995 -0.0055163      -0.0151559
38586   0.95293 -0.0244987      -0.0137916
38587   0.95027 -0.0352222      0.0093923
39311   0.975525        0.0153476       -0.00639655
39313   0.986045        0.0361385       -0.0151294
39642   0.967855        0.0136454       -0.0336096
40587   0.975345        0.0143316       -0.0024669
40589   0.98584 0.0369375       -0.00610765
63884   0.96542 0       0
* Here's the list of stimulated cells from the new metis mesh file with tetgen numbering:
39934   0.97537499999999999201  0.012290000000000000577 0.0038544999999999998534        0
40021   0.98699999999999998845  0.033452799999999997815 0.010608400000000000468 0
40650   0.9606599999999999584   -0.0092641000000000008063       0.013230600000000000346 0
39295   0.95181000000000004491  -0.022879300000000001719        0.03267509999999999859  0
40639   0.96616000000000001879  -0.00030529999999999999413      0.035901099999999998291 0
40700   0.95994999999999996998  -0.0055163000000000000228       -0.015155899999999999928    0
39299   0.95293000000000005478  -0.024498700000000001725        -0.013791599999999999346    0
39360   0.95026999999999994806  -0.035222200000000002118        0.0093922999999999992604    0
39925   0.97552499999999997549  0.015347599999999999437 -0.0063965499999999999789       0
39900   0.9860449999999999493   0.036138499999999997014 -0.015129399999999999446        0
40620   0.96785500000000002085  0.013645400000000000237 -0.033609600000000003361        0
39991   0.97534500000000001751  0.014331599999999999895 -0.0024669000000000001691       0
40020   0.9858400000000000496   0.036937499999999998113 -0.0061076500000000000581       0
39993   0.96541999999999994486  0       0       0
*/
        int stimulated_cells[] = {  39295, 
                                    39299,  
                                    39360,
                                    39900,
                                    39925,
                                    39934,
                                    39991,
                                    39993,
                                    40020,
                                    40021,
                                    40620,
                                    40639,
                                    40650,
                                    40700
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
    
    ~PointStimulusHeartCellFactoryMetis(void)
    {
        delete mpStimulus;
    }
};

class TestBidomainHeart : public CxxTest::TestSuite
{

public:

    /*
     * 
    void xTestPermuteWithMetisBinaries()
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/halfheart");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        mesh.PermuteNodesWithMetisBinaries();
        
        TrianglesMeshWriter<3,3> mesh_writer("","halfheart_metis");
        mesh_writer.WriteFilesUsingMesh(mesh);
    }
    *
    */
    
    void TestBidomainDg0Heart() throw (Exception)
    {
        double pde_time_step = 0.005;  // ms
        double ode_time_step = 0.0025; // ms
        double end_time = 100;        // ms
        double printing_time_step = end_time/1000;
        
        PointStimulusHeartCellFactory cell_factory(ode_time_step);
        BidomainProblem<3> bidomain_problem(&cell_factory);
        
        bidomain_problem.SetMeshFilename("heart/test/data/halfheart");
        bidomain_problem.SetOutputDirectory("BiDg0Heart");
        bidomain_problem.SetOutputFilenamePrefix("BidomainLR91_Heart");
        
        bidomain_problem.SetEndTime(end_time);
        bidomain_problem.SetPdeTimeStep(pde_time_step);
        bidomain_problem.SetPrintingTimeStep(printing_time_step);
        
        bidomain_problem.SetLinearSolverRelativeTolerance(5e-7);
        PetscOptionsSetValue("-ksp_type", "symmlq");
        PetscOptionsSetValue("-pc_type", "bjacobi");
        //PetscOptionsSetValue("-log_summary", "");
        //PetscOptionsSetValue("-ksp_monitor", "");
        PetscOptionsSetValue("-options_table", "");
        
        bidomain_problem.SetWriteInfo();

        bidomain_problem.SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        bidomain_problem.SetExtracellularConductivities(Create_c_vector(7.0, 7.0, 7.0));
        
        bidomain_problem.Initialise();
        bidomain_problem.Solve();
    }
    
       void xTestBidomainDg0HeartMetis() throw (Exception)
    {
        double pde_time_step = 0.005;  // ms
        double ode_time_step = 0.0025; // ms
        double end_time = 100;        // ms
        double printing_time_step = end_time/1000;
        
        PointStimulusHeartCellFactoryMetis cell_factory(ode_time_step);
        BidomainProblem<3> bidomain_problem(&cell_factory);
        
        bidomain_problem.SetMeshFilename("heart/test/data/halfheart_metis");
        bidomain_problem.SetOutputDirectory("BiDg0HeartMetis");
        bidomain_problem.SetOutputFilenamePrefix("BidomainLR91_HeartMetis");
        
        bidomain_problem.SetEndTime(end_time);
        bidomain_problem.SetPdeTimeStep(pde_time_step);
        bidomain_problem.SetPrintingTimeStep(printing_time_step);
        
        //bidomain_problem.SetLinearSolverRelativeTolerance(5e-7);
        //PetscOptionsSetValue("-ksp_type", "symmlq");
        //PetscOptionsSetValue("-pc_type", "bjacobi");
        //PetscOptionsSetValue("-log_summary", "");
        //PetscOptionsSetValue("-ksp_monitor", "");
        PetscOptionsSetValue("-options_table", "");
        
        bidomain_problem.SetWriteInfo();
        
        bidomain_problem.Initialise();
        bidomain_problem.Solve();
    }
};

#endif //_TESTBIDOMAINHEART_HPP_
