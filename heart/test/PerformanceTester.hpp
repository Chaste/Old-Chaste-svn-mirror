#ifndef PERFORMANCETESTER_HPP_
#define PERFORMANCETESTER_HPP_

#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "ConformingTetrahedralMesh.cpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshWriter.cpp"
#include "PropagationPropertiesCalculator.hpp"
#include "Hdf5DataReader.hpp"
#include "GeneralPlaneStimulusCellFactory.hpp"
#include "CuboidMeshConstructor.hpp"

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM>
class PerformanceTester
{
private:
    void ConstructHyperCube(ConformingTetrahedralMesh<1,1> &rMesh, unsigned width)
    {
        rMesh.ConstructLinearMesh(width);
    }
    void ConstructHyperCube(ConformingTetrahedralMesh<2,2> &rMesh, unsigned width)
    {
        rMesh.ConstructRectangularMesh(width, width);
    }
    void ConstructHyperCube(ConformingTetrahedralMesh<3,3> &rMesh, unsigned width)
    {
        rMesh.ConstructCuboid(width, width, width);
    }
protected:
    const static double mMeshWidth=0.2;//cm
public:    
    PerformanceTester()
    : OdeTimeStep(0.0025),
      PdeTimeStep(0.0025),
      MeshNum(2u),
      KspRtol(1e-8),
      RelativeConvergenceCriterion(1e-4),
      SimTime(8.0),
      PrintingTimeStep(0.04)
    {
    }
    
    void Run()
    {
        // Create the meshes on which the test will be based
        const std::string mesh_dir = "ConvergenceMesh";
        OutputFileHandler output_file_handler(mesh_dir);

        unsigned prev_mesh_num=9999;
        std::string mesh_pathname;
        std::string mesh_filename;

        CuboidMeshConstructor<DIM> constructor;

        if (MeshNum!=prev_mesh_num)
        {
            mesh_pathname = constructor.Construct(MeshNum, mMeshWidth);
            mNumElements = constructor.NumElements;
            mNumNodes = constructor.NumNodes;
            prev_mesh_num = MeshNum;
        }
        
        unsigned num_ele_across = (unsigned) pow(2, this->MeshNum+2);
        GeneralPlaneStimulusCellFactory<CELL, DIM> cell_factory(OdeTimeStep, num_ele_across, mMeshWidth);
        CARDIAC_PROBLEM cardiac_problem(&cell_factory);
        
        cardiac_problem.SetMeshFilename(mesh_pathname);
        cardiac_problem.SetOutputDirectory ("Convergence");
        cardiac_problem.SetOutputFilenamePrefix ("Results");
        
        cardiac_problem.SetEndTime(SimTime);   // ms
        cardiac_problem.SetLinearSolverRelativeTolerance(KspRtol);

        assert(fabs(0.04/PdeTimeStep - round(0.04/PdeTimeStep)) <1e-15 );
        cardiac_problem.SetPdeAndPrintingTimeSteps(PdeTimeStep, PrintingTimeStep);
        cardiac_problem.Initialise();

        //// use this to get some info printed out
        //cardiac_problem.SetWriteInfo();
        
        try
        {
            cardiac_problem.Solve();
        }
        catch (Exception e)
        {
            std::cout<<"Warning - this run threw an exception.  Check convergence results\n";  
            throw(e);               
        }
        
        DisplayRun();
    }
    
    static void DisplayHeadings()
    {
        const unsigned NUM_HEADINGS=7;
        const char* heading[NUM_HEADINGS]={"Dimen", "Elts", "Nodes", "PdeStp", "OdeStp", "PriStp", "SimTim"} ;
        for (unsigned i=0; i<NUM_HEADINGS; i++)
        {
            printf("%6s\t", heading[i]);
        }
    }
    
    void DisplayRun()
    {
        printf("%6u\t%6u\t%6u\t%2.1e\t%2.1e\t%2.1e\t%2.1e\t",
               DIM, mNumElements, mNumNodes, PdeTimeStep, OdeTimeStep, PrintingTimeStep, SimTime);
    }
    

public:
    double OdeTimeStep;
    double PdeTimeStep;
    unsigned MeshNum;
    double KspRtol;
    double RelativeConvergenceCriterion;
    double SimTime;
    double PrintingTimeStep;
    
private:
    unsigned mNumNodes;
    unsigned mNumElements;

};
#endif /*PERFORMANCETESTER_HPP_*/
