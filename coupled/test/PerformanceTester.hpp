#ifndef ABSTRACTCONVERGENCETESTER_HPP_
#define ABSTRACTCONVERGENCETESTER_HPP_

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
#include "ColumnDataReader.hpp"

const double mesh_width = 0.2; // cm
// meshes to use for convergence testing.
// n-dimensional cubes with 2^(mesh_num+2) elements in each dimension


template <class CELL, unsigned DIM>
class PointStimulusCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    // define a new stimulus
    InitialStimulus* mpStimulus;
    
public:
    PointStimulusCellFactory(double timeStep, double numElements) : AbstractCardiacCellFactory<DIM>(timeStep)
    {
        // scale stimulus depending on space_step of elements
        //\todo It looks like the value of the stimulus is specific to 3D

        switch(DIM)
        {
            case 1:
            {
                mpStimulus = new InitialStimulus(-10000000*numElements/64.0, 0.5);
                break;
            }
            case 2:
            {
                mpStimulus = new InitialStimulus(-5000*numElements, 0.5);
                break;
            }
            case 3:
            {
                assert(0);
                break;
            }
            default:
                assert(0);
        }
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        double x = this->mpMesh->GetNode(node)->GetPoint()[0];
        if (x*x<=1e-10)
        {
            return new CELL(this->mpSolver, this->mTimeStep, this->mpStimulus, this->mpZeroStimulus);
        }
        else
        {
            return new CELL(this->mpSolver, this->mTimeStep, this->mpZeroStimulus, this->mpZeroStimulus);
        }
    }
    
    ~PointStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};


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

        double scaling;
        unsigned prev_mesh_num=9999;
        std::string mesh_pathname;
        std::string mesh_filename;


        if (MeshNum!=prev_mesh_num)
        {
            // create the mesh
            unsigned mesh_size = (unsigned) pow(2, MeshNum+2); // number of elements in each dimension
            scaling = mesh_width/(double) mesh_size;
            ConformingTetrahedralMesh<DIM,DIM> mesh;
            ConstructHyperCube(mesh, mesh_size);
            mesh.Scale(scaling, scaling, scaling);
            mNumElements = mesh.GetNumElements();
            mNumNodes = mesh.GetNumNodes();
            std::stringstream file_name_stream;
            file_name_stream<< "cube_" << DIM << "D_2mm_"<< mNumElements <<"_elements";
            mesh_filename = file_name_stream.str();
            TrianglesMeshWriter<DIM,DIM> mesh_writer(mesh_dir, mesh_filename, false);           
            mesh_writer.WriteFilesUsingMesh(mesh);
            mesh_pathname = output_file_handler.GetTestOutputDirectory()
                                        + mesh_filename;
            prev_mesh_num = MeshNum;
        }                            
        
        PointStimulusCellFactory<CELL, DIM> cell_factory(OdeTimeStep, mNumElements);
        CARDIAC_PROBLEM cardiac_problem(&cell_factory);
        
        cardiac_problem.SetMeshFilename(mesh_pathname);
        cardiac_problem.SetOutputDirectory ("Convergence");
        cardiac_problem.SetOutputFilenamePrefix ("Results");
        
        cardiac_problem.SetEndTime(SimTime);   // ms
        cardiac_problem.SetLinearSolverRelativeTolerance(KspRtol);

        cardiac_problem.SetPdeTimeStep(PdeTimeStep);
        
        assert(fabs(0.04/PdeTimeStep - round(0.04/PdeTimeStep)) <1e-15 );
        cardiac_problem.SetPrintingTimeStep(PrintingTimeStep);  //Otherwise we can't take the timestep down to machine precision without generating thousands of output files
        cardiac_problem.Initialise();


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
    
//    virtual ~AbstractConvergenceTester() {}

};
#endif /*ABSTRACTCONVERGENCETESTER_HPP_*/
