#ifndef CONVERGENCETESTER_HPP_
#define CONVERGENCETESTER_HPP_

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
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "BackwardEulerLuoRudyIModel1991.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshWriter.cpp"
#include "PropagationPropertiesCalculator.hpp"
#include "ColumnDataReader.hpp"

const double initial_pde_time_step = 0.04; //ms
const double mesh_width = 0.2; // cm
const double simulation_time = 8.0; //ms
// meshes to use for convergence testing.
// n-dimensional cubes with 2^(mesh_num+2) elements in each dimension 
const unsigned first_mesh=0;
const unsigned last_mesh=4;

const unsigned first_quadrant_nodes_3d[5]={61, 362, 2452, 17960, 137296};
const unsigned third_quadrant_nodes_3d[5]={63, 366, 2460, 17976, 137328};


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
        mpStimulus = new InitialStimulus(-1000000, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        double x = this->mpMesh->GetNode(node)->GetPoint()[0];
        //double y = mpMesh->GetNode(node)->GetPoint()[1];
        //double z = mpMesh->GetNode(node)->GetPoint()[2];
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
class ConvergenceTester
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
    bool convergedInSpace;
    double odeTimeStep;
    double pdeTimeStep;
    double spaceStep;

public:    
    ConvergenceTester()
    {
        
        // Create the meshes on which the test will be based
        const std::string mesh_dir = "ConvergenceMesh";
        OutputFileHandler output_file_handler(mesh_dir);
        
        // To ensure that the first test fails
        double prev_velocity_for_space = -999;
        bool converging_in_space = false;
        bool failed_to_converge_in_space = false;
        
        double pde_time_step;   // ms
        
        double probe_velocity;
        ReplicatableVector voltage_replicated;
        
        unsigned mesh_num = first_mesh;
        double ode_time_step;
        double scaling;

        do //do while: space_step
        {
            
            // create the mesh
            unsigned mesh_size = (unsigned) pow(2, mesh_num+2); // number of elements in each dimension
            scaling = mesh_width/(double) mesh_size;
            ConformingTetrahedralMesh<DIM,DIM> mesh;
            ConstructHyperCube(mesh, mesh_size);
            mesh.Scale(scaling, scaling, scaling);
            unsigned num_elements = mesh.GetNumElements();
            std::stringstream file_name_stream;
            file_name_stream<< "cube_" << DIM << "D_2mm_"<< num_elements <<"_elements";
            std::string file_name = file_name_stream.str();
            TrianglesMeshWriter<DIM,DIM> mesh_writer(mesh_dir, file_name, false);           
            mesh_writer.WriteFilesUsingMesh(mesh);
            
            bool converging_in_time = false;
            // To ensure that the first test fails
            double prev_velocity_for_time = -999;
            pde_time_step = initial_pde_time_step;  // ms
            //Following line is redundant
            ode_time_step=pde_time_step;
            
            std::string mesh_pathname = output_file_handler.GetTestOutputDirectory()
                                        + file_name;
                                        
            std::cout<<"================================================================================"<<std::endl  << std::flush;
            std::cout<<"Solving with a space step of "<< scaling << " cm - mesh " << mesh_num <<std::endl  << std::flush;
            
            do  //do while: pde_time_step
            {

                /*******/
                //Over-ride the new ode_time_step tweaking code so that we might expect rapid 
                //pde_time_step convergence
                ode_time_step=0.0025;
                /*******/
                
                PointStimulusCellFactory<CELL, DIM> cell_factory(ode_time_step, num_elements);
                CARDIAC_PROBLEM cardiac_problem(&cell_factory);
                
                cardiac_problem.SetMeshFilename(mesh_pathname);
                cardiac_problem.SetOutputDirectory ("Convergence");
                cardiac_problem.SetOutputFilenamePrefix ("Results");
                
                cardiac_problem.SetEndTime(simulation_time);   // ms
                cardiac_problem.SetLinearSolverRelativeTolerance(1e-6);

                cardiac_problem.SetPdeTimeStep(pde_time_step);
                cardiac_problem.SetPrintingTimeStep(pde_time_step);
                cardiac_problem.Initialise();
                
                std::cout<<"   Solving with a time step of "<<pde_time_step<<" ms"<<std::endl  << std::flush;
                std::cout<<"   Solving with an ode time step of "<<ode_time_step<<" ms"<<std::endl  << std::flush;
                
                cardiac_problem.Solve();
                
                // Calculate conduction velocity between 1/4 and 3/4 through the mesh
                unsigned third_quadrant_node;
                unsigned first_quadrant_node;
                switch(DIM)
                {
                    case 1:
                    {
                        first_quadrant_node = (unsigned) (0.25*num_elements);
                        third_quadrant_node = (unsigned) (0.75*num_elements);
                        assert(cardiac_problem.rGetMesh().GetNode(first_quadrant_node)->rGetLocation()[0]==0.25*mesh_width);
                        assert(cardiac_problem.rGetMesh().GetNode(third_quadrant_node)->rGetLocation()[0]==0.75*mesh_width);
                        break;
                    }
                    case 2:
                    {
                        unsigned n= (unsigned) pow (2, mesh_num+2);
                        first_quadrant_node =   (n+1)*(n/2)+  n/4 ;
                        third_quadrant_node =   (n+1)*(n/2)+3*n/4 ;
                        break;
                    }
                    case 3:
                    {
                        first_quadrant_node = first_quadrant_nodes_3d[mesh_num];
                        third_quadrant_node = third_quadrant_nodes_3d[mesh_num];
                        break;
                    }
                    
                    default:
                        assert(0);
                }
                
                Node<DIM>* fqn = cardiac_problem.rGetMesh().GetNode(first_quadrant_node);
                Node<DIM>* tqn = cardiac_problem.rGetMesh().GetNode(third_quadrant_node);
                assert(fqn->rGetLocation()[0]==0.25*mesh_width);
                assert(tqn->rGetLocation()[0]==0.75*mesh_width);
                for (unsigned coord=1; coord<DIM; coord++)
                {
                    assert(fqn->rGetLocation()[coord]==0.5*mesh_width);
                    assert(tqn->rGetLocation()[coord]==0.5*mesh_width);
                }
                
                
                OutputFileHandler results_handler("Convergence", false);
                ColumnDataReader results_reader(results_handler.GetTestOutputDirectory(), "Results", false);
                PropagationPropertiesCalculator propagation_calc(&results_reader);
                probe_velocity = propagation_calc.CalculateConductionVelocity(first_quadrant_node, third_quadrant_node, 0.5*mesh_width);
                
                double relerr = fabs ((probe_velocity - prev_velocity_for_time) / prev_velocity_for_time);
                std::cout<<"   >>> Convergence test: conduction velocity= "<<probe_velocity<<" cm/msec | prev_velocity_for_time = "<<prev_velocity_for_time
                <<" cm/msec | relerr = "<<relerr<<std::endl  << std::flush;
                
                if (relerr < 1e-2)
                {
                    converging_in_time = true;
                }
                else
                {
                    // Get ready for the next test by halving the time step
                    pde_time_step *= 0.5;
                    /******* Currently redundant code*/
                    ode_time_step = pde_time_step;
                }
                
                if (pde_time_step < 1e-4)
                {
                    std::cout << "**** NO TIMESTEP GREATER THAN 1e-4 FOUND WHICH WORKS, MOVING ONTO NEXT MESH...****\n" << std::flush;
                    converging_in_time = true;
                }
                prev_velocity_for_time = probe_velocity;
                
            }
            while (!converging_in_time);   //do while: pde_time_step
            
            double relerr = fabs ((probe_velocity - prev_velocity_for_space) / prev_velocity_for_space);
            std::cout<<">>> Convergence test: conduction velocity = "<<probe_velocity<<" cm/msec | prev_velocity_for_space = "<<prev_velocity_for_space
            <<" cm/msec | relerr = "<<relerr<<std::endl << std::flush;
            
            
            if (relerr < 1e-2)
            {
                converging_in_space = true;
            }
            else
            {
                // Use the next mesh next time
                mesh_num++;
                if (mesh_num>=last_mesh)
                {
                    convergedInSpace=false;
                    
                    failed_to_converge_in_space = true;
                    return;
                }
            }
            prev_velocity_for_space = probe_velocity;
        }
        while (!converging_in_space && !failed_to_converge_in_space);   //do while: space_step
        
        if (converging_in_space)
        {
            std::cout<<"================================================================================"<<std::endl << std::flush;
            
            std::cout << "Converged both in space ("<< scaling <<" cm) and time ("<< pde_time_step << " ms)" << std::endl << std::flush;
        }
        
        TS_ASSERT(converging_in_space);
        
        odeTimeStep=ode_time_step;
        pdeTimeStep=pde_time_step;
        convergedInSpace=true;
    }    
};
#endif /*CONVERGENCETESTER_HPP_*/
