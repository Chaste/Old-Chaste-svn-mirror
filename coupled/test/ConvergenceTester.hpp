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
//    double kspRtol;
//    double conductionVelocity;
    double spaceStep;
//    double stimulusMagnitude;??

public:    
    ConvergenceTester()
    {
        const unsigned number_of_meshes = 5;
        
        double num_elements[number_of_meshes];
        std::string file_name[number_of_meshes];
        
        unsigned opposite_corner_node[number_of_meshes];
        
        double space_steps[number_of_meshes];
        
        // Create the meshes on which the test will be based
        const std::string mesh_dir = "BidomainConvergenceMesh";
        OutputFileHandler output_file_handler(mesh_dir);
        
        for (unsigned i = 0; i < number_of_meshes; i++)
        {
            // start with 4 elements along each edge
            unsigned mesh_size = (unsigned) pow(2, i+2);
            double scaling = mesh_width/(double) mesh_size;
            
            ConformingTetrahedralMesh<DIM,DIM> mesh;
            
            ConstructHyperCube(mesh, mesh_size);
            mesh.Scale(scaling, scaling, scaling);
             
            num_elements[i] = mesh.GetNumElements();
            opposite_corner_node[i] = mesh.GetNumNodes()-1;
            
//            for(unsigned dim=0; dim<DIM; dim++)
//            {
//                TS_ASSERT_DELTA(mesh.GetNode(opposite_corner_node[i])->rGetLocation()[dim],mesh_width,1e-6);
//            }
            
            space_steps[i] = scaling;
            
            std::stringstream file_name_stream;
            file_name_stream<< "cube_" << DIM << "D_2mm_"<< num_elements[i] <<"_elements";
            file_name[i]=file_name_stream.str();
            
            TrianglesMeshWriter<DIM,DIM> mesh_writer(mesh_dir, file_name[i], false);
            
            mesh_writer.WriteFilesUsingMesh(mesh);
        }
        
        // To ensure that the first test fails
        double prev_velocity_for_space = -999;
        bool converging_in_space = false;
        bool failed_to_converge_in_space = false;
        
        double pde_time_step;   // ms
        
        double probe_velocity;
        ReplicatableVector voltage_replicated;
        
        unsigned current_file_num = 0;
        double ode_time_step;

        do //do while: space_step
        {
            bool converging_in_time = false;
            // To ensure that the first test fails
            double prev_velocity_for_time = -999;
            
            pde_time_step = initial_pde_time_step;  // ms
            //Following line is redundant
            ode_time_step=pde_time_step;
            
            std::string mesh_pathname = output_file_handler.GetTestOutputDirectory()
                                        + file_name[current_file_num];
                                        
            std::cout<<"================================================================================"<<std::endl  << std::flush;
            std::cout<<"Solving with a space step of "<< space_steps[current_file_num] << " cm - mesh " << current_file_num <<std::endl  << std::flush;
            
            do  //do while: pde_time_step
            {

                /*******/
                //Over-ride the new ode_time_step tweaking code so that we might expect rapid 
                //pde_time_step convergence
                ode_time_step=0.0025;
                /*******/
                
                PointStimulusCellFactory<CELL, DIM> cell_factory(ode_time_step, num_elements[current_file_num]);
                CARDIAC_PROBLEM cardiac_problem(&cell_factory);
                
                cardiac_problem.SetMeshFilename(mesh_pathname);
                cardiac_problem.SetOutputDirectory ("Convergence");
                cardiac_problem.SetOutputFilenamePrefix ("Results");
                
                
                cardiac_problem.SetEndTime(simulation_time);   // ms
                
                cardiac_problem.SetLinearSolverRelativeTolerance(1e-6);

                cardiac_problem.SetPdeTimeStep(pde_time_step);
                cardiac_problem.SetPrintingTimeStep(initial_pde_time_step);
                cardiac_problem.Initialise();
                
                std::cout<<"   Solving with a time step of "<<pde_time_step<<" ms"<<std::endl  << std::flush;
                std::cout<<"   Solving with an ode time step of "<<ode_time_step<<" ms"<<std::endl  << std::flush;
                
//                try
                {
                    cardiac_problem.Solve();
//                    Vec voltage=cardiac_problem.GetVoltage();
//                    voltage_replicated.ReplicatePetscVector(voltage);
                    

                    
                    // Calculate conduction velocity between 1/4 and 3/4 through the mesh
                    
                    unsigned third_quadrant_node;
                    unsigned first_quadrant_node;
                    
                    switch(DIM)
                    {
                        case 1:
                        {
                            first_quadrant_node = (unsigned) (0.25*num_elements[current_file_num]);
                            third_quadrant_node = (unsigned) (0.75*num_elements[current_file_num]);
                            assert(cardiac_problem.rGetMesh().GetNode(first_quadrant_node)->rGetLocation()[0]==0.25*mesh_width);
                            assert(cardiac_problem.rGetMesh().GetNode(third_quadrant_node)->rGetLocation()[0]==0.75*mesh_width);
                            break;
                        }
                        case 2:
                        {
                            unsigned n= (unsigned) pow (2, current_file_num+2);
                            first_quadrant_node =   (n+1)*(n/2)+  n/4 ;
                            third_quadrant_node =   (n+1)*(n/2)+3*n/4 ;
                            break;
                        }
                        case 3:
                        {
                            int n= (int) pow (2, current_file_num+2);
                            first_quadrant_node =   61*(n-8)*(n-16)/48 + 362*(n-4)*(n-16)/-32 + 2452*(n-4)*(n-8)/96;
                            third_quadrant_node =   63*(n-8)*(n-16)/48 + 366*(n-4)*(n-16)/-32 + 2460*(n-4)*(n-8)/96;
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
//                catch (Exception e)
//                {
//                    // An exception has been caught, meaning that the time step is too big, so halve it
//                    std::cout << "   >>> Convergence test: an exception was thrown (" << e.GetMessage() << ")" << std::endl  << std::flush;
//
//                    
//                    // try halving ode_time_step 
////                    if (ode_time_step> pde_time_step/33.0)
////                    {
////                        /******* Currently redundant code*/
////                        ode_time_step *= 0.5;
////                        std::cout << "   >>>                       We assume that the ode time step was too big" << std::endl << std::flush;
////                    }
////                    // unless its very small in which case reduce pde time step 
////                    else
//                    {
//                        std::cout << "   >>>                   We assume that the pde time step was too big" << std::endl << std::flush;
//                        pde_time_step *= 0.5;
//                        ode_time_step = pde_time_step;
//                    }
//                }
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
                current_file_num++;
                if (current_file_num==number_of_meshes)
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
            
            std::cout << "Converged both in space ("<< space_steps[current_file_num] <<" cm) and time ("<< pde_time_step << " ms)" << std::endl << std::flush;
        }
        
        // TS_ASSERT_DELTA(space_steps[current_file_num], 0.005, 0.0);
        // TS_ASSERT_DELTA(pde_time_step, 0.005, 0.0);
        // TS_ASSERT_DELTA(probe_voltage, -10.3432, 0.0001);
        // Note: the delta is because of floating point issues (!!)
        odeTimeStep=ode_time_step;
        pdeTimeStep=pde_time_step;
        spaceStep=space_steps[current_file_num];
        convergedInSpace=true;
    }    
};
#endif /*CONVERGENCETESTER_HPP_*/
