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
#include "OutputFileHandler.hpp"
#include "TrianglesMeshWriter.cpp"
#include "PropagationPropertiesCalculator.hpp"
#include "ColumnDataReader.hpp"

const double initial_pde_time_step = 0.04; //ms
const double mesh_width = 0.2; // cm
const double simulation_time = 8.0; //ms
// meshes to use for convergence testing.
// n-dimensional cubes with 2^(mesh_num+2) elements in each dimension 
const unsigned first_mesh=4;
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
        mpStimulus = new InitialStimulus(-10000000*numElements/64.0, 0.5);
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
        bool converging_in_space = false;
        double pde_time_step;   // ms
        ReplicatableVector voltage_replicated;

        unsigned file_num=0;
        double ksp_rtol=1e-8;
        double scaling;
        unsigned mesh_num = first_mesh;
        double ode_time_step;
        bool converged=false;
            
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
            
            pde_time_step = initial_pde_time_step;  // ms
            //Following line is redundant
            ode_time_step=pde_time_step;
            
            std::string mesh_pathname = output_file_handler.GetTestOutputDirectory()
                                        + file_name;
                                        
            std::cout<<"================================================================================"<<std::endl  << std::flush;
            std::cout<<"Solving with a space step of "<< scaling << " cm - mesh " << mesh_num <<std::endl  << std::flush;
            
            
            unsigned time_step_write_increment=1;
            ode_time_step = pde_time_step;
            double prev_voltage[201];

            do  //do while: pde_time_step
            {
                ode_time_step=pde_time_step;
                
                PointStimulusCellFactory<CELL, DIM> cell_factory(ode_time_step, num_elements);
                CARDIAC_PROBLEM cardiac_problem(&cell_factory);
                
                cardiac_problem.SetMeshFilename(mesh_pathname);
                cardiac_problem.SetOutputDirectory ("Convergence");
                cardiac_problem.SetOutputFilenamePrefix ("Results");
                
                cardiac_problem.SetEndTime(simulation_time);   // ms
                cardiac_problem.SetLinearSolverRelativeTolerance(ksp_rtol);

                cardiac_problem.SetPdeTimeStep(pde_time_step);
                cardiac_problem.SetPrintingTimeStep(pde_time_step);
                cardiac_problem.Initialise();
                
                std::cout<<"   Solving with a time step of "<<pde_time_step<<" ms"<<std::endl  << std::flush;
                std::cout<<"   Solving with an ode time step of "<<ode_time_step<<" ms"<<std::endl  << std::flush;

                cardiac_problem.Solve();
                // Calculate positions of nodes 1/4 and 3/4 through the mesh
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
                
                // Write out the time series for the node at first and third quadrant
                {
                    std::vector<double> transmembrane_potential=results_reader.GetValues("V", third_quadrant_node);
                    std::vector<double> time_series = results_reader.GetUnlimitedDimensionValues();
                    OutputFileHandler plot_file_handler("ConvergencePlots", false);
                    std::stringstream plot_file_name_stream;
                    plot_file_name_stream<< "Node1_"<< file_num << "_timestep.csv";
                    out_stream p_plot_file = plot_file_handler.OpenOutputFile(plot_file_name_stream.str());
                    for (unsigned data_point = 0; data_point<time_series.size(); data_point+=time_step_write_increment)
                    {
                        (*p_plot_file) << time_series[data_point] << "\t" << transmembrane_potential[data_point] << "\n";                 
                    }
                    p_plot_file->close();
                    // calculate l2norm
                    double *p_prev_voltage = prev_voltage;
                    double max_abs_error = 0;
                    double sum_sq_abs_error =0;
                    double sum_sq_prev_voltage = 0;
                    for (unsigned data_point = 0; data_point<time_series.size(); data_point+=time_step_write_increment)
                    {
                        if (file_num!=0)
                        {
                            double abs_error = fabs(transmembrane_potential[data_point]-*p_prev_voltage);
                            max_abs_error = (abs_error > max_abs_error) ? abs_error : max_abs_error;
                            sum_sq_abs_error += abs_error*abs_error;
                            sum_sq_prev_voltage += *p_prev_voltage * *p_prev_voltage;
                        } 
                        
                        *p_prev_voltage= transmembrane_potential[data_point];
                        p_prev_voltage++;
                    }                 
                    if (file_num!=0)
                    {
                        std::cout << "log10 timestep= " << log10(pde_time_step) << "\n";
                        std::cout << "max_abs_error = " << max_abs_error << " log10 = " << log10(max_abs_error) << "\n";
                        std::cout << "l2 error = " << sum_sq_abs_error/sum_sq_prev_voltage << " log10 = " << log10(sum_sq_abs_error/sum_sq_prev_voltage) << "\n";
                        std::cout << log10(pde_time_step) << "\t"<<log10(max_abs_error)<<"\t"
                            <<log10(sum_sq_abs_error/sum_sq_prev_voltage) <<"\t#Logs for Gnuplot\n";
                            
                        // convergence criterion
                        converged = sum_sq_abs_error/sum_sq_prev_voltage<1e-4;
                    }
                }
                
                {
                    std::vector<double> transmembrane_potential=results_reader.GetValues("V", first_quadrant_node);
                    std::vector<double> time_series = results_reader.GetUnlimitedDimensionValues();
                    OutputFileHandler plot_file_handler("ConvergencePlots", false);
                    std::stringstream plot_file_name_stream;
                    plot_file_name_stream<< "Node2_"<< file_num << "_timestep.csv";
                    out_stream p_plot_file = plot_file_handler.OpenOutputFile(plot_file_name_stream.str());
                    for (unsigned data_point = 0; data_point<time_series.size(); data_point+=time_step_write_increment)
                    {
                        (*p_plot_file) << time_series[data_point] << "\t" << transmembrane_potential[data_point] << "\n";                 
                    }
                    p_plot_file->close();
                }                    
                
                // Get ready for the next test by halving the time step
                if (!converged)
                {
                    pde_time_step *= 0.5;
                    time_step_write_increment *= 2;
                    file_num++;
                }
            }
            while (pde_time_step> 1e-8 && !converged);   //do while: pde_time_step
            mesh_num++;
        }
        while (!converging_in_space && mesh_num<=last_mesh);   //do while: space_step
        
        TS_ASSERT_DELTA(pde_time_step, 2.5e-3, 1e-10);
    }    
};
#endif /*CONVERGENCETESTER_HPP_*/
