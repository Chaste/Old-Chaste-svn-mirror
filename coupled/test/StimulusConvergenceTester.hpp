#ifndef STIMULUSCONVERGENCETESTER_HPP_
#define STIMULUSCONVERGENCETESTER_HPP_


#include "AbstractConvergenceTester.hpp"

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM>
class StimulusConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM>
{
public:
    double ExpectedStimulus;
    double Increment;
    unsigned FirstMesh;
    unsigned PlotPoints;
    void PopulateStandardResult(double prev_voltage[])
    {
        assert(this->PopulatedResult==false);
  
        this->MeshNum=FirstMesh;
        CuboidMeshConstructor<DIM> constructor;
        std::string mesh_pathname;
 
        mesh_pathname = constructor.Construct(this->MeshNum);
        
        
        assert (this->UseAbsoluteStimulus==true);
        GeneralPlaneStimulusCellFactory<CELL, DIM> cell_factory(this->OdeTimeStep, this->AbsoluteStimulus, true);
        CARDIAC_PROBLEM cardiac_problem(&cell_factory);
        
        cardiac_problem.SetMeshFilename(mesh_pathname);
        cardiac_problem.SetOutputDirectory ("Convergence");
        cardiac_problem.SetOutputFilenamePrefix ("Results");
        
        cardiac_problem.SetEndTime(simulation_time);   // ms
        cardiac_problem.SetLinearSolverRelativeTolerance(this->KspRtol);

        cardiac_problem.SetPdeTimeStep(this->PdeTimeStep);
        
        assert(fabs(0.04/this->PdeTimeStep - round(0.04/this->PdeTimeStep)) <1e-15 );
        cardiac_problem.SetPrintingTimeStep(0.04);  //Otherwise we can't take the timestep down to machine precision without generating thousands of output files
        cardiac_problem.Initialise();
        
        this->DisplayRun();
        try
        {
            cardiac_problem.Solve();
        }
        catch (Exception e)
        {
            std::cout<<"Warning - this run threw an exception.  Check convergence results\n"; 
            assert(0); //Can't continue                
        }
        // Calculate positions of nodes 1/4 and 3/4 through the mesh
        unsigned third_quadrant_node;
        unsigned first_quadrant_node;
        switch(DIM)
        {
            case 1:
            {
                first_quadrant_node = (unsigned) (0.25*constructor.NumElements);
                third_quadrant_node = (unsigned) (0.75*constructor.NumElements);
                assert(cardiac_problem.rGetMesh().GetNode(first_quadrant_node)->rGetLocation()[0]==0.25*mesh_width);
                assert(cardiac_problem.rGetMesh().GetNode(third_quadrant_node)->rGetLocation()[0]==0.75*mesh_width);
                break;
            }
            case 2:
            {
                unsigned n= (unsigned) pow (2, this->MeshNum+2);
                first_quadrant_node =   (n+1)*(n/2)+  n/4 ;
                third_quadrant_node =   (n+1)*(n/2)+3*n/4 ;
                break;
            }
            case 3:
            {
                const unsigned first_quadrant_nodes_3d[5]={61, 362, 2452, 17960, 137296};
                const unsigned third_quadrant_nodes_3d[5]={63, 366, 2460, 17976, 137328};
                assert(this->MeshNum<5);
                first_quadrant_node = first_quadrant_nodes_3d[this->MeshNum];
                third_quadrant_node = third_quadrant_nodes_3d[this->MeshNum];
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
            plot_file_name_stream<< "Node1_"<< "base" << "_timestep.csv";
            out_stream p_plot_file = plot_file_handler.OpenOutputFile(plot_file_name_stream.str());
            for (unsigned data_point = 0; data_point<time_series.size(); data_point++)
            {
                (*p_plot_file) << time_series[data_point] << "\t" << transmembrane_potential[data_point] << "\n";                 
            }
            p_plot_file->close();
            
            
            for (unsigned data_point = 0; data_point<time_series.size(); data_point++)
            {
                prev_voltage[data_point] = transmembrane_potential[data_point];
            }
        }
        
        {
            std::vector<double> transmembrane_potential=results_reader.GetValues("V", first_quadrant_node);
            std::vector<double> time_series = results_reader.GetUnlimitedDimensionValues();
            OutputFileHandler plot_file_handler("ConvergencePlots", false);
            std::stringstream plot_file_name_stream;
            plot_file_name_stream<< "Node2_"<< "base" << "_timestep.csv";
            out_stream p_plot_file = plot_file_handler.OpenOutputFile(plot_file_name_stream.str());
            for (unsigned data_point = 0; data_point<time_series.size(); data_point++)
            {
                (*p_plot_file) << time_series[data_point] << "\t" << transmembrane_potential[data_point] << "\n";                 
            }
            p_plot_file->close();
        }                    
        this->PopulatedResult=true;
        
        this->MeshNum++;
        ExpectedStimulus = this->AbsoluteStimulus*2;//For next refinement
        Increment=-ExpectedStimulus/10;
        PlotPoints=8;
        this->AbsoluteStimulus = ExpectedStimulus-PlotPoints*Increment/2;
    }
    
    void SetInitialConvergenceParameters()
    {
        //this->MeshNum=1; //This is set later using FirstMesh
        this->FixedResult=true;
        this->UseAbsoluteStimulus=true;
        
        double num_elements=pow(2,FirstMesh+2);
        this->AbsoluteStimulus=-10000000*num_elements/64.0;
        
        //Real picky for now
        this->RelativeConvergenceCriterion=-1; //Never converge
    }
    
    void UpdateConvergenceParameters()
    {
        this->AbsoluteStimulus += Increment;
    }
    
    bool GiveUpConvergence()
    {
        return (this->AbsoluteStimulus>ExpectedStimulus+PlotPoints*Increment/2); //3.5e7
    }
    double Abscissa()
    {
        return this->AbsoluteStimulus;
    }
};

#endif /*STIMULUSCONVERGENCETESTER_HPP_*/
