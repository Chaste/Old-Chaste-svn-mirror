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
#include "GeneralPlaneStimulusCellFactory.hpp"
#include "CuboidMeshConstructor.hpp"
#include "OutputFileHandler.hpp"

const double simulation_time = 8.0; //ms

template <class CELL, unsigned DIM>
class QuarterStimulusCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    // define a new stimulus
    InitialStimulus* mpStimulus;
    double mMeshWidth;
public:
    QuarterStimulusCellFactory(double timeStep, double meshWidth) : AbstractCardiacCellFactory<DIM>(timeStep)
    {
        mpStimulus = new InitialStimulus(-1000000, 0.5);
        mMeshWidth=meshWidth;
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        double x = this->mpMesh->GetNode(node)->GetPoint()[0];
        if (x<=mMeshWidth*0.25+1e-10)
        {
            return new CELL(this->mpSolver, this->mTimeStep, this->mpStimulus, this->mpZeroStimulus);
        }
        else
        {
            return new CELL(this->mpSolver, this->mTimeStep, this->mpZeroStimulus, this->mpZeroStimulus);
        }
    }
    
    ~QuarterStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};


class AbstractUntemplatedConvergenceTester
{
protected:
    double mMeshWidth;   
public:
    double OdeTimeStep;
    double PdeTimeStep;
    unsigned MeshNum;
    double KspRtol;
    double RelativeConvergenceCriterion;
    double AbsoluteStimulus;
    bool PopulatedResult;
    bool FixedResult;
    bool UseAbsoluteStimulus;
    bool Converged;
    bool StimulateRegion;
    
    AbstractUntemplatedConvergenceTester()   
    : mMeshWidth(0.2),//cm
      OdeTimeStep(0.0025),//Justification from 1D test with this->PdeTimeStep held at 0.01 (allowing two hits at convergence)
      PdeTimeStep(0.005),//Justification from 1D test with this->OdeTimeStep held at 0.0025
      MeshNum(5u),//Justification from 1D test
      KspRtol(5e-7),//Justification from overlayed 1D time/space convergence plots with varied KSP tolerances
      RelativeConvergenceCriterion(1e-4),
      AbsoluteStimulus(-1e7),
      PopulatedResult(false),
      FixedResult(false),
      UseAbsoluteStimulus(false),
      Converged(false),
      StimulateRegion(false)
    {
    }
    
    virtual void Converge()=0;   
    
    virtual ~AbstractUntemplatedConvergenceTester()
    {
    }
    
};

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM>
class AbstractConvergenceTester : public AbstractUntemplatedConvergenceTester
{

public:    
    void Converge()
    {
        
        // Create the meshes on which the test will be based
        const std::string mesh_dir = "ConvergenceMesh";
        OutputFileHandler output_file_handler(mesh_dir);
        ReplicatableVector voltage_replicated;

        unsigned file_num=0;
        
        SetInitialConvergenceParameters();
        
        unsigned prev_mesh_num=9999;
        std::string mesh_pathname;
        std::string mesh_filename;
        
        double prev_voltage[201];
        PopulateStandardResult(prev_voltage);
        do
        {
            CuboidMeshConstructor<DIM> constructor;

            if (this->MeshNum!=prev_mesh_num)
            {
                mesh_pathname = constructor.Construct(this->MeshNum, mMeshWidth);
                prev_mesh_num = this->MeshNum;
            }                            
            unsigned num_ele_across = (unsigned) pow(2, this->MeshNum+2); // number of elements in each dimension
            
            AbstractCardiacCellFactory<DIM>* p_cell_factory;
            if (!this->StimulateRegion)
            {
                //\todo The UseAbsoluteStimulus is temporary, while we are sorting out 
                //3D stimulus.  It is to be removed later (along with StimulusConvergenceTester)
          
                if (this->UseAbsoluteStimulus)
                {
                    #define COVERAGE_IGNORE
                    p_cell_factory = new GeneralPlaneStimulusCellFactory<CELL, DIM>(this->OdeTimeStep, 0, this->AbsoluteStimulus, true);
                    #undef COVERAGE_IGNORE                
                }
                else
                {
                    p_cell_factory = new GeneralPlaneStimulusCellFactory<CELL, DIM>(this->OdeTimeStep, num_ele_across, constructor.GetWidth());                
                }
            }
            else
            {
                p_cell_factory = new QuarterStimulusCellFactory<CELL, DIM>(this->OdeTimeStep, constructor.GetWidth());
            }
            
            CARDIAC_PROBLEM cardiac_problem(p_cell_factory);
            
            cardiac_problem.SetMeshFilename(mesh_pathname);
            cardiac_problem.SetOutputDirectory ("Convergence");
            cardiac_problem.SetOutputFilenamePrefix ("Results");
            
            cardiac_problem.SetEndTime(simulation_time);   // ms
            cardiac_problem.SetLinearSolverRelativeTolerance(this->KspRtol);
    
            cardiac_problem.SetPdeTimeStep(this->PdeTimeStep);
            
            assert(fabs(0.04/this->PdeTimeStep - round(0.04/this->PdeTimeStep)) <1e-15 );
            cardiac_problem.SetPrintingTimeStep(0.04);  //Otherwise we can't take the timestep down to machine precision without generating thousands of output files
            cardiac_problem.Initialise();
            
      	    DisplayRun();
            
            //// use this to get some info printed out
            //cardiac_problem.SetWriteInfo();
            
            try
            {
                cardiac_problem.Solve();
            }
            catch (Exception e)
            {
                #define COVERAGE_IGNORE
                //\todo Cover this
                std::cout<<"Warning - this run threw an exception.  Check convergence results\n";
                std::cout<<e.GetMessage() << std::endl;                 
                #undef COVERAGE_IGNORE
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
                    assert(this->PdeTimeStep<5);
                    first_quadrant_node = first_quadrant_nodes_3d[this->MeshNum];
                    third_quadrant_node = third_quadrant_nodes_3d[this->MeshNum];
                    break;
                }
                
                default:
                    assert(0);
            }
            
            #ifndef NDEBUG
            Node<DIM>* fqn = cardiac_problem.rGetMesh().GetNode(first_quadrant_node);
            Node<DIM>* tqn = cardiac_problem.rGetMesh().GetNode(third_quadrant_node);
            double mesh_width=constructor.GetWidth();
            assert(fqn->rGetLocation()[0]==0.25*mesh_width);
            assert(fabs(tqn->rGetLocation()[0] - 0.75*mesh_width) < 1e-10);
            for (unsigned coord=1; coord<DIM; coord++)
            {
                assert(fqn->rGetLocation()[coord]==0.5*mesh_width);
                assert(tqn->rGetLocation()[coord]==0.5*mesh_width);
            }
            #endif
            
            OutputFileHandler results_handler("Convergence", false);
            ColumnDataReader results_reader(results_handler.GetOutputDirectoryFullPath(), "Results", false);
            
            
            {
                std::vector<double> transmembrane_potential=results_reader.GetValues("V", third_quadrant_node);
                std::vector<double> time_series = results_reader.GetUnlimitedDimensionValues();
                
                // Write out the time series for the node at third quadrant
                if (results_handler.IsMaster())
                {
                    OutputFileHandler plot_file_handler("ConvergencePlots", false);
                    std::stringstream plot_file_name_stream;
                    plot_file_name_stream<< "Node1_"<< file_num << "_timestep.csv";
                    out_stream p_plot_file = plot_file_handler.OpenOutputFile(plot_file_name_stream.str());
                    for (unsigned data_point = 0; data_point<time_series.size(); data_point++)
                    {
                        (*p_plot_file) << time_series[data_point] << "\t" << transmembrane_potential[data_point] << "\n";                 
                    }
                    p_plot_file->close();
                }

                // calculate l2norm
                //double *p_prev_voltage = prev_voltage;
                double max_abs_error = 0;
                double sum_sq_abs_error =0;
                double sum_sq_prev_voltage = 0;
                
                
                for (unsigned data_point = 0; data_point<time_series.size(); data_point++)
                {
                    if (this->PopulatedResult)
                    {
                        double abs_error = fabs(transmembrane_potential[data_point]-prev_voltage[data_point]);
                        max_abs_error = (abs_error > max_abs_error) ? abs_error : max_abs_error;
                        sum_sq_abs_error += abs_error*abs_error;
                        sum_sq_prev_voltage += prev_voltage[data_point] * prev_voltage[data_point];
                    } 
                    
                    if (!this->PopulatedResult || !FixedResult)
                    {
                        prev_voltage[data_point] = transmembrane_potential[data_point];
                    }
                }
                if (this->PopulatedResult)
                {
                    std::cout << "max_abs_error = " << max_abs_error << " log10 = " << log10(max_abs_error) << "\n";
                    std::cout << "l2 error = " << sum_sq_abs_error/sum_sq_prev_voltage << " log10 = " << log10(sum_sq_abs_error/sum_sq_prev_voltage) << "\n";
                    //std::cout << log10(Abscissa()) << "\t" << log10(sum_sq_abs_error/sum_sq_prev_voltage) <<"\t#Logs for Gnuplot\n";
                    //Use "set logscale x; set logscale y" to get loglog plots in Gnuplot
                    std::cout << Abscissa() << "\t" << sum_sq_abs_error/sum_sq_prev_voltage <<"\t#Gnuplot raw data\n";
                    // convergence criterion
                    this->Converged = sum_sq_abs_error/sum_sq_prev_voltage<this->RelativeConvergenceCriterion;
                }
                if (!this->PopulatedResult)
                {
                    this->PopulatedResult=true;
                }
            }
            
            // Write time series for first quadrant node
            if (results_handler.IsMaster())
            {
                std::vector<double> transmembrane_potential=results_reader.GetValues("V", first_quadrant_node);
                std::vector<double> time_series = results_reader.GetUnlimitedDimensionValues();
                OutputFileHandler plot_file_handler("ConvergencePlots", false);
                std::stringstream plot_file_name_stream;
                plot_file_name_stream<< "Node2_"<< file_num << "_timestep.csv";
                out_stream p_plot_file = plot_file_handler.OpenOutputFile(plot_file_name_stream.str());
                for (unsigned data_point = 0; data_point<time_series.size(); data_point++)
                {
                    (*p_plot_file) << time_series[data_point] << "\t" << transmembrane_potential[data_point] << "\n";                 
                }
                p_plot_file->close();
            }                    
            
            // Get ready for the next test by halving the time step
            if (!this->Converged)
            {
                UpdateConvergenceParameters();
                file_num++;
            }
            delete p_cell_factory;
        }
        while (!GiveUpConvergence() && !this->Converged);
    }
    
    void DisplayRun()
    {
        unsigned num_ele_across = (unsigned) pow(2, this->MeshNum+2);// number of elements in each dimension
        double scaling = mMeshWidth/(double) num_ele_across;
        
        std::cout<<"================================================================================"<<std::endl;
        std::cout<<"Solving with a space step of "<< scaling << " cm (mesh " << this->MeshNum << ")" << std::endl;
        std::cout<<"Solving with a time step of "<<this->PdeTimeStep<<" ms"<<std::endl;
        std::cout<<"Solving with an ode time step of "<<this->OdeTimeStep<<" ms"<<std::endl;
        std::cout<<"Solving with a KSP relative tolerance of "<<this->KspRtol<<std::endl;
        std::cout<<"Solving with stimulating a quarter of the mesh? " << this->StimulateRegion<<std::endl;
        system("date");//To keep track of what Nightly things are doing
        //\todo The UseAbsoluteStimulus is temporary, while we are sorting out 
        //3D stimulus.  It is to be removed later (along with StimulusConvergenceTester)
        if (this->UseAbsoluteStimulus)
        {
            #define COVERAGE_IGNORE
            std::cout<<"Using absolute stimulus of "<<this->AbsoluteStimulus<<std::endl;
            #undef COVERAGE_IGNORE                
        }
        std::cout << std::flush;
        
    }
    
public:
    virtual ~AbstractConvergenceTester() {}
    
    virtual void SetInitialConvergenceParameters()=0;
    virtual void UpdateConvergenceParameters()=0;
    virtual bool GiveUpConvergence()=0;
    virtual double Abscissa()=0;
    virtual void PopulateStandardResult(double result[])
    {
        assert(this->PopulatedResult==false);
    }
    
    bool IsConverged()
    {
        return Converged;
    }
    
    void SetMeshWidth(double meshWidth)
    {
    	mMeshWidth=meshWidth;
    }
};
#endif /*ABSTRACTCONVERGENCETESTER_HPP_*/
