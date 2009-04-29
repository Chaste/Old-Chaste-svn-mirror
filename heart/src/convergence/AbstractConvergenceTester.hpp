/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef ABSTRACTCONVERGENCETESTER_HPP_
#define ABSTRACTCONVERGENCETESTER_HPP_

#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include "AbstractMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshWriter.hpp"
#include "PropagationPropertiesCalculator.hpp"
#include "Hdf5DataReader.hpp"
#include "GeneralPlaneStimulusCellFactory.hpp"
#include "CuboidMeshConstructor.hpp"
#include "OutputFileHandler.hpp"
#include "ZeroStimulusCellFactory.hpp"
#include "SimpleStimulus.hpp"
#include "ConstBoundaryCondition.hpp"
#include "StimulusBoundaryCondition.hpp"


typedef enum StimulusType_
{
    PLANE=0,
    REGION,
    NEUMANN
} StimulusType;

/**
 * \todo Documentation...
 */
template <class CELL, unsigned DIM>
class QuarterStimulusCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    // define a new stimulus
    SimpleStimulus* mpStimulus;
    double mMeshWidth;
public:
    QuarterStimulusCellFactory(double meshWidth) : AbstractCardiacCellFactory<DIM>()
    {
        mpStimulus = new SimpleStimulus(-1000000, 0.5);
        mMeshWidth=meshWidth;
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        double x = this->GetMesh()->GetNode(node)->GetPoint()[0];
        if (x<=mMeshWidth*0.25+1e-10)
        {
            return new CELL(this->mpSolver, this->mpStimulus);
        }
        else
        {
            return new CELL(this->mpSolver, this->mpZeroStimulus);
        }
    }

    ~QuarterStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};


/**
 * \todo Documentation...
 */
class AbstractUntemplatedConvergenceTester
{
protected:
    double mMeshWidth;
    double mKspTolerance;
    bool mUseKspAbsoluteTolerance;
public:
    double OdeTimeStep;
    double PdeTimeStep;
    unsigned MeshNum;
    double RelativeConvergenceCriterion;
    double LastDifference;
    double AbsoluteStimulus;
    bool PopulatedResult;
    bool FixedResult;
    bool UseAbsoluteStimulus;
    //bool UseNeumannStimulus;
    bool Converged;
    //bool StimulateRegion;
    StimulusType Stimulus;
    double NeumannStimulus;

    AbstractUntemplatedConvergenceTester();

    virtual void Converge(std::string nameOfTest)=0;

    void SetKspRelativeTolerance(const double relativeTolerance);

    void SetKspAbsoluteTolerance(const double absoluteTolerance);

    double GetKspAbsoluteTolerance();

    double GetKspRelativeTolerance();

    virtual ~AbstractUntemplatedConvergenceTester();
};

/**
 * Use template specialization to set the appropriate conductivities for the problem.
 */
template<class CARDIAC_PROBLEM, unsigned DIM>
void SetConductivities(CARDIAC_PROBLEM& rCardiacProblem);

template<unsigned DIM>
void SetConductivities(BidomainProblem<DIM>& rProblem)
{
    c_vector<double, DIM> conductivities;
    for (unsigned i=0; i<DIM; i++)
    {
        conductivities[i] = 1.75;
    }
    HeartConfig::Instance()->SetIntracellularConductivities(conductivities);

    for (unsigned i=0; i<DIM; i++)
    {
        conductivities[i] = 7.0;
    }
    HeartConfig::Instance()->SetExtracellularConductivities(conductivities);
}

template<unsigned DIM>
void SetConductivities(MonodomainProblem<DIM>& rProblem)
{
    c_vector<double, DIM> conductivities;
    for (unsigned i=0; i<DIM; i++)
    {
        conductivities[i] = 1.75;
    }
    HeartConfig::Instance()->SetIntracellularConductivities(conductivities);
}



/**
 * \todo Documentation...
 */
template<class CELL, class CARDIAC_PROBLEM, unsigned DIM, unsigned PROBLEM_DIM>
class AbstractConvergenceTester : public AbstractUntemplatedConvergenceTester
{
public:
    /**
     * \todo This is a scarily long method; could do with some parts extracted?
     */
    void Converge(std::string nameOfTest)
    {
        std::cout << "=========================== Beginning Test...==================================\n";
        // Create the meshes on which the test will be based
        const std::string mesh_dir = "ConvergenceMesh";
        OutputFileHandler output_file_handler(mesh_dir);
        ReplicatableVector voltage_replicated;

        unsigned file_num=0;

        // Create a file for storing conduction velocity and AP data and write the header
        OutputFileHandler conv_info_handler("ConvergencePlots", false);
        out_stream p_conv_info_file;
        if (conv_info_handler.IsMaster())
        {
            p_conv_info_file = conv_info_handler.OpenOutputFile(nameOfTest+"_info.csv");
            (*p_conv_info_file) << "#Abcisa\t"
                                << "(l2-norm)^2\t"
                                << "l2-norm\t"
                                << "Max absolute err\t"
                                << "APD90_1st_quad\t"
                                << "APD90_3rd_quad\t"
                                << "Conduction velocity (relative errors)" << std::endl;
        }
        SetInitialConvergenceParameters();

        double prev_apd90_first_qn=0.0;
        double prev_apd90_third_qn=0.0;
        double prev_cond_velocity=0.0;
        double prev_voltage[201];
        PopulateStandardResult(prev_voltage);

        do
        {
            CuboidMeshConstructor<DIM> constructor;

            assert(fabs(0.04/this->PdeTimeStep - round(0.04/this->PdeTimeStep)) <1e-15 );
            //Otherwise we can't take the timestep down to machine precision without generating thousands of output files
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(this->OdeTimeStep, this->PdeTimeStep, 0.04);
            HeartConfig::Instance()->SetSimulationDuration(8.0);
            HeartConfig::Instance()->SetOutputDirectory ("Convergence");
            HeartConfig::Instance()->SetOutputFilenamePrefix ("Results");

            //Don't use parallel mesh for now
            //HeartConfig::Instance()->SetMeshFileName( constructor.Construct(this->MeshNum, mMeshWidth) );

            TetrahedralMesh<DIM, DIM> mesh;
            TrianglesMeshReader<DIM, DIM> mesh_reader(constructor.Construct(this->MeshNum, mMeshWidth) );
            mesh.ConstructFromMeshReader(mesh_reader);

            unsigned num_ele_across = (unsigned) pow(2, this->MeshNum+2); // number of elements in each dimension

            AbstractCardiacCellFactory<DIM>* p_cell_factory=NULL;

            switch (this->Stimulus)
            {
                case NEUMANN:
                {
                    p_cell_factory = new ZeroStimulusCellFactory<CELL, DIM>();
                    break;
                }
                case PLANE:
                {
                    if (this->UseAbsoluteStimulus)
                    {
                        #define COVERAGE_IGNORE
                        p_cell_factory = new GeneralPlaneStimulusCellFactory<CELL, DIM>(0, this->AbsoluteStimulus, true);
                        #undef COVERAGE_IGNORE
                    }
                    else
                    {
                        p_cell_factory = new GeneralPlaneStimulusCellFactory<CELL, DIM>(num_ele_across, constructor.GetWidth());
                    }
                    break;
                }
                case REGION:
                {
                    p_cell_factory = new QuarterStimulusCellFactory<CELL, DIM>(constructor.GetWidth());
                    break;
                }
            }


            CARDIAC_PROBLEM cardiac_problem(p_cell_factory);
            ///\todo this is a sequential mesh
            cardiac_problem.SetMesh(&mesh);
            if (mUseKspAbsoluteTolerance)
            {
                HeartConfig::Instance()->SetUseAbsoluteTolerance(this->mKspTolerance);

//                cardiac_problem.SetLinearSolverAbsoluteTolerance(this->mKspTolerance);
            }
            else
            {
                HeartConfig::Instance()->SetUseRelativeTolerance(this->mKspTolerance);

//                cardiac_problem.SetLinearSolverRelativeTolerance(this->mKspTolerance);
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

            double mesh_width=constructor.GetWidth();

            // We only need the output of these two nodes
            std::vector<unsigned> nodes_to_be_output;
            nodes_to_be_output.push_back(first_quadrant_node);
            nodes_to_be_output.push_back(third_quadrant_node);
            cardiac_problem.SetOutputNodes(nodes_to_be_output);
            ///\todo #606 Added back in about r4058

            // The results of the tests were originally obtained with the following conductivity
            // values. After implementing fibre orientation the defaults changed. Here we set
            // the former ones to be used.
            SetConductivities(cardiac_problem);

            cardiac_problem.Initialise();

            #ifndef NDEBUG
            Node<DIM>* fqn = cardiac_problem.rGetMesh().GetNode(first_quadrant_node);
            Node<DIM>* tqn = cardiac_problem.rGetMesh().GetNode(third_quadrant_node);
            assert(fqn->rGetLocation()[0]==0.25*mesh_width);
            assert(fabs(tqn->rGetLocation()[0] - 0.75*mesh_width) < 1e-10);
            for (unsigned coord=1; coord<DIM; coord++)
            {
                assert(fqn->rGetLocation()[coord]==0.5*mesh_width);
                assert(tqn->rGetLocation()[coord]==0.5*mesh_width);
            }
            #endif

            BoundaryConditionsContainer<DIM,DIM,PROBLEM_DIM> bcc;
            SimpleStimulus stim(NeumannStimulus, 0.5);
            if (Stimulus==NEUMANN)
            {

                StimulusBoundaryCondition<DIM> *p_bc_stim = new StimulusBoundaryCondition<DIM>(&stim);

                // get mesh
                AbstractMesh<DIM, DIM> &r_mesh = cardiac_problem.rGetMesh();
                // loop over boundary elements
                typename AbstractMesh<DIM, DIM>::BoundaryElementIterator iter;
                iter = r_mesh.GetBoundaryElementIteratorBegin();
                while (iter != r_mesh.GetBoundaryElementIteratorEnd())
                {
                    double x = ((*iter)->CalculateCentroid())[0];
                    if (x*x<=1e-10)
                    {
                        bcc.AddNeumannBoundaryCondition(*iter, p_bc_stim);
                    }
                    iter++;
                }
                // pass the bcc to the problem
                cardiac_problem.SetBoundaryConditionsContainer(&bcc);
            }

              DisplayRun();
            double time_before=MPI_Wtime();
            //// use this to get some info printed out
            //cardiac_problem.SetWriteInfo();

            try
            {
                cardiac_problem.Solve();
            }
            catch (Exception e)
            {
                #define COVERAGE_IGNORE
                ///\todo Cover this
                std::cout<<"Warning - this run threw an exception.  Check convergence results\n";
                std::cout<<e.GetMessage() << std::endl;
                #undef COVERAGE_IGNORE
            }

            std::cout << "Time to solve = "<<MPI_Wtime()-time_before<<" seconds\n";

            OutputFileHandler results_handler("Convergence", false);
            Hdf5DataReader results_reader = cardiac_problem.GetDataReader();

            {
                std::vector<double> transmembrane_potential=results_reader.GetVariableOverTime("V", third_quadrant_node);
                std::vector<double> time_series = results_reader.GetUnlimitedDimensionValues();

                // Write out the time series for the node at third quadrant
                if (results_handler.IsMaster())
                {
                    OutputFileHandler plot_file_handler("ConvergencePlots", false);
                    std::stringstream plot_file_name_stream;
                    plot_file_name_stream<< nameOfTest << "_Third_quadrant_node_run_"<< file_num << ".csv";
                    out_stream p_plot_file = plot_file_handler.OpenOutputFile(plot_file_name_stream.str());
                    for (unsigned data_point = 0; data_point<time_series.size(); data_point++)
                    {
                        (*p_plot_file) << time_series[data_point] << "\t" << transmembrane_potential[data_point] << "\n";
                    }
                    p_plot_file->close();
                }

                // Write time series for first quadrant node
                if (results_handler.IsMaster())
                {
                    std::vector<double> transmembrane_potential_1qd=results_reader.GetVariableOverTime("V", first_quadrant_node);
                    std::vector<double> time_series_1qd = results_reader.GetUnlimitedDimensionValues();
                    OutputFileHandler plot_file_handler("ConvergencePlots", false);
                    std::stringstream plot_file_name_stream;
                    plot_file_name_stream<< nameOfTest << "_First_quadrant_node_run_"<< file_num << ".csv";
                    out_stream p_plot_file = plot_file_handler.OpenOutputFile(plot_file_name_stream.str());
                    for (unsigned data_point = 0; data_point<time_series.size(); data_point++)
                    {
                        (*p_plot_file) << time_series_1qd[data_point] << "\t" << transmembrane_potential_1qd[data_point] << "\n";
                    }
                    p_plot_file->close();
                }

                // calculate conduction velocity and APD90 error
                PropagationPropertiesCalculator ppc(&results_reader);

                double cond_velocity=0.0, apd90_first_qn=0.0, apd90_third_qn=0.0;
                try
                {
                    apd90_first_qn = ppc.CalculateActionPotentialDuration(0.9, first_quadrant_node);
                    apd90_third_qn = ppc.CalculateActionPotentialDuration(0.9, third_quadrant_node);

                    cond_velocity  = ppc.CalculateConductionVelocity(first_quadrant_node,third_quadrant_node,0.5*mesh_width);
                }
                catch (Exception e)
                {
                    #define COVERAGE_IGNORE
                    std::cout<<"Warning - this run threw an exception in calculating propagation.  Check convergence results\n";
                    std::cout<<e.GetMessage() << std::endl;
                    #undef COVERAGE_IGNORE
                }

                double cond_velocity_error = 0.0;
                double apd90_first_qn_error = 0.0;
                double apd90_third_qn_error = 0.0;

                if (this->PopulatedResult)
                {
                    //std::cout << "APD90\n"
                    //          << "Current: " << apd90_first_qn << "\t" << apd90_third_qn << "\n"
                    //          << "Previous: " << prev_apd90_first_qn << "\t" << prev_apd90_third_qn << "\n";

                    cond_velocity_error = fabs(cond_velocity - prev_cond_velocity) / prev_cond_velocity;
                    apd90_first_qn_error = fabs(apd90_first_qn - prev_apd90_first_qn) / prev_apd90_first_qn;
                    apd90_third_qn_error = fabs(apd90_third_qn - prev_apd90_third_qn) / prev_apd90_third_qn;
                }

                prev_cond_velocity = cond_velocity;
                prev_apd90_first_qn = apd90_first_qn;
                prev_apd90_third_qn = apd90_third_qn;

                // calculate l2norm
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
                    //std::cout << Abscissa() << "\t" << sum_sq_abs_error/sum_sq_prev_voltage <<"\t#Gnuplot raw data\n";

                    if (conv_info_handler.IsMaster())
                    {
                        (*p_conv_info_file) << std::setprecision(8)
                                            << Abscissa() << "\t"
                                            << sum_sq_abs_error/sum_sq_prev_voltage << "\t"
                                            << sqrt(sum_sq_abs_error/sum_sq_prev_voltage) << "\t"
                                            << max_abs_error << "\t"
                                            << apd90_first_qn_error << "\t"
                                            << apd90_third_qn_error << "\t"
                                            << cond_velocity_error  << std::endl;
                    }
                    // convergence criterion
                    this->Converged = sum_sq_abs_error/sum_sq_prev_voltage<this->RelativeConvergenceCriterion;
                    this->LastDifference=sum_sq_abs_error/sum_sq_prev_voltage;
                }

                if (!this->PopulatedResult)
                {
                    this->PopulatedResult=true;

                }
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

        if (conv_info_handler.IsMaster())
        {
            p_conv_info_file->close();

            std::cout << "Results: " << std::endl;
            EXPECT0(system, "cat " + conv_info_handler.GetOutputDirectoryFullPath() + nameOfTest + "_info.csv");
        }

    }

    void DisplayRun()
    {
        unsigned num_ele_across = (unsigned) pow(2, this->MeshNum+2);// number of elements in each dimension
        double scaling = mMeshWidth/(double) num_ele_across;

        std::cout<<"================================================================================"<<std::endl;
        std::cout<<"Solving in "<<DIM<<"D\n";
        std::cout<<"Solving with a space step of "<< scaling << " cm (mesh " << this->MeshNum << ")" << std::endl;
        std::cout<<"Solving with a time step of "<<this->PdeTimeStep<<" ms"<<std::endl;
        std::cout<<"Solving with an ode time step of "<<this->OdeTimeStep<<" ms"<<std::endl;
        if (mUseKspAbsoluteTolerance)
        {
            std::cout<<"Solving with a KSP absolute tolerance of "<<this->mKspTolerance<<std::endl;
        }
        else
        {
            std::cout<<"Solving with a KSP relative tolerance of "<<this->mKspTolerance<<std::endl;
        }
        switch (this->Stimulus)
        {
            case PLANE:
            std::cout<<"Stimulus = Plane\n";
            break;

            case REGION:
            std::cout<<"Stimulus = Region\n";
            break;

            case NEUMANN:
            std::cout<<"Stimulus = Neumann\n";
            break;

        }
        EXPECT0(system, "date");//To keep track of what Nightly things are doing
        ///\todo The UseAbsoluteStimulus is temporary, while we are sorting out
        ///3D stimulus.  It is to be removed later (along with StimulusConvergenceTester)
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
