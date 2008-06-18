/*

Copyright (C) University of Oxford, 2008

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


#ifndef TEST1DFIBREWITHHETEROGENEITY_HPP_
#define TEST1DFIBREWITHHETEROGENEITY_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>

#include "ConformingTetrahedralMesh.hpp"
#include "Hdf5DataReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"

#include "FaberRudy2000Version3.cpp"
#include "AbstractCardiacCellFactory.hpp"

class HeterogeneousCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    // define a new stimulus
    SimpleStimulus* mpStimulus;

public:
    HeterogeneousCellFactory() : AbstractCardiacCellFactory<1>(0.005)//Ode timestep
    {
        // set the new stimulus
        mpStimulus = new SimpleStimulus(-600, 0.5);
    }

    HeterogeneousCellFactory(double timeStep, double stimulusMagnitude) : AbstractCardiacCellFactory<1>(timeStep)
    {
        // set the new stimulus
        mpStimulus = new SimpleStimulus(stimulusMagnitude, 0.5);
    }

    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        FaberRudy2000Version3 *cell;

        if (this->mpMesh->GetNode(node)->GetPoint()[0] == 0.0)
        {
            cell = new FaberRudy2000Version3(this->mpSolver,
                                             this->mTimeStep,
                                             mpStimulus,
                                             this->mpZeroStimulus);

        }
        else
        {
            cell = new FaberRudy2000Version3(this->mpSolver,
                                             this->mTimeStep,
                                             this->mpZeroStimulus,
                                             this->mpZeroStimulus);
        }

        if (this->mpMesh->GetNode(node)->GetPoint()[0] < 0.3333)
        {
            cell->SetScaleFactorGks(0.462);
            cell->SetScaleFactorIto(0.0);
        }
        else if (this->mpMesh->GetNode(node)->GetPoint()[0] < 0.6666)
        {
            cell->SetScaleFactorGks(1.154);
            cell->SetScaleFactorIto(0.85);
        }
        else //this->mpMesh->GetNode(node)->GetPoint()[0] < 1
        {
            cell->SetScaleFactorGks(1.154);
            cell->SetScaleFactorIto(1.0);
        }

        return cell;
    }

    ~HeterogeneousCellFactory(void)
    {
        delete mpStimulus;
    }
};


class Test1dFibreWithHeterogeneity : public CxxTest::TestSuite
{
public:
    // Solve on a 1D string of cells, 1cm long with a space step of 0.1mm and heterogeneous cell types.
    void TestFibreHeterogeneity()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));

        HeterogeneousCellFactory cell_factory;
        MonodomainProblem<1> monodomain_problem(&cell_factory);

        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        monodomain_problem.SetEndTime(300);   // ms
        monodomain_problem.SetPdeAndPrintingTimeSteps(0.01, 0.1);
        monodomain_problem.SetOutputDirectory("FibreWithHeterogeneity");
        monodomain_problem.SetOutputFilenamePrefix("Monodomain1d");

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        monodomain_problem.GetMonodomainPde()->SetCapacitance(1.0);

        monodomain_problem.Solve();


        // write out results for node 20 (and 50 and 80)
        OutputFileHandler results_handler("FibreWithHeterogeneity", false);
        Hdf5DataReader results_reader(results_handler.GetOutputDirectoryFullPath(), "Monodomain1d", false);

        unsigned relevant_nodes[3]={20,50,80};

        for (unsigned i=0; i<3; i++)
        {
            std::vector<double> transmembrane_potential=results_reader.GetVariableOverTime("V", relevant_nodes[i]);
            std::vector<double> time_series = results_reader.GetUnlimitedDimensionValues();

            // Write out the time series for the node at third quadrant
            if (results_handler.IsMaster())
            {
                OutputFileHandler plot_file_handler("HeterogeneityPlots", false);
                std::stringstream plot_file_name_stream;
                plot_file_name_stream<< "Node_" << relevant_nodes[i] << ".csv";
                out_stream p_plot_file = plot_file_handler.OpenOutputFile(plot_file_name_stream.str());
                for (unsigned data_point = 0; data_point<time_series.size(); data_point++)
                {
                    (*p_plot_file) << time_series[data_point] << "\t" << transmembrane_potential[data_point] << "\n";
                }
                p_plot_file->close();

                std::stringstream cmd;
                cmd << "ndiff " << plot_file_handler.GetChasteTestOutputDirectory() << "HeterogeneityPlots/" << plot_file_name_stream.str() << " heart/test/data/HeterogeneityPlots/Node_" << relevant_nodes[i] << ".csv";
                TS_ASSERT_EQUALS(system(cmd.str().c_str()), 0);
            }
        }

    }
};

#endif /*TEST1DFIBREWITHHETEROGENEITY_HPP_*/
