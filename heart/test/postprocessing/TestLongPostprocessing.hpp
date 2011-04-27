/*

Copyright (C) University of Oxford, 2005-2010

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

#ifndef TESTLONGPOSTPROCESSING_HPP_
#define TESTLONGPOSTPROCESSING_HPP_


#include <cxxtest/TestSuite.h>
#include <ctime>

#include "CheckpointArchiveTypes.hpp" // Needs to be before other Chaste code
#include "CardiacSimulationArchiver.hpp"
#include "ArchiveOpener.hpp"
#include "Exception.hpp"
#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/shared_ptr.hpp>

#include "AbstractCardiacCell.hpp"
#include "AbstractCvodeCell.hpp"
#include "AbstractCardiacCellWithModifiers.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "Exception.hpp"

#include "CellProperties.hpp"

#include "SimpleStimulus.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudy1991.hpp"

#include "MonodomainProblem.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "PostProcessingWriter.hpp"

#include "CommandLineArguments.hpp"

#include <sstream>

template <unsigned DIM>
class PointStimulusCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    boost::shared_ptr<RegularStimulus> mpStimulus;
    double mAreaOrVolume;

public:
    PointStimulusCellFactory(const double& rStimMag, const double& rStimDuration, const double& rPacingCycleLength, const double& rAreaOrVolume )
        : AbstractCardiacCellFactory<DIM>(),
          mpStimulus(new RegularStimulus(rStimMag, rStimDuration, rPacingCycleLength, 0)),// Default stimulus introduced at t=0.
          mAreaOrVolume(rAreaOrVolume)
    {
        assert(DIM==1 || DIM==2);
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x=0, y=0, z=0, z_target=0;
        x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];
        y = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[1];

        ChasteCuboid<DIM> extremes = this->GetMesh()->CalculateBoundingBox();
        // double x_middle = (extremes.rGetLowerCorner()[0]+extremes.rGetUpperCorner()[0])/2.0;
        // double y_middle = (extremes.rGetLowerCorner()[1]+extremes.rGetUpperCorner()[1])/2.0;
        //double endo_proportion = x/extremes.GetWidth(0);

        // Get a threshold for each dimension by getting sqrt or cube root of the area/volume.
        // double threshold = pow(mAreaOrVolume,1.0/(double)(DIM));

        if (DIM==3)
        {
            z = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[2];
            z_target = extremes.rGetUpperCorner()[2];
        }
        // std::cout << "x location "<< x << " y location "<< y;
        // We want a stimulus in the middle of x and y and on the top surface of the the mesh.
        if (x < 0.03)
        {
        	return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpStimulus);
        	// std::cout << "Cell with stimulus created\n" << std::flush;
        }
        else
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpZeroStimulus);
        }
    }
};

class TestLongPostprocessing : public CxxTest::TestSuite
{

public:


    void Test2DSimulations() throw(Exception)
    {
//        // Get command line arguments
//        CommandLineArguments* p_args = CommandLineArguments::Instance();
//        unsigned argc = *(p_args->p_argc); // has the number of arguments.
//        std::cout << "# " << argc-1 << " arguments supplied.\n" << std::flush;

        // Define conductivity scale

        double conductivity_scale = 1;
//
//        bool Conductivity_Scale_Option = CommandLineArguments::Instance()->OptionExists("--conductivity_scale");
//        if (Conductivity_Scale_Option == true) {
//            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--conductivity_scale");
//            conductivity_scale = atof(val);
//        }
//
//        std::cout << "conductivity_scale is "<< conductivity_scale << "\n";

        double h = 0.01; // cm

//        bool Spatial_Discretization_Option = CommandLineArguments::Instance()->OptionExists("--spatial_discretization");
//        if (Spatial_Discretization_Option == true) {
//            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--spatial_discretization");
//            h = atof(val);
//        }
//
//        std::cout << "spatial_discretization is "<< h << "\n";

//        // This option is not used at present.
//        double simulation_duration = 5.0; //ms
//
//        bool Simulation_Duration_Option = CommandLineArguments::Instance()->OptionExists("--simulation_duration");
//        if (Simulation_Duration_Option == true) {
//            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--simulation_duration");
//            simulation_duration = atof(val);
//        }
//
//        std::cout << "simulation_duration is "<< simulation_duration << "\n";

        double ode_time_step = 0.005; //ms

//        bool ODE_Time_Step_Option = CommandLineArguments::Instance()->OptionExists("--ode_time_step");
//        if (ODE_Time_Step_Option == true) {
//            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--ode_time_step");
//            ode_time_step = atof(val);
//        }
//
//        std::cout << "ode_time_step is "<< ode_time_step << "\n";

        double pde_time_step = 0.01; //ms

//        bool PDE_Time_Step_Option = CommandLineArguments::Instance()->OptionExists("--pde_time_step");
//        if (PDE_Time_Step_Option == true) {
//            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--pde_time_step");
//            pde_time_step = atof(val);
//        }
//
//        std::cout << "pde_time_step is "<< pde_time_step << "\n";

        unsigned num_stims = 1;

//        bool Num_Stim_Option = CommandLineArguments::Instance()->OptionExists("--num_stims");
//        if (Num_Stim_Option == true) {
//            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--num_stims");
//            num_stims = atof(val);
//        }
        // Check at least one stimulus is applied
        assert(num_stims >= 1);

        std::cout << "num_stims is "<< num_stims << "\n";

        TetrahedralMesh<2,2> mesh;
        unsigned num_elem_x = (unsigned)(1.0/h);  // num elements to make 1cm
        unsigned num_elem_y = (unsigned)(1.0/h);  // num elements to make 1cm
        //unsigned num_elem_z = (unsigned)(0.15/h);// Num elements to make 0.3cm
        double pacing_cycle_length = 250;
        double stim_mag = -500;
        double stim_dur = 3;
        double area = 0.005;

        mesh.ConstructRectangularMesh(num_elem_x, num_elem_y);
        mesh.Scale(h,h); // Get mesh into units of cm.

        std::string archive_dir_base("LongPostprocessing_archives/archive");
        std::string archive_dir_current;

        // Setup
        HeartConfig::Instance()->SetSimulationDuration(pacing_cycle_length); //ms
        HeartConfig::Instance()->SetOutputDirectory("LongPostprocessing");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        
        // These lines make postprocessing fast or slow.
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(ode_time_step,pde_time_step,1);
        //HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(ode_time_step,pde_time_step,0.01);
        
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.4*conductivity_scale*1.171, 1.4*conductivity_scale*1.171));
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400.0); // 1/cm
        HeartConfig::Instance()->SetCapacitance(1.0); // uF/cm^2

        std::vector<std::pair<double,double> > apds_requested;
        apds_requested.push_back(std::pair<double, double>(90,-30)); //repolarisation percentage and threshold
        HeartConfig::Instance()->SetApdMaps(apds_requested);
//        std::vector<double> excitation_threshold;
//        excitation_threshold.push_back(double (-30));
//        HeartConfig::Instance()->SetUpstrokeTimeMaps(excitation_threshold);
//        HeartConfig::Instance()->SetMaxUpstrokeVelocityMaps(excitation_threshold);

        for (unsigned stim_counter=0; stim_counter < num_stims; stim_counter++ )
        {
            // Load problem
            MonodomainProblem<2> *p_monodomain_problem;
            if (stim_counter==0)
            {
                PointStimulusCellFactory<2> cell_factory(stim_mag, stim_dur, pacing_cycle_length, area);
                p_monodomain_problem = new MonodomainProblem<2>( &cell_factory );
                p_monodomain_problem->SetMesh(&mesh);
                p_monodomain_problem->Initialise();
            }
            else
            {
                p_monodomain_problem = CardiacSimulationArchiver<MonodomainProblem<2> >::Load(archive_dir_current);
            }

            HeartConfig::Instance()->SetSimulationDuration((double) (stim_counter+1)*pacing_cycle_length); //ms

            // set new directories to work from
            std::stringstream stringoutput;
            stringoutput << stim_counter;
            std::string stim_counter_string = stringoutput.str();

            archive_dir_current = archive_dir_base + "_" + stim_counter_string;
            HeartConfig::Instance()->SetOutputFilenamePrefix("results_" + stim_counter_string);

            // Solve problem (this does the postprocessing too when HeartConfig options are set).
            p_monodomain_problem->Solve();

            // Save problem to archive
            CardiacSimulationArchiver<MonodomainProblem<2> >::Save(*p_monodomain_problem, archive_dir_current, false);
            std::cout << "Archived to " << archive_dir_current << "\n" << std::flush;

            if (PetscTools::AmMaster()) // Best to have only the master processor messing with files
            {
                // Copy the postprocessing results into the archive folders so they aren't wiped.
                std::vector<std::string> files;
                files.push_back("Apd_90_-30_Map");
//                files.push_back("MaxUpstrokeVelocityMap_-30");
//                files.push_back("UpstrokeTimeMap_-30");
                for (unsigned i=0; i<files.size(); i++)
                {
                    std::string command = "mv " +  OutputFileHandler::GetChasteTestOutputDirectory() + HeartConfig::Instance()->GetOutputDirectory() + "/output/" + files[i] + ".dat "
                                          + OutputFileHandler::GetChasteTestOutputDirectory() + archive_dir_current + "/" + files[i] + ".dat";

                    // If we are in an PetscTools::AmMaster() block we need to use this instead of EXPECT0.
                    MPIABORTIFNON0(system, command);
                }
            }
            PetscTools::Barrier();
        }// close for loop
    }//close void Test2dSimulations

};

#endif /* TESTLONGPOSTPROCESSING_HPP_ */

