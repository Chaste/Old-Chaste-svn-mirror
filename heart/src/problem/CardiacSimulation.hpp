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

#ifndef CARDIACSIMULATION_HPP_
#define CARDIACSIMULATION_HPP_

// Must go first
#include "CardiacSimulationArchiver.hpp"
// Note that since we include this header, we can't (easily) create a .cpp file for this class.
// Doing so would break the build with chaste_libs=0 on linking, since it would mean that no test
// (or executable) could include CardiacSimulation.hpp, since Boost archiving GUIDs would then be
// defined in two object files (CardiacSimulation.o and the test runner .o).  You get errors such
// as: multiple definition of `boost::archive::detail::guid_initializer<HeartConfig>::instance'.

#include <vector>
#include <ctime>
#include <memory>

#include "UblasIncludes.hpp"

#include "MonodomainProblem.hpp"
#include "BidomainProblem.hpp"
#include "PetscTools.hpp"
#include "TimeStepper.hpp"
#include "Exception.hpp"

#include "HeartConfig.hpp"
#include "HeartConfigRelatedCellFactory.hpp"

#include "TetrahedralMesh.hpp"
#include "NonCachedTetrahedralMesh.hpp"
#include "ChastePoint.hpp"
#include "ChasteCuboid.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include "TrianglesMeshWriter.hpp"

#include "OrthotropicConductivityTensors.hpp"

#include "Hdf5ToMeshalyzerConverter.hpp"
#include "PostProcessingWriter.hpp"

#include "OutputDirectoryFifoQueue.hpp"

/**
 * A class which encapsulates the executable functionality.
 *
 * Takes in a chaste parameters xml file and runs the relevant simulation.
 *
 * \todo High level user documentation.
 * This should describe the functionality available from the XML file.
 * It should include information on the structure of the output directory, especially when checkpointing,
 * and how to resume a simulation.
 */
class CardiacSimulation
{
private:
    /**
     * Read parameters from the HeartConfig XML file.
     *
     * @param parameterFileName a string containing the chaste simulation parameters XML file name.
     */
    void ReadParametersFromFile(std::string parameterFileName);

    /**
     * Templated method which creates and runs a cardiac simulation, based on the
     * XML file passed to our constructor.
     */
    template<class Problem, unsigned SPACE_DIM>
    void CreateAndRun()
    {
        Problem* p_problem;

        if (HeartConfig::Instance()->IsSimulationDefined())
        {
            HeartConfigRelatedCellFactory<SPACE_DIM> cell_factory;
            p_problem = new Problem(&cell_factory);

            p_problem->Initialise();
        }
        else // (HeartConfig::Instance()->IsSimulationResumed())
        {
            p_problem = CardiacSimulationArchiver<Problem>::Load(HeartConfig::Instance()->GetArchivedSimulationDir());
        }

        if (HeartConfig::Instance()->GetCheckpointSimulation())
        {
            // Create the checkpoints directory and set up a fifo queue of subdirectory names
            OutputDirectoryFifoQueue directory_queue(HeartConfig::Instance()->GetOutputDirectory() + "_checkpoints/",
                                                     HeartConfig::Instance()->GetMaxCheckpointsOnDisk());

            TimeStepper checkpoint_stepper(0.0, HeartConfig::Instance()->GetSimulationDuration(), HeartConfig::Instance()->GetCheckpointTimestep());
            while ( !checkpoint_stepper.IsTimeAtEnd() )
            {
                // Solve checkpoint timestep
                HeartConfig::Instance()->SetSimulationDuration(checkpoint_stepper.GetNextTime());
                p_problem->Solve();

                // Create directory that will contain archive and partial results for this checkpoint timestep.
                std::stringstream checkpoint_id;
                checkpoint_id << HeartConfig::Instance()->GetSimulationDuration() << "ms/";
                std::string checkpoint_dir_basename = directory_queue.CreateNextDir(checkpoint_id.str());

                // Archive simulation (in a subdirectory of checkpoint_dir_basename).
                std::stringstream archive_foldername;
                archive_foldername << HeartConfig::Instance()->GetOutputDirectory() << "_" << HeartConfig::Instance()->GetSimulationDuration() << "ms";
                CardiacSimulationArchiver<Problem>::Save(*p_problem, checkpoint_dir_basename + archive_foldername.str(), false);

                // Put a copy of the partial results aside (in a subdirectory of checkpoint_dir_basename).
                OutputFileHandler checkpoint_dir_basename_handler(checkpoint_dir_basename, false);
                OutputFileHandler partial_output_dir_handler(HeartConfig::Instance()->GetOutputDirectory(), false);
                if (PetscTools::AmMaster())
                {
                    EXPECT0(system, "cp -r " + partial_output_dir_handler.GetOutputDirectoryFullPath() + " " + checkpoint_dir_basename_handler.GetOutputDirectoryFullPath());
                }

                // Advance time stepper
                checkpoint_stepper.AdvanceOneTimeStep();
            }
        }
        else
        {
            p_problem->Solve();
        }

        delete p_problem;
    }

    /**
     * Run the simulation.
     * This method basically contains switches on the problem type and space dimension,
     * and calls CreateAndRun() to do the work.
     */
    void Run();

public:
    /**
     * Constructor
     *
     * This also runs the simulation immediately.
     *
     * @param parameterFileName  The name of the chaste parameters xml file to use to run a simulation (not mandatory since HeartConfig may be set by hand)
     */
    CardiacSimulation(std::string parameterFileName);
};

//
// Implementation must remain in this file (see comment by #include "CardiacSimulationArchiver.hpp").
//

CardiacSimulation::CardiacSimulation(std::string parameterFileName)
{
    // If we have been passed an XML file then parse the XML file, otherwise throw
    if (parameterFileName == "")
    {
        EXCEPTION("No XML file name given");
    }
    ReadParametersFromFile(parameterFileName);
    Run();
    HeartEventHandler::Headings();
    HeartEventHandler::Report();
}


void CardiacSimulation::ReadParametersFromFile(std::string parameterFileName)
{
    try
    {
        // Try the hardcoded schema location first
        HeartConfig::Instance()->SetUseFixedSchemaLocation(true);
        HeartConfig::Instance()->SetParametersFile(parameterFileName);
    }
    catch (Exception& e)
    {
        if (e.CheckShortMessageContains("Missing file parsing configuration") == "")
        {
            // Try using the schema location given in the XML
            HeartConfig::Instance()->Reset();
            HeartConfig::Instance()->SetUseFixedSchemaLocation(false);
            HeartConfig::Instance()->SetParametersFile(parameterFileName);
        }
        else
        {
            throw e;
        }
    }
}


void CardiacSimulation::Run()
{
    switch (HeartConfig::Instance()->GetDomain())
    {
        case cp::domain_type::Mono :
        {
            switch (HeartConfig::Instance()->GetSpaceDimension())
            {
                case 3:
                {
                    CreateAndRun<MonodomainProblem<3>,3>();
                    break;
                }

                case 2:
                {
                    CreateAndRun<MonodomainProblem<2>,2>();
                    break;
                }

                case 1:
                {
                    CreateAndRun<MonodomainProblem<1>,1>();
                    break;
                }
                default :
                    EXCEPTION("Monodomain space dimension not supported: should be 1, 2 or 3");
            }
            break;
        }

        case cp::domain_type::Bi :
        {
            switch (HeartConfig::Instance()->GetSpaceDimension())
            {
                case 3:
                {
                    CreateAndRun<BidomainProblem<3>,3>();
                    break;
                }
                case 2:
                {
                    CreateAndRun<BidomainProblem<2>,2>();
                    break;
                }
                case 1:
                {
                    CreateAndRun<BidomainProblem<1>,1>();
                    break;
                }
                default :
                {
                    EXCEPTION("Bidomain space dimension not supported: should be 1, 2 or 3");
                }
            }
            break;
        }
        default :
        {
            // If the domain is not set correctly then the XML parser will have picked it up before now!
            NEVER_REACHED;
        }
    }
}



#endif /*CARDIACSIMULATION_HPP_*/
