/*

Copyright (C) University of Oxford, 2005-2011

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

#include "CardiacSimulationArchiver.hpp"  // Must go first
#include "CardiacSimulation.hpp"


boost::shared_ptr<AbstractUntemplatedCardiacProblem> CardiacSimulation::GetSavedProblem()
{
    return mSavedProblem;
}

std::string CardiacSimulation::BoolToString(bool yesNo)
{
    std::string result;
    if (yesNo)
    {
        result = "yes";
    }
    else
    {
        result = "no";
    }
    return result;
}

void CardiacSimulation::CreateResumeXmlFile(const std::string& rOutputDirectory, const std::string& rArchiveDirectory)
{
    OutputFileHandler handler(rOutputDirectory, false);
    if (PetscTools::AmMaster())
    {
        out_stream p_file = handler.OpenOutputFile("ResumeParameters.xml");
        (*p_file) << "<?xml version='1.0' encoding='UTF-8'?>" << std::endl;
        (*p_file) << "<ChasteParameters xmlns='https://chaste.comlab.ox.ac.uk/nss/parameters/2_4' "
                  << "xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance' "
                  << "xsi:schemaLocation='https://chaste.comlab.ox.ac.uk/nss/parameters/2_4 ChasteParameters_2_4.xsd'>" << std::endl;
        (*p_file) << std::endl;
        (*p_file) << "    <ResumeSimulation>" << std::endl;
        (*p_file) << "        <ArchiveDirectory relative_to='this_file'>" << rArchiveDirectory << "</ArchiveDirectory>" << std::endl;
        (*p_file) << "        <SpaceDimension>" << HeartConfig::Instance()->GetSpaceDimension() << "</SpaceDimension>" << std::endl;
        (*p_file) << "        <SimulationDuration unit='ms'>0.0</SimulationDuration> <!-- Edit with new simulation duration. Please "
                  << "note that the simulation does not restart at t=0 but at the time where the checkpoint was created.-->" << std::endl;
        (*p_file) << "        <Domain>" << HeartConfig::Instance()->GetDomain() << "</Domain>" << std::endl;
        (*p_file) << "        <CheckpointSimulation timestep='" << HeartConfig::Instance()->GetCheckpointTimestep()
                  << "' unit='ms' max_checkpoints_on_disk='" << HeartConfig::Instance()->GetMaxCheckpointsOnDisk()
                  << "'/> <!-- This is optional; if not given, the loaded simulation will NOT itself be checkpointed -->" << std::endl;
        (*p_file) << "        <OutputVisualizer meshalyzer='" << BoolToString(HeartConfig::Instance()->GetVisualizeWithMeshalyzer())
                  << "' vtk='" << BoolToString(HeartConfig::Instance()->GetVisualizeWithVtk())
                  << "' parallel_vtk='" << BoolToString(HeartConfig::Instance()->GetVisualizeWithParallelVtk())
                  << "' cmgui='" << BoolToString(HeartConfig::Instance()->GetVisualizeWithCmgui()) << "'/>" << std::endl;
        (*p_file) << "    </ResumeSimulation>" << std::endl;
        (*p_file) << std::endl;
        (*p_file) << "    <!-- These elements must exist, but their contents are ignored -->" << std::endl;
        (*p_file) << "    <Physiological/>" << std::endl;
        (*p_file) << "    <Numerical/>" << std::endl;
        (*p_file) << "</ChasteParameters>" << std::endl;
        p_file->close();
    }
    HeartConfig::Instance()->CopySchema(handler.GetOutputDirectoryFullPath());
}

CardiacSimulation::CardiacSimulation(std::string parameterFileName,
                                     bool writeProvenanceInfo,
                                     bool saveProblemInstance)
    : mSaveProblemInstance(saveProblemInstance)
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
    if (writeProvenanceInfo)
    {
        ExecutableSupport::SetOutputDirectory(HeartConfig::Instance()->GetOutputDirectory());
        ExecutableSupport::WriteProvenanceInfoFile();
        ExecutableSupport::WriteMachineInfoFile("machine_info");
    }
}

void CardiacSimulation::ReadParametersFromFile(std::string parameterFileName)
{
    // Ensure the singleton is in a clean state
    HeartConfig::Reset();
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
            HeartConfig::Reset();
            HeartConfig::Instance()->SetUseFixedSchemaLocation(false);
            HeartConfig::Instance()->SetParametersFile(parameterFileName);
        }
        else
        {
            throw e;
        }
    }
}


#define DOMAIN_CASE(VALUE, CLASS, DIM)    \
    case VALUE:                           \
    {                                     \
        CreateAndRun<CLASS<DIM>, DIM>();  \
        break;                            \
    }

#define DOMAIN_SWITCH(DIM)                                                     \
    switch (HeartConfig::Instance()->GetDomain())                              \
    {                                                                          \
        DOMAIN_CASE(cp::domain_type::Mono, MonodomainProblem, DIM)             \
        DOMAIN_CASE(cp::domain_type::Bi, BidomainProblem, DIM)                 \
        DOMAIN_CASE(cp::domain_type::BiWithBath, BidomainWithBathProblem, DIM) \
        default:                                                               \
            NEVER_REACHED;                                                     \
    }                                                                          \
    break
// Note that if the domain is not set correctly then the XML parser will have picked it up before now!
// Missing semi-colon after break so we can put it after the macro call.

void CardiacSimulation::Run()
{
    switch (HeartConfig::Instance()->GetSpaceDimension())
    {
        case 3:
        {
            DOMAIN_SWITCH(3);
        }
        case 2:
        {
            DOMAIN_SWITCH(2);
        }
        case 1:
        {
            DOMAIN_SWITCH(1);
        }
        default:
            // We could check for this too in the XML Schema...
            EXCEPTION("Space dimension not supported: should be 1, 2 or 3");
    }
}

// These aren't needed externally
#undef DOMAIN_SWITCH
#undef DOMAIN_CASE

