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

#ifndef TISSUESIMULATIONARCHIVER_HPP_
#define TISSUESIMULATIONARCHIVER_HPP_

#include <climits> // work around boost bug

// Must be included before any other serialisation headers
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <string>
#include <iostream>

#include "OutputFileHandler.hpp"
#include "SimulationTime.hpp"
#include "WntConcentration.hpp"
#include "CellwiseData.hpp"
#include "ArchiveLocationInfo.hpp"
#include "ArchiveOpener.hpp"

/**
 * TissueSimulationArchiver handles the checkpointing (saving and loading)
 * of all the various TissueSimulation objects. It has no explicit constructor
 * (just uses a default one) and no member variables.
 */
template<unsigned DIM, class SIM>
class TissueSimulationArchiver
{
public:

    /**
     * Loads a saved tissue simulation to run further.
     *
     * @param rArchiveDirectory  the name of the simulation to load
     *   (specified originally by simulation.SetOutputDirectory("wherever"); )
     * @param rTimeStamp  the time at which to load the simulation (this must
     *   be one of the times at which simulation.Save() was called)
     */
    static SIM* Load(const std::string& rArchiveDirectory, const double& rTimeStamp);

    /**
     * Saves the whole tissue simulation for restarting later.
     *
     * Puts it in the archive folder under the simulation's OutputDirectory,
     * in the file "tissue_sim_at_time_<SIMULATION TIME>.arch".
     * The mesh is written to files in the same folder.
     *
     * First archives simulation time (and other singletons, if used)
     * then the simulation itself.
     *
     * @param pSim pointer to the simulation
     */
    static void Save(SIM* pSim);
};


template<unsigned DIM, class SIM>
SIM* TissueSimulationArchiver<DIM, SIM>::Load(const std::string& rArchiveDirectory, const double& rTimeStamp)
{
    /**
     * Find the right archive (and mesh) to load.  The files are contained within
     * the 'archive' folder in rArchiveDirectory, with the archive itself called
     * 'tissue_sim_at_time_`rTimeStamp`.arch'.  The path to this file is returned.
     *
     * The path to the mesh is stored in ArchiveLocationInfo for use by the
     * Tissue de-serialization routines.
     */
    std::ostringstream time_stamp;
    time_stamp << rTimeStamp;
    std::string archive_filename = "tissue_sim_at_time_" + time_stamp.str() + ".arch";
    std::string mesh_filename = "mesh_" + time_stamp.str();
    ArchiveLocationInfo::SetMeshPathname(OutputFileHandler::GetChasteTestOutputDirectory()
                                         + rArchiveDirectory + "/archive/", mesh_filename);

    // Create an input archive
    ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(rArchiveDirectory + "/archive/", archive_filename);
    boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

    // Load any data that isn't the simulation itself, mainly singletons
    // - simulation time
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    assert(p_simulation_time->IsStartTimeSetUp());
    (*p_arch) & *p_simulation_time;

    // - Wnt concentration (if used)
    bool archive_wnt;
    (*p_arch) & archive_wnt;
    if (archive_wnt)
    {
        WntConcentration<DIM>* p_wnt = WntConcentration<DIM>::Instance();
        (*p_arch) & *p_wnt;
    }

    // - CellwiseData (if used)
    bool archive_cellwise_data;
    (*p_arch) & archive_cellwise_data;
    if (archive_cellwise_data)
    {
        CellwiseData<DIM>* p_cellwise_data = CellwiseData<DIM>::Instance();
        (*p_arch) & *p_cellwise_data;
    }

    // Load the simulation
    SIM* p_sim;
    (*p_arch) >> p_sim;
    return p_sim;
}

template<unsigned DIM, class SIM>
void TissueSimulationArchiver<DIM, SIM>::Save(SIM* pSim)
{
    // Get the simulation time as a string
    const SimulationTime* p_sim_time = SimulationTime::Instance();
    assert(p_sim_time->IsStartTimeSetUp());
    std::ostringstream time_stamp;
    time_stamp << p_sim_time->GetTime();

    // Set up folder and filename of archive
    std::string archive_directory = pSim->GetOutputDirectory() + "/archive/";
    std::string archive_filename = "tissue_sim_at_time_" + time_stamp.str() + ".arch";
    ArchiveLocationInfo::SetMeshFilename(std::string("mesh_") + time_stamp.str());

    // Create output archive
    ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_directory, archive_filename);
    boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

    // Save the simulation.  We save the time directly first to maintain its
    // singleton-ness on load.
    (*p_arch) << *p_sim_time;

    // Archive the Wnt concentration if it's used
    bool archive_wnt = WntConcentration<DIM>::Instance()->IsWntSetUp();
    (*p_arch) & archive_wnt;
    if (archive_wnt)
    {
        WntConcentration<DIM>* p_wnt = WntConcentration<DIM>::Instance();
        (*p_arch) & *p_wnt;
    }

    // Archive the CellwiseData if it's used
    bool archive_cellwise_data = CellwiseData<DIM>::Instance()->IsSetUp();
    (*p_arch) & archive_cellwise_data;
    if (archive_cellwise_data)
    {
        CellwiseData<DIM>* p_cellwise_data = CellwiseData<DIM>::Instance();
        (*p_arch) & *p_cellwise_data;
    }

    // Archive the simulation itself
    (*p_arch) & pSim; // const-ness would be a pain here
}

#endif /*TISSUESIMULATIONARCHIVER_HPP_*/
