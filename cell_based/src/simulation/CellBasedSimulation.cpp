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
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <set>

#include "CellBasedSimulation.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"
#include "Version.hpp"

template<unsigned DIM>
CellBasedSimulation<DIM>::CellBasedSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                        std::vector<AbstractForce<DIM>*> forceCollection,
                                        bool deleteCellPopulationAndForceCollection,
                                        bool initialiseCells)
    : mEndTime(0.0),  // hours - this is set later on
      mrCellPopulation(rCellPopulation),
      mDeleteCellPopulation(deleteCellPopulationAndForceCollection),
      mAllocatedMemoryForForceCollection(deleteCellPopulationAndForceCollection),
      mInitialiseCells(initialiseCells),
      mNoBirth(false),
      mUpdateCellPopulation(true),
      mOutputDirectory(""),
      mSimulationOutputDirectory(mOutputDirectory),
      mNumBirths(0),
      mNumDeaths(0),
      mSamplingTimestepMultiple(1),
      mForceCollection(forceCollection),
      mOutputNodeVelocities(false)
{
    mpConfig = CellBasedConfig::Instance();

    // This line sets a random seed of 0 if it wasn't specified earlier.
    mpRandomGenerator = RandomNumberGenerator::Instance();

    // Different time steps are used for cell-centre and vertex-based simulations
    if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&rCellPopulation))
    {
        mDt = 1.0/120.0; // 30 seconds
    }
    else
    {
        mDt = 0.002; // smaller time step required for convergence/stability
    }

    if (mInitialiseCells)
    {
        mrCellPopulation.InitialiseCells();
    }
}


template<unsigned DIM>
CellBasedSimulation<DIM>::~CellBasedSimulation()
{
    if (mAllocatedMemoryForForceCollection)
    {
        for (typename std::vector<AbstractForce<DIM>*>::iterator force_iter = mForceCollection.begin();
             force_iter != mForceCollection.end();
             ++force_iter)
        {
            delete *force_iter;
        }
    }

    if (mDeleteCellPopulation)
    {
        for (typename std::vector<AbstractCellKiller<DIM>*>::iterator it=mCellKillers.begin();
             it != mCellKillers.end();
             ++it)
        {
            delete *it;
        }
        delete &mrCellPopulation;
    }
}


template<unsigned DIM>
unsigned CellBasedSimulation<DIM>::DoCellBirth()
{
    if (mNoBirth)
    {
        return 0;
    }

    unsigned num_births_this_step = 0;

    // Iterate over all cells, seeing if each one can be divided
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        // Check if this cell is ready to divide - if so create a new cell etc.
        if (cell_iter->GetAge() > 0.0)
        {
            if (cell_iter->ReadyToDivide())
            {
                // Create a new cell
                CellPtr p_new_cell = cell_iter->Divide();

                // Call method that determines how cell division occurs and returns a vector
                c_vector<double, DIM> new_location = CalculateCellDivisionVector(*cell_iter);

                // Add new cell to the cell population
                mrCellPopulation.AddCell(p_new_cell, new_location, *cell_iter);

                // Update counter
                num_births_this_step++;
            }
        }
    }
    return num_births_this_step;
}


template<unsigned DIM>
unsigned CellBasedSimulation<DIM>::DoCellRemoval()
{
    unsigned num_deaths_this_step = 0;

    // This labels cells as dead or apoptosing. It does not actually remove the cells,
    // mrCellPopulation.RemoveDeadCells() needs to be called for this.
    for (typename std::vector<AbstractCellKiller<DIM>*>::iterator killer_iter = mCellKillers.begin();
         killer_iter != mCellKillers.end();
         ++killer_iter)
    {
        (*killer_iter)->TestAndLabelCellsForApoptosisOrDeath();
    }

    num_deaths_this_step += mrCellPopulation.RemoveDeadCells();

    return num_deaths_this_step;
}


template<unsigned DIM>
const std::vector<AbstractForce<DIM>*> CellBasedSimulation<DIM>::rGetForceCollection() const
{
    return mForceCollection;
}


template<unsigned DIM>
c_vector<double, DIM> CellBasedSimulation<DIM>::CalculateCellDivisionVector(CellPtr pParentCell)
{
    /**
     * \todo Could remove this dynamic_cast by moving the code block below into
     * AbstractCentreBasedCellPopulation::AddCell(), allowing it to be overruled by
     * this method when overridden in subclasses. See also comment on #1093.
     */
    if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&mrCellPopulation))
    {
        // Location of parent and daughter cells
        c_vector<double, DIM> parent_coords = mrCellPopulation.GetLocationOfCellCentre(pParentCell);
        c_vector<double, DIM> daughter_coords;

        // Get separation parameter
        double separation = static_cast<AbstractCentreBasedCellPopulation<DIM>*>(&mrCellPopulation)->GetMeinekeDivisionSeparation();

        // Make a random direction vector of the required length
        c_vector<double, DIM> random_vector;

        /*
         * Pick a random direction and move the parent cell backwards by 0.5*separation
         * in that direction and return the position of the daughter cell 0.5*separation
         * forwards in that direction.
         */
        switch (DIM)
        {
            case 1:
            {
                double random_direction = -1.0 + 2.0*(RandomNumberGenerator::Instance()->ranf() < 0.5);

                random_vector(0) = 0.5*separation*random_direction;
                break;
            }
            case 2:
            {
                double random_angle = 2.0*M_PI*RandomNumberGenerator::Instance()->ranf();

                random_vector(0) = 0.5*separation*cos(random_angle);
                random_vector(1) = 0.5*separation*sin(random_angle);
                break;
            }
            case 3:
            {
                double random_zenith_angle = M_PI*RandomNumberGenerator::Instance()->ranf(); // phi
                double random_azimuth_angle = 2*M_PI*RandomNumberGenerator::Instance()->ranf(); // theta

                random_vector(0) = 0.5*separation*cos(random_azimuth_angle)*sin(random_zenith_angle);
                random_vector(1) = 0.5*separation*sin(random_azimuth_angle)*sin(random_zenith_angle);
                random_vector(2) = 0.5*separation*cos(random_zenith_angle);
                break;
            }
            default:
                // This can't happen
                NEVER_REACHED;
        }

        parent_coords = parent_coords - random_vector;
        daughter_coords = parent_coords + random_vector;

        // Set the parent to use this location
        ChastePoint<DIM> parent_coords_point(parent_coords);
        unsigned node_index = mrCellPopulation.GetLocationIndexUsingCell(pParentCell);
        mrCellPopulation.SetNode(node_index, parent_coords_point);

        return daughter_coords;
    }
    else
    {
        return zero_vector<double>(DIM);
    }
}


template<unsigned DIM>
void CellBasedSimulation<DIM>::UpdateNodePositions(const std::vector< c_vector<double, DIM> >& rNodeForces)
{
    unsigned num_nodes = mrCellPopulation.GetNumNodes();

    /*
     * Get the previous node positions (these may be needed when applying boundary conditions,
     * e.g. in the case of immotile cells)
     */
    std::vector<c_vector<double, DIM> > old_node_locations;
    old_node_locations.reserve(num_nodes);
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        old_node_locations[node_index] = mrCellPopulation.GetNode(node_index)->rGetLocation();
    }

    // Update node locations
    mrCellPopulation.UpdateNodeLocations(rNodeForces, mDt);

    // Apply any boundary conditions
    ApplyCellPopulationBoundaryConditions(old_node_locations);

    // Write node velocities to file if required
    if (mOutputNodeVelocities)
    {
        if (SimulationTime::Instance()->GetTimeStepsElapsed()%mSamplingTimestepMultiple==0)
        {
            *mpNodeVelocitiesFile << SimulationTime::Instance()->GetTime() << "\t";
            for (unsigned node_index=0; node_index<num_nodes; node_index++)
            {
                // We should never encounter deleted nodes due to where this method is called by Solve()
                assert(!mrCellPopulation.GetNode(node_index)->IsDeleted());

                // Check that results should be written for this node
                bool is_real_node = true;
                if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&mrCellPopulation))
                {
                    if (static_cast<AbstractCentreBasedCellPopulation<DIM>*>(&mrCellPopulation)->IsGhostNode(node_index))
                    {
                        // If this node is a ghost node then don't record its velocity
                        is_real_node = false;
                    }
                    else
                    {
                        // We should never encounter nodes associated with dead cells due to where this method is called by Solve()
                        assert(!mrCellPopulation.GetCellUsingLocationIndex(node_index)->IsDead());
                    }
                }

                // Write node data to file
                if (is_real_node)
                {
                    const c_vector<double,DIM>& position = mrCellPopulation.GetNode(node_index)->rGetLocation();
                    c_vector<double, 2> velocity = this->mDt * rNodeForces[node_index] / this->mrCellPopulation.GetDampingConstant(node_index);

                    *mpNodeVelocitiesFile << node_index  << " ";
                    for (unsigned i=0; i<DIM; i++)
                    {
                        *mpNodeVelocitiesFile << position[i] << " ";
                    }
                    for (unsigned i=0; i<DIM; i++)
                    {
                        *mpNodeVelocitiesFile << velocity[i] << " ";
                    }
                }
            }
            *mpNodeVelocitiesFile << "\n";
        }
    }
}


template<unsigned DIM>
void CellBasedSimulation<DIM>::SetDt(double dt)
{
    assert(dt > 0);
    mDt = dt;
}


template<unsigned DIM>
double CellBasedSimulation<DIM>::GetDt()
{
    return mDt;
}


template<unsigned DIM>
unsigned CellBasedSimulation<DIM>::GetNumBirths()
{
    return mNumBirths;
}


template<unsigned DIM>
unsigned CellBasedSimulation<DIM>::GetNumDeaths()
{
    return mNumDeaths;
}


template<unsigned DIM>
void CellBasedSimulation<DIM>::SetEndTime(double endTime)
{
    assert(endTime > 0);
    mEndTime = endTime;
}


template<unsigned DIM>
void CellBasedSimulation<DIM>::SetOutputDirectory(std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;
    mSimulationOutputDirectory = mOutputDirectory;
}


template<unsigned DIM>
std::string CellBasedSimulation<DIM>::GetOutputDirectory()
{
    return mOutputDirectory;
}

template<unsigned DIM>
void CellBasedSimulation<DIM>::SetSamplingTimestepMultiple(unsigned samplingTimestepMultiple)
{
    assert(samplingTimestepMultiple > 0);
    mSamplingTimestepMultiple = samplingTimestepMultiple;
}


template<unsigned DIM>
AbstractCellPopulation<DIM>& CellBasedSimulation<DIM>::rGetCellPopulation()
{
    return mrCellPopulation;
}


template<unsigned DIM>
const AbstractCellPopulation<DIM>& CellBasedSimulation<DIM>::rGetCellPopulation() const
{
    return mrCellPopulation;
}


template<unsigned DIM>
void CellBasedSimulation<DIM>::SetUpdateCellPopulationRule(bool updateCellPopulation)
{
    mUpdateCellPopulation = updateCellPopulation;
}


template<unsigned DIM>
void CellBasedSimulation<DIM>::SetNoBirth(bool noBirth)
{
    mNoBirth = noBirth;
}


template<unsigned DIM>
void CellBasedSimulation<DIM>::AddCellKiller(AbstractCellKiller<DIM>* pCellKiller)
{
    mCellKillers.push_back(pCellKiller);
}


template<unsigned DIM>
std::vector<double> CellBasedSimulation<DIM>::GetNodeLocation(const unsigned& rNodeIndex)
{
    std::vector<double> location;
    for (unsigned i=0; i<DIM; i++)
    {
        location.push_back(mrCellPopulation.GetNode(rNodeIndex)->rGetLocation()[i]);
    }
    return location;
}


template<unsigned DIM>
void CellBasedSimulation<DIM>::Solve()
{
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::EVERYTHING);
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::SETUP);

    // Set up the simulation time
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();

    unsigned num_time_steps = (unsigned) ((mEndTime-current_time)/mDt+0.5);

    if (current_time > 0) // use the reset function if necessary
    {
        p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
    }
    else
    {
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
    }

    if (mOutputDirectory=="")
    {
        EXCEPTION("OutputDirectory not set");
    }

    double time_now = p_simulation_time->GetTime();
    std::ostringstream time_string;
    time_string << time_now;

    std::string results_directory = mOutputDirectory +"/results_from_time_" + time_string.str();
    mSimulationOutputDirectory = results_directory;

    ///////////////////////////////////////////////////////////
    // Set up Simulation
    ///////////////////////////////////////////////////////////

    // Create output files for the visualizer
    OutputFileHandler output_file_handler(results_directory+"/", true);

    mrCellPopulation.CreateOutputFiles(results_directory+"/", false);

    mpVizSetupFile = output_file_handler.OpenOutputFile("results.vizsetup");

    if (mOutputNodeVelocities)
    {
        OutputFileHandler output_file_handler2(results_directory+"/", false);
        mpNodeVelocitiesFile = output_file_handler2.OpenOutputFile("nodevelocities.dat");
    }

    SetupSolve();

    // Age the cells to the correct time. Note that cells are created with
    // negative birth times so that some are initially almost ready to divide.
    LOG(1, "Setting up cells...");
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        // We don't use the result; this call is just to force the cells to age
        // to the current time running their cell cycle models to get there
        cell_iter->ReadyToDivide();
    }
    LOG(1, "\tdone\n");

    // Write initial conditions to file for the visualizer
    WriteVisualizerSetupFile();

    *mpVizSetupFile << std::flush;

    mrCellPopulation.WriteResultsToFiles();

    OutputSimulationSetup();

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::SETUP);

    // Initialise a vector of forces on node
    std::vector<c_vector<double, DIM> > forces(mrCellPopulation.GetNumNodes(), zero_vector<double>(DIM));

    /////////////////////////////////////////////////////////////////////
    // Main time loop
    /////////////////////////////////////////////////////////////////////

    while ((p_simulation_time->GetTimeStepsElapsed() < num_time_steps) && !(StoppingEventHasOccurred()) )
    {
        LOG(1, "--TIME = " << p_simulation_time->GetTime() << "\n");

        /**
         * This function calls:
         * DoCellRemoval()
         * DoCellBirth()
         * CellPopulation::Update()
         */
        UpdateCellPopulation();

        /////////////////////////
        // Calculate Forces
        /////////////////////////
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::FORCE);

        // First set all the forces to zero
        for (unsigned i=0; i<forces.size(); i++)
        {
             forces[i].clear();
        }

        // Then resize the std::vector if the number of cells has increased or decreased
        // (note this should be done after the above zeroing)
        unsigned num_nodes = mrCellPopulation.GetNumNodes();
        if (num_nodes != forces.size())
        {
            forces.resize(num_nodes, zero_vector<double>(DIM));
        }

        // Now add force contributions from each AbstractForce
        for (typename std::vector<AbstractForce<DIM>*>::iterator iter = mForceCollection.begin();
             iter != mForceCollection.end();
             ++iter)
        {
            (*iter)->AddForceContribution(forces, mrCellPopulation);
        }
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::FORCE);

        ////////////////////////////
        // Update node positions
        ////////////////////////////
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::POSITION);
        UpdateNodePositions(forces);
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::POSITION);

        //////////////////////////////////////////
        // PostSolve, which may be implemented by
        // child classes (e.g. to solve PDEs)
        //////////////////////////////////////////
        PostSolve();

        // Increment simulation time here, so results files look sensible
        p_simulation_time->IncrementTimeOneStep();

        ////////////////////////////
        // Output current results
        ////////////////////////////
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);

        // Write results to file
        if (p_simulation_time->GetTimeStepsElapsed()%mSamplingTimestepMultiple==0)
        {
            mrCellPopulation.WriteResultsToFiles();
        }

        CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);
    }

    LOG(1, "--END TIME = " << SimulationTime::Instance()->GetTime() << "\n");
    /* Carry out a final update so that cell population is coherent with new cell positions.
     * NB cell birth/death still need to be checked because they may be spatially-dependent.*/
    UpdateCellPopulation();

    AfterSolve();

    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);

    mrCellPopulation.CloseOutputFiles();

    if (mOutputNodeVelocities)
    {
        mpNodeVelocitiesFile->close();
    }

    *mpVizSetupFile << "Complete\n";
    mpVizSetupFile->close();

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::EVERYTHING);
}


template<unsigned DIM>
bool CellBasedSimulation<DIM>::StoppingEventHasOccurred()
{
    return false;
}


template<unsigned DIM>
void CellBasedSimulation<DIM>::UpdateCellPopulation()
{
    /////////////////////////
    // Remove dead cells
    /////////////////////////
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::DEATH);
    unsigned deaths_this_step = DoCellRemoval();
    mNumDeaths += deaths_this_step;
    LOG(1, "\tNum deaths = " << mNumDeaths << "\n");
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::DEATH);

    /////////////////////////
    // Divide cells
    /////////////////////////
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::BIRTH);
    unsigned births_this_step = DoCellBirth();
    mNumBirths += births_this_step;
    LOG(1, "\tNum births = " << mNumBirths << "\n");
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::BIRTH);

    // This allows NodeBasedCellPopulation::Update() to do the minimum amount of work
    bool births_or_death_occurred = ((births_this_step>0) || (deaths_this_step>0));

    ////////////////////////////
    // Update topology of cell population
    ////////////////////////////
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::UPDATECELLPOPULATION);
    if (mUpdateCellPopulation)
    {
        LOG(1, "\tUpdating cell population...");
        mrCellPopulation.Update(births_or_death_occurred);
        LOG(1, "\tdone.\n");
    }
    else if (births_or_death_occurred)
    {
        EXCEPTION("CellPopulation has had births or deaths but mUpdateCellPopulation is set to false, please set it to true.");
    }
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::UPDATECELLPOPULATION);
}

template<unsigned DIM>
void CellBasedSimulation<DIM>::OutputSimulationSetup()
{
    // Create Output file
    OutputFileHandler output_file_handler(this->mSimulationOutputDirectory + "/", false);

    out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");

    *parameter_file << "<Chaste>\n" ;

    // Output Chaste provenance information
    *parameter_file << "\t<ChasteInfo>\n" ;

    *parameter_file << "\t\t<VersionString>"<< ChasteBuildInfo::GetVersionString() << "</VersionString>\n";
    *parameter_file << "\t\t<IsWorkingCopyModified>"<< ChasteBuildInfo::IsWorkingCopyModified() << "</IsWorkingCopyModified>\n";
    *parameter_file << "\t\t<BuildInformation>"<< ChasteBuildInfo::GetBuildInformation() << "</BuildInformation>\n";
    *parameter_file << "\t\t<BuildTime>"<< ChasteBuildInfo::GetBuildTime() << "</BuildTime>\n";
    *parameter_file << "\t\t<CurrentTime>"<< ChasteBuildInfo::GetCurrentTime() << "</CurrentTime>\n";
    *parameter_file << "\t\t<BuilderUnameInfo>"<< ChasteBuildInfo::GetBuilderUnameInfo() << "</BuilderUnameInfo>\n";

    *parameter_file << "\t</ChasteInfo>\n";

    // Output CellBasedSimulation details
    ///\todo This should be independent of boost version (#1453)
	std::string simulation_type = "Should be simulation type here see #1453";
//	#if BOOST_VERSION >= 103700
//		simulation_type = GetIdentifier();
//	#endif

    *parameter_file <<  "\n\t<" << simulation_type << ">" "\n";
	OutputSimulationParameters(parameter_file);
    *parameter_file <<  "\t</" << simulation_type << ">" "\n";

    *parameter_file << "\n";

    // Output cell population details
    mrCellPopulation.OutputCellPopulationInfo(parameter_file);

    // Loop over forces
    *parameter_file << "\n\t<Forces>\n" ;
    for (typename std::vector<AbstractForce<DIM>*>::iterator iter = mForceCollection.begin();
                 iter != mForceCollection.end();
                 ++iter)
    {
    	// Output force details
    	(*iter)->OutputForceInfo(parameter_file);
    }
    *parameter_file << "\t</Forces>\n";

    // Loop over Cell Killers
	*parameter_file << "\n\t<CellKillers>\n" ;

	for (typename std::vector<AbstractCellKiller<DIM>*>::iterator iter = mCellKillers.begin();
				iter != mCellKillers.end();
				++iter)
	{
		// Output cell killer details
		(*iter)->OutputCellKillerInfo(parameter_file);
	}
	*parameter_file << "\t</CellKillers>\n";

	// Loop over ccm's


    // Output Extra Parameters from CellBasedConfig
    *parameter_file <<  "\n\t<CellBasedConfig>\n";

    CellBasedConfig* p_inst = CellBasedConfig::Instance();
    *parameter_file << "\t\t<SG2MDuration>"<< p_inst->GetSG2MDuration() << "</SG2MDuration>\n";
    *parameter_file << "\t\t<SDuration>"<< p_inst->GetSDuration() << "</SDuration>\n";
	*parameter_file << "\t\t<G2Duration>"<< p_inst->GetG2Duration() << "</G2Duration>\n";
	*parameter_file << "\t\t<Duration>"<< p_inst->GetMDuration() << "</Duration>\n";
	*parameter_file << "\t\t<StemCellG1Duration>"<< p_inst->GetStemCellG1Duration() << "</StemCellG1Duration>\n";
	*parameter_file << "\t\t<TransitCellG1Duration>"<< p_inst->GetTransitCellG1Duration() << "</TransitCellG1Duration>\n";
	*parameter_file << "\t\t<CryptLength>"<< p_inst->GetCryptLength() << "</CryptLength>\n";
	*parameter_file << "\t\t<CryptWidth>"<< p_inst->GetCryptWidth() << "</CryptWidth>\n";
	*parameter_file << "\t\t<MechanicsCutOffLength>"<< p_inst->GetMechanicsCutOffLength() << "</MechanicsCutOffLength>\n";
	*parameter_file << "\t\t<DampingConstantNormal>"<< p_inst->GetDampingConstantNormal() << "</DampingConstantNormal>\n";
	*parameter_file << "\t\t<DampingConstantMutant>"<< p_inst->GetDampingConstantMutant() << "</DampingConstantMutant>\n";
	//*parameter_file << "\t\t<CryptProjectionParameterA>"<< p_inst->GetCryptProjectionParameterA() << "</CryptProjectionParameterA>\n";
	//*parameter_file << "\t\t<CryptProjectionParameterB>"<< p_inst->GetCryptProjectionParameterB() << "</CryptProjectionParameterB>\n";

    *parameter_file <<  "\t</CellBasedConfig>\n";

    *parameter_file << "</Chaste>\n" ;

    parameter_file->close();
}

template<unsigned DIM>
bool CellBasedSimulation<DIM>::GetOutputNodeVelocities()
{
    return mOutputNodeVelocities;
}

template<unsigned DIM>
void CellBasedSimulation<DIM>::SetOutputNodeVelocities(bool outputNodeVelocities)
{
    mOutputNodeVelocities = outputNodeVelocities;
}

template<unsigned DIM>
void CellBasedSimulation<DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
	  *rParamsFile << "\t\t<Dt>"<< mDt << "</Dt>\n";
	  *rParamsFile << "\t\t<EndTime>"<< mEndTime << "</EndTime>\n";
	  *rParamsFile << "\t\t<SamplingTimestepMultiple>"<< mSamplingTimestepMultiple << "</SamplingTimestepMultiple>\n";
	  *rParamsFile << "\t\t<OutputNodeVelocities>"<< mOutputNodeVelocities << "</OutputNodeVelocities>\n";
}

////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class CellBasedSimulation<1>;
template class CellBasedSimulation<2>;
template class CellBasedSimulation<3>;
