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
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <set>

#include "TissueSimulation.hpp"
#include "AbstractCellCentreBasedTissue.hpp"
#include "CancerEventHandler.hpp"
#include "LogFile.hpp"


template<unsigned DIM>
TissueSimulation<DIM>::TissueSimulation(AbstractTissue<DIM>& rTissue,
                                        std::vector<AbstractForce<DIM>*> forceCollection,
                                        bool deleteTissueAndForceCollection,
                                        bool initialiseCells)
  :  mrTissue(rTissue)
{
    #define COVERAGE_IGNORE
    assert(DIM==2 || DIM==3); // there are no instances of TissueSimulation<1>
    #undef COVERAGE_IGNORE

    mDeleteTissue = deleteTissueAndForceCollection;

    mInitialiseCells = initialiseCells;

    mpParams = CancerParameters::Instance();

    // This line sets a random seed of 0 if it wasn't specified earlier.
    mpRandomGenerator = RandomNumberGenerator::Instance();

    mDt = 1.0/120.0; // Timestep of 30 seconds (as per Meineke)
    mEndTime = 0.0; // hours - this is set later on.

    // Defaults
    mOutputDirectory = "";
    mSimulationOutputDirectory = mOutputDirectory;

    if (mrTissue.HasMesh())
    {
        mUpdateTissue = true;
    }
    else
    {
        mUpdateTissue = false;
    }

    mOutputCellMutationStates = false;
    mOutputCellAncestors = false;
    mOutputCellTypes = false;
    mOutputCellVariables = false;
    mOutputCellCyclePhases = false;
    mNoBirth = false;
    mNumBirths = 0;
    mNumDeaths = 0;
    mSamplingTimestepMultiple = 1;

    mAllocatedMemoryForForceCollection = deleteTissueAndForceCollection;

    mForceCollection = forceCollection;
    if (mInitialiseCells)
    {
        mrTissue.InitialiseCells();
    }
}


template<unsigned DIM>
TissueSimulation<DIM>::~TissueSimulation()
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

    if (mDeleteTissue)
    {
        for (typename std::vector<AbstractCellKiller<DIM>*>::iterator it=mCellKillers.begin();
             it != mCellKillers.end();
             ++it)
        {
            delete *it;
        }
        delete &mrTissue;
    }
}


template<unsigned DIM>
unsigned TissueSimulation<DIM>::DoCellBirth()
{
    /// \todo DoCellBirth() has not yet been tested with a vertex-based tissue - see #852

    if (mNoBirth)
    {
        return 0;
    }

    unsigned num_births_this_step = 0;

    // Iterate over all cells, seeing if each one can be divided
    for (typename AbstractTissue<DIM>::Iterator cell_iter = mrTissue.Begin();
         cell_iter != mrTissue.End();
         ++cell_iter)
    {
        // Check if this cell is ready to divide - if so create a new cell etc.
        if (cell_iter->GetAge() > 0.0)
        {
            if (cell_iter->ReadyToDivide())
            {
                // Create new cell
                TissueCell new_cell = cell_iter->Divide();

                // Add a new node to the mesh
                c_vector<double, DIM> new_location = CalculateDividingCellCentreLocations(&(*cell_iter));

                // Add a new cell to the tissue
                mrTissue.AddCell(new_cell, new_location, &(*cell_iter));

                // Update counter
                num_births_this_step++;
            }
        }
    }
    return num_births_this_step;
}


template<unsigned DIM>
unsigned TissueSimulation<DIM>::DoCellRemoval()
{
    /// \todo DoCellRemoval() has not yet been tested with a vertex-based tissue - see #853

    unsigned num_deaths_this_step = 0;

    // This labels cells as dead or apoptosing. It does not actually remove the cells,
    // tissue.RemoveDeadCells() needs to be called for this.
    for (typename std::vector<AbstractCellKiller<DIM>*>::iterator killer_iter = mCellKillers.begin();
         killer_iter != mCellKillers.end();
         ++killer_iter)
    {
        (*killer_iter)->TestAndLabelCellsForApoptosisOrDeath();
    }

    num_deaths_this_step += mrTissue.RemoveDeadCells();

    return num_deaths_this_step;
}


template<unsigned DIM>
const std::vector<AbstractForce<DIM>*> TissueSimulation<DIM>::rGetForceCollection() const
{
    return mForceCollection;
}


template<unsigned DIM>
c_vector<double, DIM> TissueSimulation<DIM>::CalculateDividingCellCentreLocations(TissueCell* pParentCell)
{
    /// \todo CalculateDividingCellCentreLocations() has not yet been tested with a vertex-based tissue - see #852

    double separation = CancerParameters::Instance()->GetDivisionSeparation();
    c_vector<double, DIM> parent_coords = dynamic_cast<AbstractCellCentreBasedTissue<DIM>*>(&mrTissue)->GetLocationOfCell(pParentCell);
    c_vector<double, DIM> daughter_coords;

    // Pick a random direction and move the parent cell backwards by 0.5*sep in that
    // direction and return the position of the daughter cell (0.5*sep forwards in the
    // random vector direction

    // Make a random direction vector of the required length
    c_vector<double, DIM> random_vector;

    if (DIM==2)
    {
        double random_angle = RandomNumberGenerator::Instance()->ranf();
        random_angle *= 2.0*M_PI;

        random_vector(0) = 0.5*separation*cos(random_angle);
        random_vector(1) = 0.5*separation*sin(random_angle);

        parent_coords = parent_coords - random_vector;
        daughter_coords = parent_coords + random_vector;
    }
    else if (DIM==3)
    {
        double random_zenith_angle = RandomNumberGenerator::Instance()->ranf(); // phi
        random_zenith_angle *= M_PI;
        double random_azimuth_angle = RandomNumberGenerator::Instance()->ranf(); // theta
        random_azimuth_angle *= 2*M_PI;

        random_vector(0) = 0.5*separation*cos(random_azimuth_angle)*sin(random_zenith_angle);
        random_vector(1) = 0.5*separation*sin(random_azimuth_angle)*sin(random_zenith_angle);
        random_vector(2) = 0.5*separation*cos(random_zenith_angle);

        daughter_coords = parent_coords + random_vector;
        parent_coords = parent_coords - random_vector;
    }

    // Set the parent to use this location
    ChastePoint<DIM> parent_coords_point(parent_coords);
    unsigned node_index = (static_cast<AbstractCellCentreBasedTissue<DIM>*>(&mrTissue))->GetNodeCorrespondingToCell(pParentCell)->GetIndex();
    mrTissue.SetNode(node_index, parent_coords_point);

    return daughter_coords;
}


template<unsigned DIM>
void TissueSimulation<DIM>::UpdateNodePositions(const std::vector< c_vector<double, DIM> >& rNodeForces)
{
    // Get the previous node positions (these may be needed
    // when applying boundary conditions, e.g. in the case
    // of immotile cells)
    std::vector<c_vector<double, DIM> > old_node_locations;
    old_node_locations.reserve(mrTissue.GetNumNodes());
    for (unsigned node_index=0; node_index<mrTissue.GetNumNodes(); node_index++)
    {
        old_node_locations[node_index] = mrTissue.GetNode(node_index)->rGetLocation();
    }

    // Update node locations
    mrTissue.UpdateNodeLocations(rNodeForces, mDt);

    // Apply any boundary conditions
    ApplyTissueBoundaryConditions(old_node_locations);
}


template<unsigned DIM>
void TissueSimulation<DIM>::SetDt(double dt)
{
    assert(dt > 0);
    mDt = dt;
}


template<unsigned DIM>
double TissueSimulation<DIM>::GetDt()
{
    return mDt;
}


template<unsigned DIM>
unsigned TissueSimulation<DIM>::GetNumBirths()
{
    return mNumBirths;
}


template<unsigned DIM>
unsigned TissueSimulation<DIM>::GetNumDeaths()
{
    return mNumDeaths;
}


template<unsigned DIM>
void TissueSimulation<DIM>::SetEndTime(double endTime)
{
    assert(endTime > 0);
    mEndTime = endTime;
}


template<unsigned DIM>
void TissueSimulation<DIM>::SetOutputDirectory(std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;
    mSimulationOutputDirectory = mOutputDirectory;
}


template<unsigned DIM>
std::string TissueSimulation<DIM>::GetOutputDirectory()
{
    return mOutputDirectory;
}


template<unsigned DIM>
void TissueSimulation<DIM>::SetSamplingTimestepMultiple(unsigned samplingTimestepMultiple)
{
    assert(samplingTimestepMultiple > 0);
    mSamplingTimestepMultiple = samplingTimestepMultiple;
}


template<unsigned DIM>
AbstractTissue<DIM>& TissueSimulation<DIM>::rGetTissue()
{
    return mrTissue;
}


template<unsigned DIM>
const AbstractTissue<DIM>& TissueSimulation<DIM>::rGetTissue() const
{
    return mrTissue;
}


template<unsigned DIM>
void TissueSimulation<DIM>::SetUpdateTissueRule(bool updateTissue)
{
    mUpdateTissue = updateTissue;
}


template<unsigned DIM>
void TissueSimulation<DIM>::SetNoBirth(bool noBirth)
{
    mNoBirth = noBirth;
}


template<unsigned DIM>
void TissueSimulation<DIM>::SetOutputCellMutationStates(bool outputCellMutationStates)
{
    mOutputCellMutationStates = outputCellMutationStates;
}


template<unsigned DIM>
void TissueSimulation<DIM>::SetOutputCellTypes(bool outputCellTypes)
{
    mOutputCellTypes = outputCellTypes;
}


template<unsigned DIM>
void TissueSimulation<DIM>::SetOutputCellAncestors(bool outputCellAncestors)
{
    mOutputCellAncestors = outputCellAncestors;
}


template<unsigned DIM>
void TissueSimulation<DIM>::SetOutputCellVariables(bool outputCellVariables)
{
    mOutputCellVariables = outputCellVariables;
}


template<unsigned DIM>
void TissueSimulation<DIM>::SetOutputCellCyclePhases(bool outputCellCyclePhases)
{
    mOutputCellCyclePhases = outputCellCyclePhases;
}


template<unsigned DIM>
void TissueSimulation<DIM>::AddCellKiller(AbstractCellKiller<DIM>* pCellKiller)
{
    mCellKillers.push_back(pCellKiller);
}


template<unsigned DIM>
std::vector<double> TissueSimulation<DIM>::GetNodeLocation(const unsigned& rNodeIndex)
{
    std::vector<double> location;
    for (unsigned i=0; i<DIM; i++)
    {
        location.push_back( mrTissue.GetNode(rNodeIndex)->rGetLocation()[i] );
    }
    return location;
}


template<unsigned DIM>
void TissueSimulation<DIM>::Solve()
{
    CancerEventHandler::BeginEvent(CANCER_EVERYTHING);
    CancerEventHandler::BeginEvent(SETUP);

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

    mrTissue.CreateOutputFiles(results_directory+"/",
                               false,
                               mOutputCellMutationStates,
                               mOutputCellTypes,
                               mOutputCellVariables,
                               mOutputCellCyclePhases,
                               mOutputCellAncestors);

    mpSetupFile = output_file_handler.OpenOutputFile("results.vizsetup");

    SetupSolve();

    // Age the cells to the correct time. Note that cells are created with
    // negative birth times so that some are initially almost ready to divide.
    LOG(1, "Setting up cells...");
    for (typename AbstractTissue<DIM>::Iterator cell_iter = mrTissue.Begin();
         cell_iter != mrTissue.End();
         ++cell_iter)
    {
        // We don't use the result; this call is just to force the cells to age
        // to the current time running their cell cycle models to get there
        cell_iter->ReadyToDivide();
    }
    LOG(1, "\tdone\n");

    // Write initial conditions to file for the visualizer
    if (DIM==2)
    {
        WriteVisualizerSetupFile();
    }
    *mpSetupFile << std::flush;

    mrTissue.WriteResultsToFiles(mOutputCellMutationStates,
                                 mOutputCellTypes,
                                 mOutputCellVariables,
                                 mOutputCellCyclePhases,
                                 mOutputCellAncestors);

    CancerEventHandler::EndEvent(SETUP);

    // Initialise a vector of forces on node
    std::vector<c_vector<double, DIM> > forces(mrTissue.GetNumNodes(),zero_vector<double>(DIM));

    /////////////////////////////////////////////////////////////////////
    // Main time loop
    /////////////////////////////////////////////////////////////////////

    while ((p_simulation_time->GetTimeStepsElapsed() < num_time_steps) && !(StoppingEventHasOccurred()) )
    {
        LOG(1, "--TIME = " << p_simulation_time->GetTime() << "\n");

        /////////////////////////
        // Remove dead cells
        /////////////////////////
        CancerEventHandler::BeginEvent(DEATH);
        unsigned deaths_this_step = DoCellRemoval();
        mNumDeaths += deaths_this_step;
        LOG(1, "\tNum deaths = " << mNumDeaths << "\n");
        CancerEventHandler::EndEvent(DEATH);

        /////////////////////////
        // Divide cells
        /////////////////////////
        CancerEventHandler::BeginEvent(BIRTH);
        unsigned births_this_step = DoCellBirth();
        mNumBirths += births_this_step;
        LOG(1, "\tNum births = " << mNumBirths << "\n");
        CancerEventHandler::EndEvent(BIRTH);

        ////////////////////////////
        // Update topology of tissue
        ////////////////////////////

        /**
         * If the tissue has a mesh, then we currently must call the Update()
         * method at each time step. Otherwise, we only need to call Update()
         * after there has been any cell birth or cell death.
         */
        if (mrTissue.HasMesh())
        {
            //This assertion is not necessarily true
            //assert(mUpdateTissue);
            //See TestCryptSimulation2dNightly::Test2DSpringSystem where the
            //default value is over-written
        }
        else
        {
            mUpdateTissue = false;
            if ( (births_this_step>0) || (deaths_this_step>0) )
            {
                mUpdateTissue = true;
            }
        }

        // Update the topology of the tissue (and tessellate if needed)
        CancerEventHandler::BeginEvent(UPDATE);
        if (mUpdateTissue)
        {
            LOG(1, "\tUpdating tissue...");
            mrTissue.Update();
            LOG(1, "\tdone.\n");
        }
        CancerEventHandler::EndEvent(UPDATE);

        /////////////////////////
        // Calculate Forces
        /////////////////////////
        CancerEventHandler::BeginEvent(FORCE);

        // First set all the forces to zero
        for (unsigned i=0; i<forces.size(); i++)
        {
             forces[i].clear();
        }

        // Then resize the std::vector if the number of cells has increased or decreased
        // (note this should be done after the above zeroing)
        if (mrTissue.GetNumNodes()!=forces.size())
        {
            forces.resize(mrTissue.GetNumNodes(),zero_vector<double>(DIM));
        }

        // Now add force contributions from each AbstractForce
        for (typename std::vector<AbstractForce<DIM>*>::iterator iter = mForceCollection.begin();
             iter !=mForceCollection.end();
             iter++)
        {
            (*iter)->AddForceContribution(forces, mrTissue);
        }
        CancerEventHandler::EndEvent(FORCE);

        ////////////////////////////
        // Update node positions
        ////////////////////////////
        CancerEventHandler::BeginEvent(POSITION);
        UpdateNodePositions(forces);
        CancerEventHandler::EndEvent(POSITION);

        //////////////////////////////////////////
        // PostSolve, which may be implemented by
        // child classes (eg to solve nutrient pdes)
        //////////////////////////////////////////
        PostSolve();

        // Increment simulation time here, so results files look sensible
        p_simulation_time->IncrementTimeOneStep();


        ////////////////////////////
        // Output current results
        ////////////////////////////
        CancerEventHandler::BeginEvent(OUTPUT);

        // Write results to file
        if (p_simulation_time->GetTimeStepsElapsed()%mSamplingTimestepMultiple==0)
        {
            mrTissue.WriteResultsToFiles(mOutputCellMutationStates,
                                         mOutputCellTypes,
                                         mOutputCellVariables,
                                         mOutputCellCyclePhases,
                                         mOutputCellAncestors);
        }

        CancerEventHandler::EndEvent(OUTPUT);

    }

    AfterSolve();

    CancerEventHandler::BeginEvent(OUTPUT);
    mrTissue.CloseOutputFiles(mOutputCellMutationStates,
                              mOutputCellTypes,
                              mOutputCellVariables,
                              mOutputCellCyclePhases,
                              mOutputCellAncestors);

    *mpSetupFile << "Complete\n";
    mpSetupFile->close();

    CancerEventHandler::EndEvent(OUTPUT);

    CancerEventHandler::EndEvent(CANCER_EVERYTHING);
}


template<unsigned DIM>
void TissueSimulation<DIM>::AfterSolve()
{
    LOG(1, "--TIME = " << SimulationTime::Instance()->GetTime() << "\n");

    // Remove dead cells then implement cell birth
    CancerEventHandler::BeginEvent(DEATH);
    mNumDeaths += DoCellRemoval();
    LOG(1, "\tNum deaths = " << mNumDeaths << "\n");
    CancerEventHandler::EndEvent(DEATH);

    CancerEventHandler::BeginEvent(BIRTH);
    mNumBirths += DoCellBirth();
    LOG(1, "\tNum births = " << mNumBirths << "\n");
    CancerEventHandler::EndEvent(BIRTH);

    // Carry out a final tissue update if necessary
    if (mrTissue.HasMesh())
    {
        //This assertion is not necessarily true
        //assert(mUpdateTissue);
        //See TestCryptSimulation2dNightly::Test2DSpringSystem where the
        //default value is over-written
    }
    else
    {
        mUpdateTissue = false;
        if ( (mNumBirths>0) || (mNumDeaths>0) )
        {
            mUpdateTissue = true;
        }
    }

    if (mUpdateTissue)
    {
        LOG(1, "\tUpdating tissue...");
        mrTissue.Update();
        LOG(1, "\tdone.\n");
    }
}


template<unsigned DIM>
bool TissueSimulation<DIM>::StoppingEventHasOccurred()
{
    return false;
}


template<unsigned DIM>
c_vector<unsigned, NUM_CELL_MUTATION_STATES> TissueSimulation<DIM>::GetCellMutationStateCount()
{
    if (!mOutputCellMutationStates)
    {
        EXCEPTION("Call simulator.SetOutputCellMutationStates before using this function");
    }
    return mrTissue.GetCellMutationStateCount();
}


template<unsigned DIM>
c_vector<unsigned, NUM_CELL_TYPES> TissueSimulation<DIM>::GetCellTypeCount()
{
    if (!mOutputCellTypes)
    {
        EXCEPTION("Call simulator.SetOutputCellTypes() before using this function");
    }
    return mrTissue.GetCellTypeCount();
}


template<unsigned DIM>
c_vector<unsigned, 5> TissueSimulation<DIM>::GetCellCyclePhaseCount()
{
    if (!mOutputCellCyclePhases)
    {
        EXCEPTION("Call simulator.SetOutputCellCyclePhases() before using this function");
    }
    return mrTissue.GetCellCyclePhaseCount();
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class TissueSimulation<1>;
template class TissueSimulation<2>;
template class TissueSimulation<3>;
