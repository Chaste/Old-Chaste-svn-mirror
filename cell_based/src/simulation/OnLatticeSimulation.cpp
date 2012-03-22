/*

Copyright (c) 2005-2012, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "OnLatticeSimulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"
#include "Version.hpp"
#include "ExecutableSupport.hpp"
#include "MultipleCaBasedCellPopulation.hpp"

template<unsigned DIM>
OnLatticeSimulation<DIM>::OnLatticeSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                              bool deleteCellPopulationInDestructor,
                                              bool initialiseCells)
    : AbstractCellBasedSimulation<DIM>(rCellPopulation,
                                       deleteCellPopulationInDestructor,
                                       initialiseCells),
      mOutputCellVelocities(false)
{
    if (!dynamic_cast<AbstractOnLatticeCellPopulation<DIM>*>(&rCellPopulation))
    {
        EXCEPTION("OnLatticeSimulations require a subclass of AbstractOnLatticeCellPopulation.");
    }

    this->mDt = 0.1; // 6 minutes
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::AddCaUpdateRule(boost::shared_ptr<AbstractCaUpdateRule<DIM> > pUpdateRule)
{
    if (dynamic_cast<CaBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
    {
        static_cast<CaBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->AddUpdateRule(pUpdateRule);
    }
}
template<unsigned DIM>
void OnLatticeSimulation<DIM>::AddMultipleCaUpdateRule(boost::shared_ptr<AbstractMultipleCaUpdateRule<DIM> > pUpdateRule)
{
    if (dynamic_cast<MultipleCaBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
    {
        static_cast<MultipleCaBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->AddUpdateRule(pUpdateRule);
    }
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::AddPottsUpdateRule(boost::shared_ptr<AbstractPottsUpdateRule<DIM> > pUpdateRule)
{
    if (dynamic_cast<PottsBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
    {
        static_cast<PottsBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->AddUpdateRule(pUpdateRule);
    }
}

template<unsigned DIM>
c_vector<double, DIM> OnLatticeSimulation<DIM>::CalculateCellDivisionVector(CellPtr pParentCell)
{
    ///\todo do something for Potts models here
    return zero_vector<double>(DIM);
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::WriteVisualizerSetupFile()
{
    if (dynamic_cast<PottsBasedCellPopulation<DIM>*>(&this->mrCellPopulation))
    {
       *this->mpVizSetupFile << "PottsSimulation\n";
    }
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::UpdateCellLocationsAndTopology()
{
    // Store whether we are sampling results at the current timestep
    SimulationTime* p_time = SimulationTime::Instance();
    bool at_sampling_timestep = (p_time->GetTimeStepsElapsed()%this->mSamplingTimestepMultiple == 0);

    /*
     * If required, store the current locations of cell centres. Note that we need to
     * use a std::map between cells and locations, rather than (say) a std::vector with
     * location indices corresponding to cells, since once we call UpdateCellLocations()
     * the location index of each cell may change. This is especially true in the case
     * of a CaBasedCellPopulation.
     */
    std::map<CellPtr, c_vector<double, DIM> > old_cell_locations;
    if (mOutputCellVelocities && at_sampling_timestep)
    {
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
             cell_iter != this->mrCellPopulation.End();
             ++cell_iter)
        {
            old_cell_locations[*cell_iter] = this->mrCellPopulation.GetLocationOfCellCentre(*cell_iter);
        }
    }

    // Update cell locations
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::POSITION);
    static_cast<AbstractOnLatticeCellPopulation<DIM>*>(&(this->mrCellPopulation))->UpdateCellLocations(this->mDt);
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::POSITION);

    // Now write cell velocities to file if required
    if (mOutputCellVelocities && at_sampling_timestep)
    {
        *mpCellVelocitiesFile << p_time->GetTime() << "\t";

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
             cell_iter != this->mrCellPopulation.End();
             ++cell_iter)
        {
            unsigned index = this->mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            const c_vector<double,DIM>& position = this->mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

            c_vector<double, DIM> velocity = (position - old_cell_locations[*cell_iter])/this->mDt;

            *mpCellVelocitiesFile << index  << " ";
            for (unsigned i=0; i<DIM; i++)
            {
                *mpCellVelocitiesFile << position[i] << " ";
            }

            for (unsigned i=0; i<DIM; i++)
            {
                *mpCellVelocitiesFile << velocity[i] << " ";
            }
        }
        *mpCellVelocitiesFile << "\n";
    }
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::SetupSolve()
{
    if (mOutputCellVelocities)
    {
        OutputFileHandler output_file_handler2(this->mSimulationOutputDirectory+"/", false);
        mpCellVelocitiesFile = output_file_handler2.OpenOutputFile("cellvelocities.dat");
    }
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::UpdateAtEndOfSolve()
{
    if (mOutputCellVelocities)
    {
        mpCellVelocitiesFile->close();
    }
}

template<unsigned DIM>
bool OnLatticeSimulation<DIM>::GetOutputCellVelocities()
{
    return mOutputCellVelocities;
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::SetOutputCellVelocities(bool outputCellVelocities)
{
    mOutputCellVelocities = outputCellVelocities;
}


template<unsigned DIM>
void OnLatticeSimulation<DIM>::UpdateCellPopulation()
{
    bool update_cell_population_this_timestep = true;
    if (dynamic_cast<CaBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
    {
        /*
         * If mInitialiseCells is false, then the simulation has been loaded from an archive.
         * In this case, we should not call UpdateCellPopulation() at the first time step. This is
         * because it will have already been called at the final time step prior to saving;
         * if we were to call it again now, then we would have introduced an extra call to
         * the random number generator compared to if we had not saved and loaded the simulation,
         * thus affecting results. This would be bad - we don't want saving and loading to have
         * any effect on the course of a simulation! See #1445.
         */
        if (!this->mInitialiseCells && (SimulationTime::Instance()->GetTimeStepsElapsed() == 0))
        {
            update_cell_population_this_timestep = false;
        }
    }

    if (update_cell_population_this_timestep)
    {
        AbstractCellBasedSimulation<DIM>::UpdateCellPopulation();
    }
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::OutputAdditionalSimulationSetup(out_stream& rParamsFile)
{
    // Loop over the collection of update rules and output info for each
    *rParamsFile << "\n\t<UpdateRules>\n";
    if (dynamic_cast<PottsBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
    {
        std::vector<boost::shared_ptr<AbstractPottsUpdateRule<DIM> > > collection =
            static_cast<PottsBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->rGetUpdateRuleCollection();

        for (typename std::vector<boost::shared_ptr<AbstractPottsUpdateRule<DIM> > >::iterator iter = collection.begin();
             iter != collection.end();
             ++iter)
        {
            (*iter)->OutputUpdateRuleInfo(rParamsFile);
        }
    }
    else if (dynamic_cast<CaBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
    {
        std::vector<boost::shared_ptr<AbstractCaUpdateRule<DIM> > > collection =
            static_cast<CaBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->rGetUpdateRuleCollection();

        for (typename std::vector<boost::shared_ptr<AbstractCaUpdateRule<DIM> > >::iterator iter = collection.begin();
             iter != collection.end();
             ++iter)
        {
            (*iter)->OutputUpdateRuleInfo(rParamsFile);
        }
    }
    else
    {
    	//\todo define the method for `MultipleCaBasedCellPopulation`
    }
    *rParamsFile << "\t</UpdateRules>\n";
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<OutputCellVelocities>" << mOutputCellVelocities << "</OutputCellVelocities>\n";

    // Call method on direct parent class
    AbstractCellBasedSimulation<DIM>::OutputSimulationParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class OnLatticeSimulation<1>;
template class OnLatticeSimulation<2>;
template class OnLatticeSimulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OnLatticeSimulation)
