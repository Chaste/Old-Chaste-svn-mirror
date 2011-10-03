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

#include "OffLatticeSimulationWithPdes.hpp"
#include "TrianglesMeshWriter.hpp"
#include "OutputFileHandler.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"
#include "CellBasedEventHandler.hpp"

template<unsigned DIM>
OffLatticeSimulationWithPdes<DIM>::OffLatticeSimulationWithPdes(AbstractCellPopulation<DIM>& rCellPopulation,
                                                                bool deleteCellPopulationInDestructor,
                                                                bool initialiseCells)
    : OffLatticeSimulation<DIM>(rCellPopulation,
                               deleteCellPopulationInDestructor,
                               initialiseCells),
      mpCellBasedPdeHandler(NULL)
{
}

template<unsigned DIM>
OffLatticeSimulationWithPdes<DIM>::~OffLatticeSimulationWithPdes()
{
}

template<unsigned DIM>
void OffLatticeSimulationWithPdes<DIM>::SetCellBasedPdeHandler(CellBasedPdeHandler<DIM>* pCellBasedPdeHandler)
{
    mpCellBasedPdeHandler = pCellBasedPdeHandler;
}

template<unsigned DIM>
CellBasedPdeHandler<DIM>* OffLatticeSimulationWithPdes<DIM>::GetCellBasedPdeHandler()
{
    return mpCellBasedPdeHandler;
}

template<unsigned DIM>
void OffLatticeSimulationWithPdes<DIM>::SetupSolve()
{
	if (mpCellBasedPdeHandler != NULL)
	{
        mpCellBasedPdeHandler->OpenResultsFiles(this->mSimulationOutputDirectory);
        *this->mpVizSetupFile << "PDE \n";
	}
}

template<unsigned DIM>
void OffLatticeSimulationWithPdes<DIM>::AfterSolve()
{
    if (mpCellBasedPdeHandler != NULL)
    {
        mpCellBasedPdeHandler->CloseResultsFiles();
    }
}

template<unsigned DIM>
void OffLatticeSimulationWithPdes<DIM>::PostSolve()
{
    if (mpCellBasedPdeHandler != NULL)
    {
    	CellBasedEventHandler::BeginEvent(CellBasedEventHandler::PDE);
        mpCellBasedPdeHandler->SolvePdeAndWriteResultsToFile(this->mSamplingTimestepMultiple);
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::PDE);
    }
}

template<unsigned DIM>
void OffLatticeSimulationWithPdes<DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    if (mpCellBasedPdeHandler != NULL)
    {
        mpCellBasedPdeHandler->OutputParameters(rParamsFile);
    }

    // Call method on direct parent class
    OffLatticeSimulation<DIM>::OutputSimulationParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class OffLatticeSimulationWithPdes<1>;
template class OffLatticeSimulationWithPdes<2>;
template class OffLatticeSimulationWithPdes<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OffLatticeSimulationWithPdes)
