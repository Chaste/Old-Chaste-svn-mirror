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

#include "ContactInhibitionOffLatticeSimulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CellwiseData.hpp"
#include "ContactInhibitionCellCycleModel.hpp"

template<unsigned DIM>
ContactInhibitionOffLatticeSimulation<DIM>::ContactInhibitionOffLatticeSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                                                  bool deleteCellPopulationInDestructor,
                                                                  bool initialiseCells)
    : OffLatticeSimulation<DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells)
{
}

template<unsigned DIM>
ContactInhibitionOffLatticeSimulation<DIM>::~ContactInhibitionOffLatticeSimulation()
{
}

template<unsigned DIM>
void ContactInhibitionOffLatticeSimulation<DIM>::PostSolve()
{
	assert(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)));

	// Loop over cells and set volume value in CellWiseData every 10 time step
	//\todo find a better way to make this efficient
	if (SimulationTime::Instance()->GetTimeStepsElapsed()%10==1)
	{
		MeshBasedCellPopulation<DIM>* p_static_cast_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&(this->mrCellPopulation));
	    this->mrCellPopulation.Update();

	    CellwiseData<DIM>::Instance()->ReallocateMemory();

		// Create and initialize a vector of cell volumes
		unsigned num_cells = this->mrCellPopulation.GetNumRealCells();
		std::vector<double> volume(num_cells);
		for (unsigned i=0; i<num_cells; i++)
		{
			volume[i] = 0;
		}

		p_static_cast_cell_population->CreateVoronoiTessellation();
		for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
			 cell_iter != this->mrCellPopulation.End();
			 ++cell_iter)
		{
			unsigned node_index = this->mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);
			volume[node_index] = p_static_cast_cell_population->GetVolumeOfVoronoiElement(node_index);
			CellwiseData<DIM>::Instance()->SetValue(volume[node_index], node_index, 0);
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class ContactInhibitionOffLatticeSimulation<1>;
template class ContactInhibitionOffLatticeSimulation<2>;
template class ContactInhibitionOffLatticeSimulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ContactInhibitionOffLatticeSimulation)
