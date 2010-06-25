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

#include "CellwiseData.hpp"


template<unsigned DIM>
CellwiseData<DIM>* CellwiseData<DIM>::mpInstance = NULL;


template<unsigned DIM>
CellwiseData<DIM>* CellwiseData<DIM>::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new CellwiseData<DIM>;
    }
    return mpInstance;
}


template<unsigned DIM>
CellwiseData<DIM>::CellwiseData()
 :  mpTissue(NULL),
    mAllocatedMemory(false),
    mConstantDataForTesting(0),
    mUseConstantDataForTesting(false)
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == NULL);
}


template<unsigned DIM>
CellwiseData<DIM>::~CellwiseData()
{
}


template<unsigned DIM>
void CellwiseData<DIM>::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = NULL;
    }
}


template<unsigned DIM>
double CellwiseData<DIM>::GetValue(TissueCell& rCell, unsigned variableNumber)
{
    // To test a cell and cell cycle models without a tissue
    if (mUseConstantDataForTesting)
    {
        return mConstantDataForTesting[variableNumber];
    }

    assert(IsSetUp());
    assert(mpTissue != NULL);
    assert(mAllocatedMemory);
    assert(variableNumber < mNumberOfVariables);

    unsigned location_index = mpTissue->GetLocationIndexUsingCell(rCell);
    unsigned vector_index = location_index*mNumberOfVariables + variableNumber;
    return mData[vector_index];
}


template<unsigned DIM>
void CellwiseData<DIM>::SetValue(double value, unsigned locationIndex, unsigned variableNumber)
{
    assert(IsSetUp());
    assert(variableNumber < mNumberOfVariables);

    unsigned vector_index = locationIndex*mNumberOfVariables + variableNumber;
    mData[vector_index] = value;
}


template<unsigned DIM>
void CellwiseData<DIM>::SetTissue(AbstractTissue<DIM>* pTissue)
{
    if (mAllocatedMemory == false)
    {
        EXCEPTION("SetTissue must be called after SetNumCellsAndVars()");
    }

    mpTissue = pTissue;
}


template<unsigned DIM>
AbstractTissue<DIM>& CellwiseData<DIM>::rGetTissue()
{
    return *mpTissue;
}


template<unsigned DIM>
void CellwiseData<DIM>::SetNumCellsAndVars(unsigned numCells, unsigned numberOfVariables)
{
    if (mpTissue!=NULL)
    {
        EXCEPTION("SetNumCellsAndVars() must be called before setting the Tissue (and after a Destroy)");
    }

    assert(numberOfVariables > 0);
    assert(mAllocatedMemory == false);

    mNumberOfVariables = numberOfVariables;
    mData.clear();
    mData.resize(numCells * mNumberOfVariables, 0.0);

    mAllocatedMemory = true;
}


template<unsigned DIM>
bool CellwiseData<DIM>::IsSetUp()
{
    return ((mAllocatedMemory) && (mpInstance!=NULL) && (mpTissue!=NULL));
}


template<unsigned DIM>
void CellwiseData<DIM>::ReallocateMemory()
{
    assert(mAllocatedMemory==true);
    assert(mpTissue!=NULL);

    unsigned num_cells = mpTissue->GetNumRealCells();
    if (mData.size() != num_cells*mNumberOfVariables)
    {
        mData.clear();
        mData.resize(num_cells*mNumberOfVariables, 0.0);
    }
}


template<unsigned DIM>
void CellwiseData<DIM>::SetConstantDataForTesting(std::vector<double>& rValues)
{
    mConstantDataForTesting = rValues;
    mUseConstantDataForTesting = true;
}


template<unsigned DIM>
unsigned CellwiseData<DIM>::GetNumVariables()
{
    return mNumberOfVariables;
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class CellwiseData<1>;
template class CellwiseData<2>;
template class CellwiseData<3>;
