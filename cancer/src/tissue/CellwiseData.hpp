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
#ifndef CELLWISEDATA_HPP_
#define CELLWISEDATA_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

#include "MeshBasedTissue.hpp"


/**
 *  A singleton object for storing data that certain cell cycle models
 *  need to know about, e.g. nutrient concentrations computed via some PDE
 *  for use in nutrient based cell cycle models.
 */
template<unsigned DIM>
class CellwiseData
{
    friend class TestCellwiseData;

private:

    /* the single instance of the singleton object */
    static CellwiseData* mpInstance;

    /* a reference to a Tissue so a cell's node can be found */
    MeshBasedTissue<DIM>* mpTissue;

    /*< allocated memory for mData object */
    bool mAllocatedMemory;

    /*< number of variables per node to be stored */
    unsigned mNumberOfVariables;

    /*< store of the data */
    std::vector<double> mData;

    std::vector<double> mConstantDataForTesting;
    bool mUseConstantDataForTesting;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mpTissue;
        archive & mAllocatedMemory;
        archive & mNumberOfVariables;
        archive & mData;
        archive & mConstantDataForTesting;
        archive & mUseConstantDataForTesting;
    }


protected:

    /**
     *  Protected constuctor. Not to be called, use Instance() instead
     */
    CellwiseData();


public:
    /**
     *  Get an instance of the object
     */
    static CellwiseData* Instance();

    virtual ~CellwiseData();

    /**
     *  Destroy the current instance. Should be called at the end of a
     *  simulation.
     */
    static void Destroy();

    /**
     *  Get the value for particular cell and given variable number (defaults
     *  to zero)
     */
    double GetValue(TissueCell* pCell, unsigned variableNumber=0);

    /**
     *  Set the value for particular node and given variable number (defaults
     *  to zero)
     */
    void SetValue(double value, Node<DIM>* pNode, unsigned variableNumber=0);

    /**
     *  Set the Tissue. Must be called before GetValue().
     */
    void SetTissue(MeshBasedTissue<DIM>& rTissue);

    /**
     *  Gets the tissue used in data.
     */
    MeshBasedTissue<DIM>& rGetTissue();

    /**
     *  Set the number of variables to be stored per cell. The constructor
     *  assumes 1 variable so only really needs to be called if num_vars > 1
     */
    void SetNumNodesAndVars(unsigned numNodes, unsigned numVars);

    /**
     *  Force the data to return given values for all cells (only for testing)
     */
    void SetConstantDataForTesting(std::vector<double> values);

    /**
     *  Is the instance in existence and fully set up
     */
    bool IsSetUp();

    /**
     *  Reallocate size of mData. Needed because of growth/death. Reallocates
     *  according to the number of nodes in the mesh in the Tissue member variable
     */
    void ReallocateMemory();

};

/** Pointer to the single instance */
template<unsigned DIM>
CellwiseData<DIM>* CellwiseData<DIM>::mpInstance = NULL;


/*
 * Return a pointer to the CellwiseData object.
 * The first time this is called, the object is created.
 */
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
double CellwiseData<DIM>::GetValue(TissueCell* pCell, unsigned variableNumber)
{
    // To test a cell and cell cycle models without a tissue
    if(mUseConstantDataForTesting)
    {
        return mConstantDataForTesting[variableNumber];
    }

    assert(IsSetUp());
    assert(mpTissue!=NULL);
    assert(mAllocatedMemory);
    assert(variableNumber < mNumberOfVariables);
    assert(pCell!=NULL);

    unsigned node_index = mpTissue->GetNodeCorrespondingToCell(*pCell)->GetIndex();
    unsigned vector_index = node_index*mNumberOfVariables + variableNumber;

    return mData[vector_index];
}

template<unsigned DIM>
void CellwiseData<DIM>::SetValue(double value, Node<DIM>* pNode, unsigned variableNumber)
{
    assert(IsSetUp());
    assert(variableNumber < mNumberOfVariables);
    unsigned vector_index = pNode->GetIndex()*mNumberOfVariables + variableNumber;
    mData[vector_index] = value;
}

template<unsigned DIM>
void CellwiseData<DIM>::SetTissue(MeshBasedTissue<DIM>& rTissue)
{
    if(mAllocatedMemory == false)
    {
        EXCEPTION("SetTissue must be called after SetNumNodesAndVars()");
    }

    mpTissue=&rTissue;
}

template<unsigned DIM>
MeshBasedTissue<DIM>& CellwiseData<DIM>::rGetTissue()
{
    return *mpTissue;
}


template<unsigned DIM>
void CellwiseData<DIM>::SetNumNodesAndVars(unsigned numNodes, unsigned numberOfVariables)
{
    if(mpTissue!=NULL)
    {
        EXCEPTION("SetNumNodesAndVars() must be called before setting the Tissue (and after a Destroy)");
    }

    assert(numberOfVariables>0);
    assert(mAllocatedMemory==false);

    mNumberOfVariables = numberOfVariables;
    mData.clear();
    mData.resize(numNodes * mNumberOfVariables, 0.0);

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

    unsigned num_nodes = mpTissue->rGetMesh().GetNumNodes();
    if (mData.size() != num_nodes*mNumberOfVariables)
    {
        mData.clear();
        mData.resize(num_nodes * mNumberOfVariables, 0.0);
    }
}

template<unsigned DIM>
void CellwiseData<DIM>::SetConstantDataForTesting(std::vector<double> values)
{
    mConstantDataForTesting = values;
    mUseConstantDataForTesting = true;
}

#endif /*CELLWISEDATA_HPP_*/
