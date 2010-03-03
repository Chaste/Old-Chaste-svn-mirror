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
#ifndef CELLWISEDATA_HPP_
#define CELLWISEDATA_HPP_

#include <climits> // work around boost bug
#include "ChasteSerialization.hpp"
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

    /** The single instance of the singleton object */
    static CellwiseData* mpInstance;

    /** A pointer to a Tissue so a cell's node can be found */
    MeshBasedTissue<DIM>* mpTissue;

    /** Allocated memory for mData object */
    bool mAllocatedMemory;

    /** Number of variables per node to be stored */
    unsigned mNumberOfVariables;

    /** Store of the data */
    std::vector<double> mData;

    /** Helper member storing constant data. Used in tests. */
    std::vector<double> mConstantDataForTesting;

    /** Helper member storing whether mConstantDataForTesting is used. */
    bool mUseConstantDataForTesting;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
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
     *  Protected constuctor. Not to be called, use Instance() instead.
     */
    CellwiseData();

public:

    /**
     *  Get an instance of the object
     */
    static CellwiseData* Instance();

    /**
     * Destructor.
     */
    virtual ~CellwiseData();

    /**
     *  Destroy the current instance. Should be called at the end of a
     *  simulation.
     */
    static void Destroy();

    /**
     * Get the value of CellwiseData for a given cell and variable number.
     *
     * @param rCell the cell
     * @param variableNumber the index of CellwiseData whose value is required (defaults to zero)
     *
     * @return the value of CellwiseData.
     */
    double GetValue(TissueCell& rCell, unsigned variableNumber=0);

    /**
     *  Set the value for a given node and variable number.
     *
     * @param value the value to set
     * @param pNode pointer to the Node
     * @param variableNumber the index of CellwiseData whose value is set (defaults to zero)
     */
    void SetValue(double value, Node<DIM>* pNode, unsigned variableNumber=0);

    /**
     *  Set the Tissue. Must be called before GetValue().
     *
     * @param rTissue reference to the Tissue
     */
    void SetTissue(MeshBasedTissue<DIM>& rTissue);

    /**
     *  @return reference to the Tissue.
     */
    MeshBasedTissue<DIM>& rGetTissue();

    /**
     *  Set the number of variables to be stored per cell. The constructor
     *  assumes 1 variable so only really needs to be called if numVars > 1.
     *
     * @param numNodes number of nodes in the tissue
     * @param numVars number of variables
     */
    void SetNumNodesAndVars(unsigned numNodes, unsigned numVars);

    /**
     *  Force the data to return given values for all cells (only for testing).
     *
     * @param values vector of CellwiseData values
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

#endif /*CELLWISEDATA_HPP_*/
