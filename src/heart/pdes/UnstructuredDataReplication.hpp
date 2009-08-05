/*

Copyright (C) University of Oxford, 2005-2009

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

#ifndef UNSTRUCTUREDDATAREPLICATION_HPP_
#define UNSTRUCTUREDDATAREPLICATION_HPP_

#include <vector>
#include <cassert>
#include "PetscTools.hpp"

/**
 * The unstructured data replication is filled with the following pattern.
 * First line shows slot numbering. Lines starting with Pi represent each
 * processor and the data it contributes (in the given slot).
 *
 * slot 0   1  ...  p  p+1 p+2 ... 2*p      p is the number of processors
 *
 *  P0 000             000
 *  P1     111             111
 * ...         ...             ...
 *  Pp             ppp             ppp
 *
 * Once replicated every processor should own the following data.
 *
 * Rep 000 111 ... ppp 000 111 ... ppp
 */
class UnstructuredDataReplication
{
    friend class TestUnstructuredDataReplication;

private:

    /** Total (global) number of slots (viz cells) */
    unsigned mNumSlots;
    /** Number of data values in each slot */
    unsigned mValuesPerSlot;
    /** Maximum length of a data == mNumSlots * mValuesPerSlot*/
    unsigned mBufferLength;

    /** Place (length mBufferLength) to store data before replication (with holes for data owned by remote processes)*/
    double* mLocalData;
    /** Place (length mBufferLength) to store data after replication*/
    double* mReplicatedData;
    
    /**Buffer of size mNumSlots which shows which slots are locally owned*/
    unsigned* mAccessPattern;
    
    /**Flag*/
    bool mReplicationDone;

public:

    /** Constructor
     * @param numSlots global number of slots
     * @param valuesPerSlot number of data values in each slot
     * 
     */
    UnstructuredDataReplication(unsigned numSlots, unsigned valuesPerSlot);

    /** Destructor
     */
    virtual ~UnstructuredDataReplication();

    /** Set a particular slot
     * - Mark it as owned by the local process
     * - fill it with data
     * @param slotNum
     * @param data
     */
    void SetSlot(unsigned slotNum, std::vector<double>& data);

    /**
     * Replicate the data
     * Uses All-to-all MPI_SUM reduction operation to share the data
     * in data buffers.
     * Uses All-to-all MPI_SUM reduction operation to ensure that each
     * slot is owned by exactly one process
     */
    void Replicate();

    /**
     * Clear the access pattern (local ownership)
     * and erase all local data.
     */
    void Reset();

    /** Get the data from a particular slot
     * @param slotNum
     * @param data
     */
    void GetSlot(unsigned slotNum, std::vector<double>& data);

    /** Get a single portion of the data from a particular slot
     * @param slotNum
     * @param numValue value within the slot
     * @return data
     */
    double GetValue(unsigned slotNum, unsigned numValue);
};

#endif /*UNSTRUCTUREDDATAREPLICATION_HPP_*/
