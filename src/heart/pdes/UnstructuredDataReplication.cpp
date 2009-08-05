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

#include "UnstructuredDataReplication.hpp"
#include "Exception.hpp"

UnstructuredDataReplication::UnstructuredDataReplication(unsigned numSlots, unsigned valuesPerSlot):
    mNumSlots(numSlots),
    mValuesPerSlot(valuesPerSlot),
    mBufferLength(numSlots*valuesPerSlot),
    mReplicationDone(false)
{
    mReplicatedData = (double*) malloc (mBufferLength * sizeof(double));

    // The replication algorithm needs this two structures to be initialised to 0
    mLocalData =  (double*) calloc (mBufferLength, sizeof(double));
    mAccessPattern = (unsigned*) calloc (mNumSlots, sizeof(unsigned));
}

UnstructuredDataReplication::~UnstructuredDataReplication()
{
    free(mLocalData);
    free(mReplicatedData);
    free(mAccessPattern);
}

void UnstructuredDataReplication::SetSlot(unsigned slotNum, std::vector<double>& data)
{
    assert(slotNum < mNumSlots);
    assert(data.size() == mValuesPerSlot);

    if (mReplicationDone && mAccessPattern[slotNum]==0)
    {
        EXCEPTION("Access pattern to the structure has changed since last replication. Reset() must be called");
    }
    mAccessPattern[slotNum]=1;

    unsigned slot_offset = slotNum*mValuesPerSlot;

    for (unsigned index=0; index<mValuesPerSlot; index++)
    {
        mLocalData[slot_offset+index] = data[index];
    }

}


void UnstructuredDataReplication::Replicate()
{
    // Replicate the data
    MPI_Allreduce(mLocalData, mReplicatedData, mBufferLength, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

    // Check for multiple access to the same slot or unset slots
    unsigned global_access_pattern[mNumSlots];
    MPI_Allreduce(mAccessPattern, global_access_pattern, mNumSlots, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);

    for (unsigned slot_index=0; slot_index<mNumSlots; slot_index++)
    {
        switch (global_access_pattern[slot_index])
        {
            case 1:
                break;
            case 0:
                EXCEPTION("One or more slots were not filled by any processor");
            default:
                EXCEPTION("More than one processor wrote data in the same slot");
        }
    }

    // Everything went OK
    mReplicationDone = true;
}


void UnstructuredDataReplication::Reset()
{
    for(unsigned slot_index=0; slot_index<mNumSlots; slot_index++)
    {
        if (mAccessPattern[slot_index])
        {
            mAccessPattern[slot_index] = 0U;

            unsigned slot_offset = slot_index*mValuesPerSlot;
            for (unsigned index=0; index<mValuesPerSlot; index++)
            {
                mLocalData[slot_offset+index] = 0.0;
            }
        }
    }

    // Ready for the next replication
    mReplicationDone = false;
}

void UnstructuredDataReplication::GetSlot(unsigned slotNum, std::vector<double>& data)
{
    assert(mReplicationDone);
    assert(slotNum < mNumSlots);
    assert(data.size() == mValuesPerSlot);

    unsigned slot_offset = slotNum*mValuesPerSlot;

    for (unsigned index=0; index<mValuesPerSlot; index++)
    {
        data[index] = mReplicatedData[slot_offset+index];
    }
}


double UnstructuredDataReplication::GetValue(unsigned numSlot, unsigned numValue)
{
    assert(mReplicationDone);
    assert(numSlot < mNumSlots);
    assert(numValue < mValuesPerSlot);

    unsigned slot_offset = numSlot*mValuesPerSlot;

    return mReplicatedData[slot_offset + numValue];
}

