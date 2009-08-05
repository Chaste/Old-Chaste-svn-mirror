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

#ifndef TESTUNSTRUCTUREDDATAREPLICATION_HPP_
#define TESTUNSTRUCTUREDDATAREPLICATION_HPP_

#include <cxxtest/TestSuite.h>
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"
#include "UnstructuredDataReplication.hpp"

class TestUnstructuredDataReplication : public CxxTest::TestSuite
{
public:
    void TestSimple()
    {
        // The unstructured data replication is filled with the following pattern.
        // First line shows slot numbering. Lines starting with Pi represent each
        // processor and the data it contributes (in the given slot).
        //
        //#slot 0   1  ...  p  p+1 p+2 ... 2*p      p is the number of processors
        //
        //  P0 000             000
        //  P1     111             111
        // ...         ...             ...
        //  Pp             ppp             ppp
        //
        // Once replicated every processor should own the following data.
        //
        // Rep 000 111 ... ppp 000 111 ... ppp

        unsigned num_procs = PetscTools::GetNumProcs();
        unsigned my_id = PetscTools::GetMyRank();

        unsigned num_slots = 2*num_procs;
        unsigned values_per_slot = 3;
        UnstructuredDataReplication buffer(num_slots, values_per_slot);

        std::vector<double> values(values_per_slot, my_id);

        buffer.SetSlot(my_id, values);
        buffer.SetSlot(num_procs + my_id, values);

        buffer.Replicate();

        for (unsigned slot_index=0; slot_index<num_slots; slot_index++)
        {
            for (unsigned value_index=0; value_index<values_per_slot; value_index++)
            {
                TS_ASSERT_DELTA(buffer.GetValue(slot_index,value_index), slot_index%num_procs, 1e-12);
            }
        }
    }

    void TestGetMethods()
    {
        // Tests both get methods return coherent data

        unsigned num_procs = PetscTools::GetNumProcs();
        unsigned my_id = PetscTools::GetMyRank();

        unsigned num_slots = 2*num_procs;
        unsigned values_per_slot = 3;
        UnstructuredDataReplication buffer(num_slots, values_per_slot);

        std::vector<double> values(values_per_slot, my_id);

        buffer.SetSlot(my_id, values);
        buffer.SetSlot(num_procs + my_id, values);

        buffer.Replicate();

        // Test we get the same results with the two different get methods
        std::vector<double> values_read(values_per_slot);
        for (unsigned slot_index=0; slot_index<num_slots; slot_index++)
        {
            buffer.GetSlot(slot_index, values_read);
            for (unsigned value_index=0; value_index<values_per_slot; value_index++)
            {
                TS_ASSERT_EQUALS(buffer.GetValue(slot_index,value_index), values_read[value_index]);
            }
        }

    }

    void TestMultipleSet()
    {
        EXIT_IF_SEQUENTIAL;

        // We don't allow more than one processor setting the same data
        unsigned num_procs = PetscTools::GetNumProcs();
        unsigned my_id = PetscTools::GetMyRank();

        unsigned num_slots = 2*num_procs;
        unsigned values_per_slot = 3;
        UnstructuredDataReplication buffer(num_slots, values_per_slot);

        std::vector<double> values(values_per_slot, my_id);

        buffer.SetSlot(my_id, values);
        buffer.SetSlot(num_procs, values);

        // All processors tried to set slot number num_procs
        TS_ASSERT_THROWS_ANYTHING(buffer.Replicate());
    }

    void TestUnsetSlot()
    {
        EXIT_IF_SEQUENTIAL;

        // We don't allow unset slots
        unsigned num_procs = PetscTools::GetNumProcs();
        unsigned my_id = PetscTools::GetMyRank();

        unsigned num_slots = 2*num_procs;
        unsigned values_per_slot = 3;
        UnstructuredDataReplication buffer(num_slots, values_per_slot);

        std::vector<double> values(values_per_slot, my_id);

        buffer.SetSlot(my_id, values);

        // Slots [num_procs, 2*num_procs[ were not set.
        TS_ASSERT_THROWS_ANYTHING(buffer.Replicate());
    }

    void TestInitialisation()
    {
        // Correctness of the algorithm relies on data being initialised to 0 when the object is constructed
        unsigned num_procs = PetscTools::GetNumProcs();
        unsigned num_slots = 2*num_procs;
        unsigned values_per_slot = 3;
        UnstructuredDataReplication buffer(num_slots, values_per_slot);

        for (unsigned index=0; index<buffer.mBufferLength; index++)
        {
            TS_ASSERT_EQUALS( buffer.mLocalData[index] , 0.0 );
        }

        for (unsigned index=0; index<buffer.mNumSlots; index++)
        {
            TS_ASSERT_EQUALS( buffer.mAccessPattern[index] , 0U );
        }

    }

    void TestMultipleReplications()
    {
        EXIT_IF_SEQUENTIAL;

        // Consecutive replications need the buffers to be reset in case the access pattern changes
        unsigned num_procs = PetscTools::GetNumProcs();
        unsigned my_id = PetscTools::GetMyRank();

        unsigned num_slots = 2*num_procs;
        unsigned values_per_slot = 3;
        UnstructuredDataReplication buffer(num_slots, values_per_slot);

        std::vector<double> values(values_per_slot, my_id);

        //  P0 000             000
        //  P1     111             111
        // ...         ...             ...
        //  Pp             ppp             ppp
        buffer.SetSlot(my_id, values);
        buffer.SetSlot(num_procs + my_id, values);

        buffer.Replicate();

        for (unsigned slot_index=0; slot_index<num_slots; slot_index++)
        {
            for (unsigned value_index=0; value_index<values_per_slot; value_index++)
            {
                TS_ASSERT_DELTA(buffer.GetValue(slot_index,value_index), slot_index%num_procs, 1e-12);
            }
        }

        //  P0     000              000
        //  P1         111              111
        // ...             ...              ...
        //  Pp ppp              ppp
        TS_ASSERT_THROWS_ANYTHING(buffer.SetSlot((my_id+1)%num_procs, values));

        buffer.Reset();

        buffer.SetSlot((my_id+1)%num_procs, values);
        buffer.SetSlot(num_procs + (my_id+1)%num_procs, values);
        TS_ASSERT_THROWS_NOTHING(buffer.Replicate());

        for (unsigned slot_index=0; slot_index<num_slots; slot_index++)
        {
            for (unsigned value_index=0; value_index<values_per_slot; value_index++)
            {
                TS_ASSERT_DELTA(buffer.GetValue(slot_index,value_index), (slot_index+num_procs-1)%num_procs, 1e-12);
            }
        }

    }
};

#endif /*TESTUNSTRUCTUREDDATAREPLICATION_HPP_*/
