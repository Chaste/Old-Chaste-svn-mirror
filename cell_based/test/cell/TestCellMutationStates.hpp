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
#ifndef TESTCELLMUTATIONSTATES_HPP_
#define TESTCELLMUTATIONSTATES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "ApoptoticCellMutationState.hpp"

#include "CellMutationStateRegistry.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractCellBasedTestSuite.hpp"
#include "OutputFileHandler.hpp"

class TestCellMutationStates : public AbstractCellBasedTestSuite
{
public:

    void TestCellMutationStateMethods() throw(Exception)
    {
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        TS_ASSERT_EQUALS(p_state->GetCellCount(), 0u);
        p_state->IncrementCellCount();
        TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
        p_state->DecrementCellCount();
        TS_ASSERT_EQUALS(p_state->GetCellCount(), 0u);
        TS_ASSERT_THROWS_THIS(p_state->DecrementCellCount(),
                "Cannot decrement cell count: no cells have this mutation state.");
        TS_ASSERT_EQUALS(p_state->GetColour(), 0u);

        TS_ASSERT_EQUALS(p_state->IsType<WildTypeCellMutationState>(), true);
        TS_ASSERT_EQUALS(p_state->IsType<ApcOneHitCellMutationState>(), false);

        boost::shared_ptr<AbstractCellMutationState> p_wt_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_apc2_state(new ApcTwoHitCellMutationState);
        TS_ASSERT(p_wt_state->IsSame(p_state.get()));
        TS_ASSERT(p_state->IsSame(p_wt_state));
        TS_ASSERT_EQUALS(p_wt_state->IsSame(p_apc2_state.get()), false);
        TS_ASSERT_EQUALS(p_apc2_state->IsSame(p_wt_state), false);

        // Check that const-ness doesn't matter
        TS_ASSERT(p_wt_state->IsType<const WildTypeCellMutationState>());
        const WildTypeCellMutationState const_wt_state;
        TS_ASSERT(p_wt_state->IsSame(&const_wt_state));
        TS_ASSERT(const_wt_state.IsSame(p_wt_state));
        TS_ASSERT(const_wt_state.IsSame(p_wt_state.get()));
    }

    void TestRegistry() throw(Exception)
    {
        boost::shared_ptr<AbstractCellMutationState> p_state1(
                CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellMutationState> p_state2(
                CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());
        TS_ASSERT(p_state1 == p_state2);
        TS_ASSERT_EQUALS(p_state1->GetCellCount(), 0u);
        p_state2->IncrementCellCount();
        TS_ASSERT_EQUALS(p_state1->GetCellCount(), 1u);

        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        TS_ASSERT_EQUALS(p_state->GetCellCount(), 0u);

        std::vector<boost::shared_ptr<AbstractCellMutationState> > states =
            CellMutationStateRegistry::Instance()->rGetAllMutationStates();
        TS_ASSERT_EQUALS(states.size(), 1u);
        TS_ASSERT(states[0] == p_state1);

        CellMutationStateRegistry::Instance()->Clear();
        states = CellMutationStateRegistry::Instance()->rGetAllMutationStates();
        TS_ASSERT_EQUALS(states.size(), 0u);

        // The taking-ownership functionality
        p_state1 = CellMutationStateRegistry::Instance()->Get<ApcTwoHitCellMutationState>();
        CellMutationStateRegistry::Instance()->Get<BetaCateninOneHitCellMutationState>();
        TS_ASSERT_EQUALS(CellMutationStateRegistry::Instance()->rGetAllMutationStates().size(), 2u);
        CellMutationStateRegistry* p_instance = CellMutationStateRegistry::Instance();
        CellMutationStateRegistry* p_registry = CellMutationStateRegistry::Instance()->TakeOwnership();
        TS_ASSERT_EQUALS(p_instance, p_registry);
        TS_ASSERT_EQUALS(p_registry->rGetAllMutationStates().size(), 2u);
        TS_ASSERT_DIFFERS(CellMutationStateRegistry::Instance(), p_registry);
        TS_ASSERT_EQUALS(CellMutationStateRegistry::Instance()->rGetAllMutationStates().size(), 0u);
        TS_ASSERT_EQUALS(p_registry->rGetAllMutationStates().size(), 2u);
        delete p_registry;
    }

    void TestMutationStateOrdering() throw(Exception)
    {
        // Ordering of mutation states matters for output at present, so we need a way of specifying it
        CellMutationStateRegistry* p_instance = CellMutationStateRegistry::Instance();
        p_instance->Clear();
        p_instance->Get<ApcOneHitCellMutationState>();
        std::vector<boost::shared_ptr<AbstractCellMutationState> > mutations;
        mutations.push_back(p_instance->Get<WildTypeCellMutationState>());
        mutations.push_back(p_instance->Get<ApcOneHitCellMutationState>());
        mutations.push_back(p_instance->Get<ApoptoticCellMutationState>());

        TS_ASSERT_EQUALS(p_instance->HasOrderingBeenSpecified(), false);
        p_instance->SpecifyOrdering(mutations);
        TS_ASSERT_EQUALS(p_instance->HasOrderingBeenSpecified(), true);

        TS_ASSERT_THROWS_THIS(p_instance->SpecifyOrdering(mutations),
                              "An ordering has already been specified.");

        std::vector<boost::shared_ptr<AbstractCellMutationState> > states = p_instance->rGetAllMutationStates();
        TS_ASSERT_EQUALS(states.size(), 3u);
        TS_ASSERT(states[0]->IsType<WildTypeCellMutationState>());
        TS_ASSERT(states[1]->IsType<ApcOneHitCellMutationState>());
        TS_ASSERT(states[2]->IsType<ApoptoticCellMutationState>());

        // The ordering must be complete
        TS_ASSERT_THROWS_THIS(p_instance->Get<BetaCateninOneHitCellMutationState>(),
                "Cannot add a new mutation state not specified in the ordering.");
        p_instance->Clear();
        p_instance->Get<BetaCateninOneHitCellMutationState>();
        TS_ASSERT_THROWS_THIS(p_instance->SpecifyOrdering(mutations),
                "The given ordering doesn't include all mutation states in the registry.");
    }

    void TestArchiveCellMutationState() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "mutation.arch";

        // Archive a mutation state
        {
            AbstractCellMutationState* p_state = new ApcOneHitCellMutationState();
            p_state->IncrementCellCount();

            TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_state->GetColour(), 3u);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            const AbstractCellMutationState* const p_const_state = p_state;
            output_arch << p_const_state;

            delete p_state;
        }

        // Restore mutation state
        {
            AbstractCellMutationState* p_state;

            // Restore the mutation state
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_state;

            TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_state->GetColour(), 3u);

            ApcOneHitCellMutationState* p_real_state = dynamic_cast<ApcOneHitCellMutationState*>(p_state);
            TS_ASSERT(p_real_state != NULL);

            // Tidy up
            delete p_state;
        }
    }

    void TestArchiveRegistry() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "mutation.arch";

        // Save
        {
            boost::shared_ptr<AbstractCellMutationState> p_state(CellMutationStateRegistry::Instance()->Get<WildTypeCellMutationState>());
            p_state->IncrementCellCount();

            const CellMutationStateRegistry* const p_registry = CellMutationStateRegistry::Instance()->TakeOwnership();

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            output_arch << p_registry;

            delete p_registry;
        }

        // Restore boost::shared_ptr to mutation state
        {
            // Initialize a mutation state
            CellMutationStateRegistry* p_registry;

            // Restore the mutation state
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_registry;

            boost::shared_ptr<AbstractCellMutationState> p_state = p_registry->Get<WildTypeCellMutationState>();

            TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_state->GetColour(), 0u);

            delete p_registry;
        }
    }
};

#endif /* TESTCELLMUTATIONSTATES_HPP_ */
