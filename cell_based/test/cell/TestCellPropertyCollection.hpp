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
#ifndef TESTCELLPROPERTYCOLLECTION_HPP_
#define TESTCELLPROPERTYCOLLECTION_HPP_

#include "AbstractCellBasedTestSuite.hpp"

//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>

#include <boost/shared_ptr.hpp>

#include "CellPropertyCollection.hpp"
#include "AbstractCellProperty.hpp"

#include "AbstractCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"

//#include <boost/serialization/shared_ptr.hpp>

//#include "OutputFileHandler.hpp"

#define NEW_PROP(type, name) boost::shared_ptr<AbstractCellProperty> name(new type)

class TestCellPropertyCollection : public AbstractCellBasedTestSuite
{
public:
    void TestPropertyCollection() throw (Exception)
    {
        CellPropertyCollection collection;

        // Add some properties
        NEW_PROP(WildTypeCellMutationState, wt_mutation);
        collection.AddProperty(wt_mutation);
        NEW_PROP(ApcOneHitCellMutationState, apc1_mutation);
        collection.AddProperty(apc1_mutation);
        // Can't add the same *object* twice
        TS_ASSERT_THROWS_THIS(collection.AddProperty(wt_mutation),
                              "That property object is already in the collection.");
        NEW_PROP(WildTypeCellMutationState, wt_mutation_2);
        collection.AddProperty(wt_mutation_2);
        collection.RemoveProperty(wt_mutation_2);

        // Check the contents
        // ...by object
        TS_ASSERT(collection.HasProperty(wt_mutation));
        TS_ASSERT(collection.HasProperty(apc1_mutation));
        NEW_PROP(ApcOneHitCellMutationState, apc1_mutation_2);
        TS_ASSERT(!collection.HasProperty(apc1_mutation_2));
        // ...by type
        TS_ASSERT(collection.HasProperty<WildTypeCellMutationState>());
        TS_ASSERT(collection.HasProperty<ApcOneHitCellMutationState>());
        TS_ASSERT(!collection.HasProperty<ApcTwoHitCellMutationState>());
        // ..by subclass
        TS_ASSERT(collection.HasPropertyType<AbstractCellProperty>());
        TS_ASSERT(collection.HasPropertyType<AbstractCellMutationState>());
        //TS_ASSERT(!collection.HasProperty<AbstractCellMutationState>()); <-- This won't compile

        // Remove property
        collection.RemoveProperty<WildTypeCellMutationState>();
        TS_ASSERT(!collection.HasProperty<WildTypeCellMutationState>());
        collection.RemoveProperty(apc1_mutation);
        TS_ASSERT(!collection.HasProperty<ApcOneHitCellMutationState>());
        TS_ASSERT_THROWS_THIS(collection.RemoveProperty<WildTypeCellMutationState>(),
                              "Collection does not contain the given property type.");
        TS_ASSERT_THROWS_THIS(collection.RemoveProperty(apc1_mutation),
                              "Collection does not contain the given property.");

        // Get matching properties
        collection.AddProperty(wt_mutation);
        collection.AddProperty(apc1_mutation);
//        CellPropertyCollection mutations = collection.GetPropertiesType<AbstractCellMutationState>();
//        CellPropertyCollection wild_types = collection.GetProperties<WildTypeCellMutationState>();
    }

    // TODO: archive the collection
    void xTestArchiveCellProperties() throw (Exception)
    {
//        OutputFileHandler handler("archive", false);
//        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "mutation.arch";
//
//        // Archive a mutation state
//        {
//            AbstractCellMutationState* p_state = new ApcOneHitCellMutationState();
//            p_state->IncrementCellCount();
//
//            TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
//            TS_ASSERT_EQUALS(p_state->GetColour(), 3u);
//
//            // Create an output archive
//            std::ofstream ofs(archive_filename.c_str());
//            boost::archive::text_oarchive output_arch(ofs);
//
//            // Write the cell to the archive
//            const AbstractCellMutationState* const p_const_state = p_state;
//            output_arch << p_const_state;
//
//            delete p_state;
//        }
//
//        // Restore mutation state
//        {
//            AbstractCellMutationState* p_state;
//
//            // Restore the mutation state
//            std::ifstream ifs(archive_filename.c_str());
//            boost::archive::text_iarchive input_arch(ifs);
//
//            input_arch >> p_state;
//
//            TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
//            TS_ASSERT_EQUALS(p_state->GetColour(), 3u);
//
//            ApcOneHitCellMutationState* p_real_state = dynamic_cast<ApcOneHitCellMutationState*>(p_state);
//            TS_ASSERT(p_real_state != NULL);
//
//            // Tidy up
//            delete p_state;
//        }
    }
};

#endif /* TESTCELLPROPERTYCOLLECTION_HPP_ */
