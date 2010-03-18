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

#ifndef CELLMUTATIONSTATEREGISTRY_HPP_
#define CELLMUTATIONSTATEREGISTRY_HPP_

#include <boost/shared_ptr.hpp>
#include <vector>

#include "AbstractCellMutationState.hpp"

/**
 * A singleton registry of available mutation states.
 */
class CellMutationStateRegistry
{
public:
	/**
	 * The main interface to this class: get a particular mutation state object.
	 * Use like:
	 *    boost::shared_ptr<AbstractCellMutationState> p_state(
				CellMutationStateRegistry::Get<WildTypeCellMutationState>());
	 */
	template<class SUBCLASS>
	static boost::shared_ptr<AbstractCellMutationState> Get();

    /**
     * Get the single instance of the registry.
     */
    static CellMutationStateRegistry* Instance();

    /**
     * Get a list of the mutation states registered.
     */
    static const std::vector<boost::shared_ptr<AbstractCellMutationState> >& rGetAllMutationStates();

    /**
     * Clear all registered mutation states.
     */
    static void Clear();

private:

    /**
     * Default constructor.
     */
	CellMutationStateRegistry();

    /**
     * Copy constructor.
     */
	CellMutationStateRegistry(const CellMutationStateRegistry&);

    /**
     * Overloaded assignment operator.
     */
	CellMutationStateRegistry& operator= (const CellMutationStateRegistry&);

    /**
     * A pointer to the singleton instance of this class.
     */
    static CellMutationStateRegistry* mpInstance;

    /**
     * The mutation states in the registry.
     */
    std::vector<boost::shared_ptr<AbstractCellMutationState> > mMutationStates;
};

template<class SUBCLASS>
boost::shared_ptr<AbstractCellMutationState> CellMutationStateRegistry::Get()
{
	boost::shared_ptr<AbstractCellMutationState> p_state;
	std::vector<boost::shared_ptr<AbstractCellMutationState> >& r_states = Instance()->mMutationStates;
	for (unsigned i=0; i<r_states.size(); i++)
	{
		if (r_states[i]->IsType<SUBCLASS>())
		{
			p_state = r_states[i];
			break;
		}
	}
	if (!p_state)
	{
		// Create a new mutation state
		p_state.reset(new SUBCLASS);
		r_states.push_back(p_state);
	}
	return p_state;
}


#endif /* CELLMUTATIONSTATEREGISTRY_HPP_ */
