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
#ifndef CELLPROPERTYCOLLECTION_HPP_
#define CELLPROPERTYCOLLECTION_HPP_

#include <list>
#include <boost/shared_ptr.hpp>

#include "AbstractCellProperty.hpp"

class CellPropertyCollection
{
private:
	/** The type of container used to store properties */
	typedef std::list<boost::shared_ptr<AbstractCellProperty> > CollectionType;

	/** Type of a const iterator over the container */
	typedef CollectionType::const_iterator ConstIteratorType;

	/** The properties stored in this collection. */
	CollectionType mProperties;

public:
	/**
	 * Create an empty collection of cell properties.
	 */
	CellPropertyCollection();

	/**
	 * Add a new property to this collection.
	 *
	 * @param rProp  the property to add
	 */
	void AddProperty(const boost::shared_ptr<AbstractCellProperty>& rProp);

	/**
	 * Test whether this collection contains the given property @b object.
	 *
	 * @param rProp  the property to compare against
	 */
	bool HasProperty(const boost::shared_ptr<AbstractCellProperty>& rProp) const;

	/**
	 * Test whether the collection contains a property that has the exact type CLASS.
	 *
	 * Should be used like
	 *   bool healthy = collection.HasProperty<WildTypeCellMutationState>();
	 */
	template <typename CLASS>
	bool HasProperty() const
	{
		for (ConstIteratorType it = mProperties.begin(); it != mProperties.end(); ++it)
		{
			if ((*it)->IsType<CLASS>())
			{
				return true;
			}
		}
		return false;
	}

	/**
	 * Test whether the collection contains a property that inherits from BASECLASS.
	 *
	 * Should be used like
	 *   collection.HasPropertyType<AbstractCellMutationState>();
	 */
	template <typename BASECLASS>
	bool HasPropertyType() const
	{
		for (ConstIteratorType it = mProperties.begin(); it != mProperties.end(); ++it)
		{
			if ((*it)->IsSubType<BASECLASS>())
			{
				return true;
			}
		}
		return false;
	}
};

#endif /* CELLPROPERTYCOLLECTION_HPP_ */
