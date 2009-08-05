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

#ifndef FASTSLOWODESYSTEMINFORMATION_HPP_
#define FASTSLOWODESYSTEMINFORMATION_HPP_

#include <boost/shared_ptr.hpp>

#include "AbstractOdeSystem.hpp"
#include "AbstractOdeSystemInformation.hpp"

/**
 * A concrete implementation of AbstractOdeSystemInformation, that uses templates
 * to provide an implementation for any AbstractFastSlowCardiacCell subclass.
 *
 * All ODE system developers need to do is provide a specialisation of the
 * Initialise method of this class, and set mpSystemInfo in their constructor:
 *   mpSystemInfo = FastSlowOdeSystemInformation<CLASS>::Instance();
 *
 * This class contains all the machinery to make it a singleton, hence providing
 * exactly one instance per value of the template parameter.
 */
template <bool ALL_VARS, class ODE_SYSTEM>
class FastSlowOdeSystemInformation : public AbstractOdeSystemInformation
{
private:
    /**
     * The single instance of this class, for this ODE_SYSTEM.
     */
    static boost::shared_ptr<FastSlowOdeSystemInformation<ALL_VARS, ODE_SYSTEM> > mpInstance;

protected:
    /**
     * Default constructor.
     *
     * Not user accessible - to obtain an instance of this class use the Instance
     * method.
     */
    FastSlowOdeSystemInformation();
    /**
     * Copy constructor
     */
    FastSlowOdeSystemInformation(const FastSlowOdeSystemInformation<ALL_VARS, ODE_SYSTEM>&);
    /**
     * Assignment operator
     */
    FastSlowOdeSystemInformation& operator= (const FastSlowOdeSystemInformation<ALL_VARS, ODE_SYSTEM>&);

    /**
     * Generic implementation of Initialise, which does nothing.
     *
     * Developers should specialise this method to their ODE system.  For example,

template<bool ALL_VARS>
void FastSlowOdeSystemInformation<ALL_VARS, MyNewFastSlowCell>::Initialise(void)
{
    this->mVariableNames.push_back("Variable_1");
    this->mVariableUnits.push_back("Units_1");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}
     */
    void Initialise(void);

public:
    /**
     * Return a (smart) pointer to the singleton instance, creating it if necessary.
     */
    static boost::shared_ptr<FastSlowOdeSystemInformation<ALL_VARS, ODE_SYSTEM> > Instance();
};


template<bool ALL_VARS, class ODE_SYSTEM>
boost::shared_ptr<FastSlowOdeSystemInformation<ALL_VARS, ODE_SYSTEM> > FastSlowOdeSystemInformation<ALL_VARS, ODE_SYSTEM>::Instance()
{
    if (!mpInstance)
    {
        mpInstance.reset(new FastSlowOdeSystemInformation<ALL_VARS, ODE_SYSTEM>);
        mpInstance->Initialise();
    }
    return mpInstance;
}

template<bool ALL_VARS, class ODE_SYSTEM>
FastSlowOdeSystemInformation<ALL_VARS, ODE_SYSTEM>::FastSlowOdeSystemInformation()
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == NULL);
}

template<bool ALL_VARS, class ODE_SYSTEM>
void FastSlowOdeSystemInformation<ALL_VARS, ODE_SYSTEM>::Initialise(void)
{
    // does nothing; designed to be specialised
}

/**
 * Definition of the instance static member.
 */
template<bool ALL_VARS, class ODE_SYSTEM>
boost::shared_ptr<FastSlowOdeSystemInformation<ALL_VARS, ODE_SYSTEM> > FastSlowOdeSystemInformation<ALL_VARS, ODE_SYSTEM>::mpInstance;

#endif /*FASTSLOWODESYSTEMINFORMATION_HPP_*/
