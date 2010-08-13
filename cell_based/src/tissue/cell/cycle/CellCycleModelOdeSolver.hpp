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

#ifndef CELLCYCLEMODELODESOLVER_HPP_
#define CELLCYCLEMODELODESOLVER_HPP_

#include "ChasteSerialization.hpp"

#include "AbstractCellCycleModelOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"

/**
 * A concrete implementation of AbstractCellCycleModelOdeSolver, that uses templates
 * to provide an implementation for any cell cycle model ODE solver class.
 *
 * All ODE-based cell cycle model developers need to do is provide a specialisation of
 * the Initialise method of this class, and set mpOdeSolver in their constructor:
 *   mpOdeSolver = CellCycleModelOdeSolver<CLASS>::Instance();
 *
 * This class contains all the machinery to make it a singleton, hence providing
 * exactly one instance per value of the template parameter.
 */
template <class CELL_CYCLE_MODEL, class ODE_SOLVER>
class CellCycleModelOdeSolver : public AbstractCellCycleModelOdeSolver
{
private:
    /** The single instance of this class, for this ODE_SOLVER. */
    static boost::shared_ptr<CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER> > mpInstance;

    /** Default constructor. Not user accessible; to obtain an instance of this class use the Instance method. */
    CellCycleModelOdeSolver();

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell cycle model, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModelOdeSolver>(*this);
        archive & mpInstance;
    }

public:

    /** Copy constructor. */
    CellCycleModelOdeSolver(const CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>&);

    /**
     * Overloaded assignment operator.
     * 
     * @param rOtherCellCycleModelOdeSolver another CellCycleModelOdeSolver
     */
    CellCycleModelOdeSolver& operator= (const CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>& rOtherCellCycleModelOdeSolver);

    /** Return a pointer to the singleton instance, creating it if necessary. */
    static boost::shared_ptr<CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER> > Instance();

    /** Is the instance in existence and fully set up. */
    bool IsSetUp();

    /** Initialise the ODE solver. */
    void Initialise();
};

/** Definition of the instance static member. */
template<class CELL_CYCLE_MODEL, class ODE_SOLVER>
boost::shared_ptr<CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER> > CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>::mpInstance;

template<class CELL_CYCLE_MODEL, class ODE_SOLVER>
CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>::CellCycleModelOdeSolver()
    : AbstractCellCycleModelOdeSolver()
{
    // Make sure there's only one instance; enforces correct serialization
    ///\todo #1427 I don't think this is doing what you think it is - on executing this,
    /// any objects pointing to the old instance will still point to it, and because they're
    /// using shared_ptr, will keep it alive.
    if (mpInstance)
    {
        mpInstance.reset();
    }
}

template<class CELL_CYCLE_MODEL, class ODE_SOLVER>
CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>& CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>::operator= (const CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>& rOtherCellCycleModelOdeSolver)
{
    mpInstance = rOtherCellCycleModelOdeSolver.mpInstance;
    mpOdeSolver = rOtherCellCycleModelOdeSolver.mpOdeSolver;
    mSizeOfOdeSystem = rOtherCellCycleModelOdeSolver.mSizeOfOdeSystem;
}

template<class CELL_CYCLE_MODEL, class ODE_SOLVER>
boost::shared_ptr<CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER> > CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>::Instance()
{
    if (!mpInstance)
    {
        mpInstance.reset(new CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>);
    }
    return mpInstance;
}

template<class CELL_CYCLE_MODEL, class ODE_SOLVER>
bool CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>::IsSetUp()
{
    return mpInstance && mpOdeSolver;
}

template<class CELL_CYCLE_MODEL, class ODE_SOLVER>
void CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>::Initialise()
{
    mpOdeSolver = boost::shared_ptr<AbstractIvpOdeSolver>(new ODE_SOLVER);
}


/**
 * Specialization for BackwardEulerIvpOdeSolver, whose constructor requires
 * an argument.
 * \todo there must be an easier way to deal with this peculiarity (#1427)
 */
template<class CELL_CYCLE_MODEL>
class CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver> : public AbstractCellCycleModelOdeSolver
{
private:

    /** The single instance of this class, for this ODE_SOLVER. */
    static boost::shared_ptr<CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver> > mpInstance;

    /** Default constructor. Not user accessible; to obtain an instance of this class use the Instance method. */
    CellCycleModelOdeSolver();

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell cycle model, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModelOdeSolver>(*this);
        archive & mpInstance;
    }

public:

    /** Copy constructor. */
    CellCycleModelOdeSolver(const CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver>&);

    /** Overloaded assignment operator. */
    CellCycleModelOdeSolver& operator= (const CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver>&);

    /** Return a pointer to the singleton instance, creating it if necessary. */
    static boost::shared_ptr<CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver> > Instance();

    /** Is the instance in existence and fully set up. */
    bool IsSetUp();

    /** Initialise the ODE solver. */
    void Initialise();

    /** Reset the instance. */
    void Reset();
};

template<class CELL_CYCLE_MODEL>
boost::shared_ptr<CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver> > CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver>::mpInstance;

template<class CELL_CYCLE_MODEL>
CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver>::CellCycleModelOdeSolver()
    : AbstractCellCycleModelOdeSolver()
{
    // Make sure there's only one instance - enforces correct serialization
    if (mpInstance)
    {
        mpInstance.reset();
    }
}

template<class CELL_CYCLE_MODEL>
boost::shared_ptr<CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver> > CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver>::Instance()
{
    if (!mpInstance)
    {
        mpInstance.reset(new CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver>);
    }
    return mpInstance;
}

template<class CELL_CYCLE_MODEL>
bool CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver>::IsSetUp()
{
    return (mpInstance!=NULL) && (mpOdeSolver!=NULL) && (mSizeOfOdeSystem != UNSIGNED_UNSET);
}

template<class CELL_CYCLE_MODEL>
void CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver>::Initialise()
{
    if (mSizeOfOdeSystem == UNSIGNED_UNSET)
    {
        EXCEPTION("SetSizeOfOdeSystem() must be called before calling Initialise()");
    }
    mpOdeSolver = boost::shared_ptr<AbstractIvpOdeSolver>(new BackwardEulerIvpOdeSolver(mSizeOfOdeSystem));
}

template<class CELL_CYCLE_MODEL>
void CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver>::Reset()
{
    mSizeOfOdeSystem = UNSIGNED_UNSET;
    mpOdeSolver.reset();
}

#endif /*CELLCYCLEMODELODESOLVER_HPP_*/
