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
#ifndef ABSTRACTWNTODEBASEDCELLCYCLEMODEL_HPP_
#define ABSTRACTWNTODEBASEDCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractOdeBasedCellCycleModel.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "CvodeAdaptor.hpp"
#include "WntConcentration.hpp"

/**
 * This class contains all the things common to the Wnt cell cycle ODE based models,
 * the resetting method and updating of cell types etc.
 *
 * The concrete models all need to operate with a WntConcentration
 * singleton object.
 *
 * These models have a constant length M phase, run ODEs to decide when
 * to finish G1 phase, then add time for S and G2 phases (in some classes,
 * random periods of time). The CellProliferativeType is updated dependent on the
 * concentration of beta-catenin (given by one of the ODEs).
 */
class AbstractWntOdeBasedCellCycleModel : public AbstractOdeBasedCellCycleModel
{
private:

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
        archive & boost::serialization::base_object<AbstractOdeBasedCellCycleModel>(*this);
        archive & mDimension;
    }

protected:

#ifdef CHASTE_CVODE
    /** The ODE solver. */
    static CvodeAdaptor msSolver;
#else
    /** The ODE solver. */
    static RungeKutta4IvpOdeSolver msSolver;
#endif //CHASTE_CVODE

    /**
     * The spatial dimension (needed by the templated class WntConcentration).
     */
    unsigned mDimension;

    /**
     * @return time when the ODEs reached their stopping condition.
     */
    virtual double GetOdeStopTime();

    /**
     * Get the Wnt level experienced by the cell.
     */
    double GetWntLevel();

public:

    /**
     * Default constructor
     *
     * @param dimension the spatial dimension (used by WntConcentration)
     */
    AbstractWntOdeBasedCellCycleModel();

    /**
     * Copy constructor.
     *
     * This is needed because we store and manage a pointer to an ODE system.
     * Note that this class doesn't actually copy the ODE system, because each
     * subclass will use a different type.  Hence subclasses *must* copy their
     * own ODE system in their copy constructor.
     *
     * @param rOtherModel the cell cycle model being copied.
     */
    AbstractWntOdeBasedCellCycleModel(const AbstractWntOdeBasedCellCycleModel& rOtherModel);

    /**
     * Resets the Wnt Model to the start of the cell cycle (this model does not cycle naturally)
     * Cells are given a new birth time and cell cycle proteins are reset.
     * Note that the Wnt pathway proteins maintain their current values.
     *
     * Should only be called by the TissueCell::Divide() method.
     */
    void ResetForDivision();

    /**
     * Updates the current cell type to reflect whether the
     * beta-catenin level has dropped low enough to make it stop dividing.
     * This should only be called when the cell cycle model has been
     * evaluated to the current time, or it may give misleading results.
     */
    void UpdateCellProliferativeType();

    /**
     * This must be implemented by subclasses to change cell type to reflect
     * current levels of beta-catenin.
     */
    virtual void ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel()=0;

    /**
	 * Set the spatial dimension.
	 *
	 * @param dimension
	 */
	void SetDimension(unsigned dimension);

	/**
     * Get the spatial dimension.
     *
     * @return mDimension
     */
    unsigned GetDimension();

};

CLASS_IS_ABSTRACT(AbstractWntOdeBasedCellCycleModel)

#endif /*ABSTRACTWNTODEBASEDCELLCYCLEMODEL_HPP_*/
