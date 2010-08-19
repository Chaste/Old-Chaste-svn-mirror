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
#ifndef SINGLEODEWNTCELLCYCLEMODEL_HPP_
#define SINGLEODEWNTCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <cfloat>

#include "SimpleWntCellCycleModel.hpp"
#include "Mirams2010WntOdeSystem.hpp"
#include "AbstractCellCycleModelOdeSolver.hpp"

/**
 * Wnt-dependent cell cycle model. Needs to operate with a WntConcentration
 * singleton object.
 *
 * This model has a constant length M phase, runs ODEs to decide when
 * to finish G1 phase then adds time for S and G2 phases. The CellProliferativeType is
 * updated dependent on the concentration of beta-catenin (given by one
 * of the ODEs).
 */
class SingleOdeWntCellCycleModel : public SimpleWntCellCycleModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell cycle model and member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<SimpleWntCellCycleModel>(*this);
        archive & mpOdeSystem;
        archive & mpOdeSolver;
        archive & mBetaCateninDivisionThreshold;
        archive & mLastTime;
    }

    /**
     * Pointer to the ODE system developed by Mirams et al (2010).
     */
    AbstractOdeSystem* mpOdeSystem;

    /**
     * The ODE solver.
     * 
     * Subclasses need to set this in their constructor to point to an instance
     * of a suitable class. See for example the CellCycleModelOdeSolver class.
     */
    boost::shared_ptr<AbstractCellCycleModelOdeSolver> mpOdeSolver;

    /**
     * The cell differentiates when the beta-catenin level drops
     * below this value. It is hard coded in
     * Initialise() because there are so many constructors.
     *
     * Set and Get methods are also provided.
     */
    double mBetaCateninDivisionThreshold;

    /** The last time at which the ODEs were solved up to */
    double mLastTime;

    /**
     * Called by ::Initialise() and ::UpdateCellProliferativeType() only.
     * Updates the mpCell::mCellProliferativeType to match mpOdeSystem's
     * beta-catenin levels
     *
     * This carries out the work for ::UpdateCellProliferativeType();
     * But does not check the current time so it can be used by the initialise method.
     */
    void ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();

    /**
     * This method runs the ODEs and updates the beta-catenin level.
     */
    void UpdateBetaCateninLevel();

public:

    /**
     * Default constructor.
     * 
     * @param pOdeSolver An optional pointer to a cell cycle model ODE solver object (allows the use of different ODE solvers)
     */
    SingleOdeWntCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Constructor used in archiving.
     * 
     * @param unused an unused argument
     */
    SingleOdeWntCellCycleModel(double unused)
    {}

    /**
     * Destructor.
     */
    ~SingleOdeWntCellCycleModel();

    /**
     * Initialise the cell cycle model at the start of a simulation.
     *
     * This overridden method sets up a new WntCellCycleOdeSystem,
     * sets the cell type according to the current beta catenin level
     * and sets a random G2 duration.
     */
    void Initialise();

    /**
     * This specialisation updates the beta-catenin level
     */
    void UpdateCellCyclePhase();

    /**
     * Overridden builder method to create new copies of
     * this cell cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Return the total beta-catenin concentration
     */
    double GetBetaCateninConcentration();

    /**
     * Set #mBetaCateninDivisionThreshold.
     *
     * @param betaCateninDivisionThreshold to be set
     */
    void SetBetaCateninDivisionThreshold(double betaCateninDivisionThreshold);

    /**
     * Get #mBetaCateninDivisionThreshold.
     */
    double GetBetaCateninDivisionThreshold();

     /**
      * @return mpOdeSolver (used in archiving).
      */
    const boost::shared_ptr<AbstractCellCycleModelOdeSolver> GetOdeSolver() const;

    /**
     * Set mLastTime.
     * 
     * @param lastTime the new value of mLastTime
     */
    void SetLastTime(double lastTime);

    /**
     * Set the values of the state variables in the cell cycle model's ODE system.
     *
     * @param rStateVariables vector containing values for the state variables
     */
    void SetStateVariables(const std::vector<double>& rStateVariables);

    /**
     * Set mpOdeSystem. Used in CreateCellCycleModel().
     * 
     * @param pOdeSystem the ODE system
     */
    void SetOdeSystem(AbstractOdeSystem* pOdeSystem);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(SingleOdeWntCellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(SingleOdeWntCellCycleModel)

namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a SingleOdeWntCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const SingleOdeWntCellCycleModel * t, const unsigned int file_version)
{
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a SingleOdeWntCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, SingleOdeWntCellCycleModel * t, const unsigned int file_version)
{
    double unused = 0.0;
    ::new(t)SingleOdeWntCellCycleModel(unused);
}
}
} // namespace


#endif /*SINGLEODEWNTCELLCYCLEMODEL_HPP_*/
