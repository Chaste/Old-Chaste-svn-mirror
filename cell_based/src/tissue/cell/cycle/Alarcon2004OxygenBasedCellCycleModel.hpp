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
#ifndef ALARCON2004OXYGENBASEDCELLCYCLEMODEL_HPP_
#define ALARCON2004OXYGENBASEDCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include <cfloat>

#include "AbstractOdeBasedCellCycleModelWithStoppingEvent.hpp"
#include "Alarcon2004OxygenBasedCellCycleOdeSystem.hpp"
#include "CellwiseData.hpp"
#include "CellLabel.hpp"
#include "Exception.hpp"

#include "CellCycleModelOdeSolver.hpp"
#ifdef CHASTE_CVODE
#include "CvodeAdaptor.hpp"
#endif //CHASTE_CVODE
#include "BackwardEulerIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "HeunIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"

/**
 * Oxygen-dependent ODE-based cell cycle model. Published by Alarcon et al.
 * (doi:10.1016/j.jtbi.2004.04.016).
 */
class Alarcon2004OxygenBasedCellCycleModel : public AbstractOdeBasedCellCycleModelWithStoppingEvent
{
private:

    /** Whether the cell associated with this cell cycle model is labelled (this affects the ODE system). */
    bool mIsLabelled;

    ///\todo Archiving could be tidied up for this class
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell cycle model and ODE system.
     *
     * @param archive the archive
     * @param version the archive version
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        assert(mpOdeSystem);
        archive & boost::serialization::base_object<AbstractOdeBasedCellCycleModelWithStoppingEvent>(*this);
        bool is_labelled = static_cast<Alarcon2004OxygenBasedCellCycleOdeSystem*>(mpOdeSystem)->IsLabelled();
        archive & is_labelled;
    }
    /**
     * Load the cell cycle model and ODE system from archive.
     *
     * @param archive the archive
     * @param version the archive version
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        // The ODE system is set up by the archiving constructor, so we can set the mutation state
        // here.  This is a horrible hack, but avoids having to regenerate test archives...
        assert(mpOdeSystem);
        archive & boost::serialization::base_object<AbstractOdeBasedCellCycleModelWithStoppingEvent>(*this);
        bool is_labelled;
        archive & is_labelled;
        static_cast<Alarcon2004OxygenBasedCellCycleOdeSystem*>(mpOdeSystem)->SetIsLabelled(is_labelled);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

public:

    /**
     * Default constructor.
     * 
     * @param pOdeSolver An optional pointer to a cell cycle model ODE solver object (allows the use of different ODE solvers)
     */
    Alarcon2004OxygenBasedCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Copy constructor.
     *
     * Also copies our ODE system.
     *
     * @param rOtherModel the instance being copied.
     */
    Alarcon2004OxygenBasedCellCycleModel(const Alarcon2004OxygenBasedCellCycleModel& rOtherModel);

    /**
     * A private constructor for archiving.
     * 
      *@param pOdeSolver a pointer to a cell cycle model ODE solver object (allows the use of different ODE solvers)
     * @param rParentProteinConcentrations a std::vector of doubles of the protein concentrations (see WntCellCycleOdeSystem)
     * @param rDimension the spatial dimension
     * @param isLabelled whether the cell associated with this cell cycle model is labelled (this affects the ODE system)
     */
    Alarcon2004OxygenBasedCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver,
                                         const std::vector<double>& rParentProteinConcentrations,
                                         const unsigned& rDimension,
                                         bool isLabelled);

    /**
     * Resets the oxygen-based model to the start of the cell cycle
     * (this model does not cycle naturally). Cells are given a new
     * birth time and cell cycle proteins are reset. Note that the
     * oxygen concentration maintains its current value.
     *
     * Should only be called by the TissueCell Divide() method.
     */
    virtual void ResetForDivision();

    /**
     * Overridden builder method to create new copies of
     * this cell cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Initialise the cell cycle model at the start of a simulation.
     *
     * This overridden method sets up a new ODE system.
     */
    void Initialise();

    /**
     * Solve the ODEs up to the current time and return whether a stopping event occurred.
     *
     * @param currentTime the current time
     * @return whether a stopping event occurred
     */
    bool SolveOdeToTime(double currentTime);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Alarcon2004OxygenBasedCellCycleModel)


namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a Alarcon2004OxygenBasedCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const Alarcon2004OxygenBasedCellCycleModel * t, const unsigned int file_version)
{
    const boost::shared_ptr<AbstractCellCycleModelOdeSolver> p_ode_solver = t->GetOdeSolver();
    ar & p_ode_solver;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a Alarcon2004OxygenBasedCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, Alarcon2004OxygenBasedCellCycleModel * t, const unsigned int file_version)
{
    boost::shared_ptr<AbstractCellCycleModelOdeSolver> p_ode_solver;
    ar & p_ode_solver;

    /**
     * Invoke inplace constructor to initialise an instance of Alarcon2004OxygenBasedCellCycleModel.
     * It doesn't actually matter what values we pass to our standard constructor, provided they are
     * valid parameter values, since the state loaded later from the archive will overwrite their
     * effect in this case.
     */

    std::vector<double> state_vars;
    for (unsigned i=0; i<6; i++)
    {
        state_vars.push_back(0.0);
    }
    unsigned dimension = 1;
    bool is_labelled = false;

    ::new(t)Alarcon2004OxygenBasedCellCycleModel(p_ode_solver, state_vars, dimension, is_labelled);
}
}
} // namespace ...


#ifdef CHASTE_CVODE
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, Alarcon2004OxygenBasedCellCycleModel, CvodeAdaptor)
#endif //CHASTE_CVODE
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, Alarcon2004OxygenBasedCellCycleModel, BackwardEulerIvpOdeSolver)
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, Alarcon2004OxygenBasedCellCycleModel, EulerIvpOdeSolver)
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, Alarcon2004OxygenBasedCellCycleModel, HeunIvpOdeSolver)
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, Alarcon2004OxygenBasedCellCycleModel, RungeKutta2IvpOdeSolver)
EXPORT_TEMPLATE_CLASS2(CellCycleModelOdeSolver, Alarcon2004OxygenBasedCellCycleModel, RungeKutta4IvpOdeSolver)

#endif /*ALARCON2004OXYGENBASEDCELLCYCLEMODEL_HPP_*/
