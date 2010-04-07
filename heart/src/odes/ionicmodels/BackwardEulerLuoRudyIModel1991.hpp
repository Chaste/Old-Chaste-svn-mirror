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
#ifndef _BACKWARDEULERLUORUDYIMODEL1991_HPP_
#define _BACKWARDEULERLUORUDYIMODEL1991_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractStimulusFunction.hpp"
#include "AbstractBackwardEulerCardiacCell.hpp"

#include <vector>

/**
 * This class sets up the Luo-Rudy I 1991 system of equations, and solves them
 * using a decoupled backward Euler approach.
 */
class BackwardEulerLuoRudyIModel1991 : public AbstractBackwardEulerCardiacCell<1>
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This calls serialize on the base class.
        archive & boost::serialization::base_object<AbstractBackwardEulerCardiacCell<1> >(*this);
    }
    /*
     * Constants for the LuoRudyIModel1991OdeSystem model
     */
    static const double membrane_C; /**< Membrane capacitance, uF/cm^2 */
    static const double membrane_F; /**< Faraday's constant, C/mol */
    static const double membrane_R; /**< Universal gas constant, J/kmol K */
    static const double membrane_T; /**< Temperature, K */
    static const double background_current_E_b; /**< Reversal potential for background current, mV */
    static const double background_current_g_b; /**< Maximal conductance for background current, mS/cm^2 */
    static const double fast_sodium_current_g_Na; /**< Maximal conductance for sodium current, mS/cm^2 */
    static const double ionic_concentrations_Ki; /**< Intracellular potassium concentration, mM */
    static const double ionic_concentrations_Ko; /**< Extracellular potassium concentration, mM */
    static const double ionic_concentrations_Nai; /**< Intracellular sodium concentration, mM */
    static const double ionic_concentrations_Nao; /**< Extracellular sodium concentration, mM */
    static const double plateau_potassium_current_g_Kp; /**< Maximal conductance for plateau potassium current, mS/cm^2 */
    static const double time_dependent_potassium_current_PR_NaK; /**< Permeability ratio of sodium over potassium, dimensionless */

    /** another parameter, which is a function of the above */
    double fast_sodium_current_E_Na;


public:
    /**
     * Constructor
     *
     * @param pIntracellularStimulus a pointer to the intracellular stimulus
     */
    BackwardEulerLuoRudyIModel1991(boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Constructor with the same signature as the forward cell models
     *
     * @param pIntracellularStimulus a pointer to the intracellular stimulus
     */
    BackwardEulerLuoRudyIModel1991(boost::shared_ptr<AbstractIvpOdeSolver> /* unused */,
                                   boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     *  Destructor
     */
    ~BackwardEulerLuoRudyIModel1991();

    /**
     *  Calculates the remaining parameters and calls the base class Init method.
     */
    void Init();

protected:
    /**
     * Compute the values of all state variables except the voltage for one
     * timestep.
     *
     * @param tStart  start of the timestep
     */
    void ComputeOneStepExceptVoltage(double tStart);

    /**
     * Perform a forward Euler step to update the transmembrane potential.
     *
     * @param time start of the timestep
     */
    void UpdateTransmembranePotential(double time);

public:
    /**
     * Compute the residual for the Newton iteration for the non-linear system portion of the model.
     *
     * @param rCurrentGuess  current values of the non-linear system variables
     * @param rResidual  to be filled in with the residual vector
     */
    void ComputeResidual(double var_environment__time, const double rCurrentGuess[1], double rResidual[1]);

    /**
     * Compute the Jacobian for the Newton iteration for the non-linear system portion of the model.
     *
     * @param rCurrentGuess  current values of the non-linear system variables
     * @param rJacobian  to be filled in with the jacobian matrix
     */
    void ComputeJacobian(double var_environment__time, const double rCurrentGuess[1], double rJacobian[1][1]);

    /**
     * Compute the ionic current at the current instant in time
     * (i.e. using the current values of the state variables).
     */
    double GetIIonic();

    /**
     *  Check that none of the gating variables have gone out of range. Throws an
     *  Exception if any have.
     */
    void VerifyStateVariables();

    /**
     * @returns the intracellular calcium concentration
     */
    double GetIntracellularCalciumConcentration();
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(BackwardEulerLuoRudyIModel1991)

namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a LuoRudyIModel1991OdeSystem instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const BackwardEulerLuoRudyIModel1991 * t, const unsigned int file_version)
{
    const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();
    const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();
    ar << p_solver;
    ar << p_stimulus;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a LuoRudyIModel1991OdeSystem instance (using existing constructor)
 *
 * NB this constructor allocates memory for the other member variables too.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, BackwardEulerLuoRudyIModel1991 * t, const unsigned int file_version)
{
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
    boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
    ar >> p_solver;
    ar >> p_stimulus;
    ::new(t)BackwardEulerLuoRudyIModel1991(p_solver, p_stimulus);
}
}
} // namespace ...

#endif // _BACKWARDEULERLUORUDYIMODEL1991_HPP_
