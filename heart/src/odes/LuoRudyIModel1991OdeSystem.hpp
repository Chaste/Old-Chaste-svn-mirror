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
#ifndef _LUORUDYIMODEL1991ODESYSTEM_HPP_
#define _LUORUDYIMODEL1991ODESYSTEM_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>

// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 * This class sets up the LuoRudyIModel1991OdeSystem system of equations.
 */
class LuoRudyIModel1991OdeSystem : public AbstractCardiacCell
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
        archive & boost::serialization::base_object<AbstractCardiacCell>(*this);
    }

    /** Constants for the LuoRudyIModel1991OdeSystem model */
    
    /** membrane capcaitance*/
    static const double membrane_C;
    /** Faraday constant*/
    static const double membrane_F;
    /** Universal gas constant*/
    static const double membrane_R;
    /** Temeperature*/
    static const double membrane_T;
    /** Reversal potentila for background current*/
    static const double background_current_E_b;
    /** Maximal conductance for background current*/
    static const double background_current_g_b;
    /** Maximal conductance for sodium current*/
    static const double fast_sodium_current_g_Na;
    /** Intracellular potassium concentration*/
    static const double ionic_concentrations_Ki;
    /** Extracellular potassium concentration*/
    static const double ionic_concentrations_Ko;
    /** Intracellular sodium concentration*/
    static const double ionic_concentrations_Nai;
    /** Extracellular sodium concentration*/
    static const double ionic_concentrations_Nao;
    /** Maximal conductance for plateau current*/
    static const double plateau_potassium_current_g_Kp;
    /** Permeability ratio Na/K for potassium currents*/
    static const double time_dependent_potassium_current_PR_NaK;

    /** Another parameter, which is a function of the above */
    double fast_sodium_current_E_Na;

    /**
     *  Range-checking on the current values of the state variables. Make sure
     *  all gating variables have are within zero and one, and all concentrations
     *  are positive
     */
    void VerifyStateVariables();

public:
    /**
     * Constructor
     * 
     * @param pSolver is a pointer to the ODE solver
     * @param pIntracellularStimulus is a pointer to the intracellular stimulus
     */
    LuoRudyIModel1991OdeSystem(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                               boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Destructor
     */
    ~LuoRudyIModel1991OdeSystem();

    /**
     * Fill in a vector representing the RHS of the TenTusscher2006 system
     * of Odes at each time step, y' = [y1' ... yn'].
     * Some ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time  the current time, in milliseconds
     * @param rY  current values of the state variables
     * @param rDY  to be filled in with derivatives
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

    /**
     * Returns the ionic current
     * 
     * @return the total ionic current
     */
    double GetIIonic();
    
    /**
     * Get the intracellular calcium concentration
     * 
     * @return the intracellular calcium concentration
     */
    double GetIntracellularCalciumConcentration();
};

BOOST_CLASS_EXPORT(LuoRudyIModel1991OdeSystem);

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
    Archive & ar, const LuoRudyIModel1991OdeSystem * t, const unsigned int file_version)
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
    Archive & ar, LuoRudyIModel1991OdeSystem * t, const unsigned int file_version)
{

    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
    boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
    ar >> p_solver;
    ar >> p_stimulus;
    ::new(t)LuoRudyIModel1991OdeSystem(p_solver, p_stimulus);
}
}
} // namespace ...

#endif // _LUORUDYIMODEL1991ODESYSTEM_HPP_
