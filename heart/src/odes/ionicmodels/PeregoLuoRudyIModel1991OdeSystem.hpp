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
#ifndef _PEREGOLUORUDYIMODEL1991ODESYSTEM_HPP_
#define _PEREGOLUORUDYIMODEL1991ODESYSTEM_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractCardiacCell.hpp"
#include "AbstractPeregoCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * This class sets up the PeregoLuoRudyIModel1991OdeSystem system of equations.
 */
class PeregoLuoRudyIModel1991OdeSystem : public AbstractPeregoCardiacCell
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

    /** Constants for the PeregoLuoRudyIModel1991OdeSystem model */
    
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
    
    /**
     * Helper function to compute the state variables according to the Perego scheme.
     * It operates on the state variables of the parent class ma_current and mb_current
     * 
     * @param rY a vector of state variables
     * @param currentTime the current time
     */
    void ComputeSystemParameters(const std::vector<double>& rY, double currentTime);
    
    
public:
    /**
     * Constructor
     * 
     * @param pIntracellularStimulus is a pointer to the intracellular stimulus
     */
    PeregoLuoRudyIModel1991OdeSystem(boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Destructor
     */
    ~PeregoLuoRudyIModel1991OdeSystem();

    /**
     * Returns the ionic current
     * 
     * @return the total ionic current
     */
    double GetIIonic();

};
//
//#include "TemplatedExport.hpp"
//CHASTE_CLASS_EXPORT(PeregoLuoRudyIModel1991OdeSystem);
//
//namespace boost
//{
//namespace serialization
//{
///**
// * Allow us to not need a default constructor, by specifying how Boost should
// * instantiate a PeregoLuoRudyIModel1991OdeSystem instance.
// */
//template<class Archive>
//inline void save_construct_data(
//    Archive & ar, const PeregoLuoRudyIModel1991OdeSystem * t, const unsigned int file_version)
//{
//    const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();
//    const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();
//    ar << p_solver;
//    ar << p_stimulus;
//}
//
///**
// * Allow us to not need a default constructor, by specifying how Boost should
// * instantiate a PeregoLuoRudyIModel1991OdeSystem instance (using existing constructor)
// *
// * NB this constructor allocates memory for the other member variables too.
// */
//template<class Archive>
//inline void load_construct_data(
//    Archive & ar, PeregoLuoRudyIModel1991OdeSystem * t, const unsigned int file_version)
//{
//
//    boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
//    ar >> p_stimulus;
//    ::new(t)PeregoLuoRudyIModel1991OdeSystem(p_stimulus);
//}
//}
//} // namespace ...

#endif // _PEREGOLUORUDYIMODEL1991ODESYSTEM_HPP_
