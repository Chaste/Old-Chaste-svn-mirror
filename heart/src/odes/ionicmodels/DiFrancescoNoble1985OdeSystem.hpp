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
#ifndef _DIFRANCESCONOBLE1985ODESYSTEM_HPP_
#define _DIFRANCESCONOBLE1985ODESYSTEM_HPP_

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * This class sets up the DifrancescoNoble1985 equations.
 * In this version of the model the time is in milliseconds (unlike the paper).
 */
class DiFrancescoNoble1985OdeSystem : public AbstractCardiacCell
{
private:

    /* Constants for the DifrancescoNoble1985 model */
    static const double membrane_C; /**< nanoF. It must be so because dV/dt is in mV/ms and currents are in nA. */
    static const double radius; /**< mm */
    static const double length; /**< mmm */
    static const double V_e_ratio; /**< dimensionless */
    static const double R; /**< J/(kmol*Kelvin) */
    static const double T; /**< Kelvin */
    static const double F; /**< C/mol */

    // conductances and the like
    static const double g_fna; /**< microS */
    static const double g_fk; /**< microS */
    static const double g_na_b; /**< microS */
    static const double g_ca_b; /**< microS */
    static const double g_na; /**< microS */
    static const double g_k1; /**< microS */
    static const double g_to; /**< microS/mM */
    static const double I_P; /**< nA */
    static const double i_kmax; /**< nA */
    static const double k_naca; /**< nA */
    static const double P_si; /**< nA/mM */
    //concentrations
    static const double Nao; /**< mM */
    static const double Cao; /**< mM */
    static const double Kb; /**< mM */
    static const double Kmf; /**< mM */
    static const double Km1; /**< mM */
    static const double Kmto; /**< mM */
    static const double KmK; /**< mM */
    static const double KmNa; /**< mM */
    static const double Kmf2; /**< mM */
    static const double KmCa; /**< mM */
    static const double Km_Ca; /**< mM */
    //others
    static const double n_naca; /**< dimensionless */
    static const double gamma; /**< dimensionless */
    static const double d_naca; /**< dimensionless */
    static const double rCa; /**< dimensionless */
    static const double tau_up; /**< seconds. Must be seconds to cancel with numerator of F and give nA for i_up */
    static const double tau_rep; /**< seconds. Must be seconds to cancel with numerator of F and give nA for i_rep */
    static const double tau_rel; /**< seconds. Must be seconds to cancel with numerator of F and give nA for i_rel */
    static const double Ca_up_max; /**< mM */
    static const double pf; /**< per millisecond */
    static const double Vecs; /**< dimensionless */

    /** total intracellular volume*/
    double Vi;
    /** uptake volume (of the SR)*/
    double Vup;
    /** release volume (of the SR)*/
    double Vrel;
    /** extracellular volume*/
    double Ve;
    /** cell volume*/
    double Vcell;
    /** R times T over F*/
    double RToNF;
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
    DiFrancescoNoble1985OdeSystem(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                                  boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Destructor
     */
    ~DiFrancescoNoble1985OdeSystem();

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

};

#endif // _DIFRANCESCONOBLE_HPP_
