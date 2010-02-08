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

#ifndef _TENTUSSCHER2006ODESYSTEM_HPP_
#define _TENTUSSCHER2006ODESYSTEM_HPP_

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * This class sets up the equations for the TenTusscher 2006 model. This is for an epicardial cell.
 */
class TenTusscher2006OdeSystem : public AbstractCardiacCell
{
private:

    
    double mScaleFactorGks;/**< Scale factor for Gks*/
    double mScaleFactorIto;/**< Scale factor for Gto*/
    double mScaleFactorGkr;/**< Scale factor for Gkr*/

    //////////////////////////////////////////////////////////////
    //Constants for the TenTusscher2006 model, values for epicardial cell.
    //////////////////////////////////////////////////////////////
    static const double L_type_Ca_current_g_CaL; /**< nanoS_per_picoF */
    static const double calcium_background_current_g_bca; /**< nanoS_per_picoF */
    static const double calcium_dynamics_Buf_c;   /**< millimolar */
    static const double calcium_dynamics_Buf_sr;   /**< millimolar */
    static const double calcium_dynamics_Buf_ss;   /**< millimolar */
    static const double calcium_dynamics_Ca_o;   /**< millimolar */
    static const double calcium_dynamics_EC;   /**< millimolar */
    static const double calcium_dynamics_K_buf_c;   /**< millimolar */
    static const double calcium_dynamics_K_buf_sr;   /**< millimolar */
    static const double calcium_dynamics_K_buf_ss;   /**< millimolar */
    static const double calcium_dynamics_K_up;   /**< millimolar */
    static const double calcium_dynamics_V_leak;   /**< millimolar_per_millisecond */
    static const double calcium_dynamics_V_rel;   /**< millimolar_per_millisecond */
    static const double calcium_dynamics_V_sr;   /**< micrometre3 */
    static const double calcium_dynamics_V_ss;   /**< micrometre3 */
    static const double calcium_dynamics_V_xfer;   /**< millimolar_per_millisecond */
    static const double calcium_dynamics_Vmax_up;   /**< millimolar_per_millisecond */
    static const double calcium_dynamics_k1_prime;   /**< per_millimolar2_per_millisecond */
    static const double calcium_dynamics_k2_prime;   /**< per_millimolar_per_millisecond */
    static const double calcium_dynamics_k3;   /**< per_millisecond */
    static const double calcium_dynamics_k4;   /**< per_millisecond */
    static const double calcium_dynamics_max_sr;   /**< dimensionless */
    static const double calcium_dynamics_min_sr;   /**< dimensionless */
    static const double calcium_pump_current_K_pCa;   /**< millimolar */
    static const double calcium_pump_current_g_pCa;   /**< nanoS_per_picoF */
    static const double fast_sodium_current_g_Na;   /**< nanoS_per_picoF */
    static const double inward_rectifier_potassium_current_g_K1;   /**< nanoS_per_picoF */
    static const double membrane_Cm;   /**< microF_per_cm2 */
    static const double membrane_F;   /**< coulomb_per_millimole */
    static const double membrane_R;   /**< joule_per_mole_kelvin */
    static const double membrane_T;   /**< kelvin */
    static const double membrane_V_c;   /**< micrometre3 */
    static const double potassium_dynamics_K_o;   /**< millimolar */
    static const double potassium_pump_current_g_pK;   /**< nanoS_per_picoF */
    static const double rapid_time_dependent_potassium_current_g_Kr;   /**< nanoS_per_picoF */
    static const double reversal_potentials_P_kna;   /**< nanoA_per_millimolar */
    static const double slow_time_dependent_potassium_current_g_Ks;   /**< nanoS_per_picoF */
    static const double sodium_background_current_g_bna;   /**< nanoS_per_picoF */
    static const double sodium_calcium_exchanger_current_K_NaCa;   /**< picoA_per_picoF */
    static const double sodium_calcium_exchanger_current_K_sat;   /**< dimensionless */
    static const double sodium_calcium_exchanger_current_Km_Ca;   /**< millimolar */
    static const double sodium_calcium_exchanger_current_Km_Nai;   /**< millimolar */
    static const double sodium_calcium_exchanger_current_alpha;   /**< dimensionless */
    static const double sodium_calcium_exchanger_current_gamma;   /**< dimensionless */
    static const double sodium_dynamics_Na_o;   /**< millimolar */
    static const double sodium_potassium_pump_current_K_mNa;   /**< millimolar */
    static const double sodium_potassium_pump_current_K_mk;   /**< millimolar */
    static const double sodium_potassium_pump_current_P_NaK;   /**< picoA_per_picoF */
    static const double transient_outward_current_g_to;   /**< nanoS_per_picoF */

    /**
     *  This private method will check that gates are within 0 and 1 and concentrations are positive
     */
    void VerifyStateVariables();

public:
    /**
     * Constructor
     * 
     * @param pSolver is a pointer to the ODE solver
     * @param pIntracellularStimulus is a pointer to the intracellular stimulus
     */
    TenTusscher2006OdeSystem(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                             boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Destructor
     */
    ~TenTusscher2006OdeSystem();

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
     * Set the scale factor for Gks in order to differentiate epi M and endo cells
     * 
     * @param sfgks is the scale factor for Gks conductance
     */
    void SetScaleFactorGks(double sfgks);

    /**
     * Set the scale factor for Gks in order to differentiate epi M and endo cells
     * 
     * @param sfito is the scale factor for Ito current
     */
    void SetScaleFactorIto(double sfito);

    /**
     * Set the scale factor for Gks in order to differentiate epi M and endo cells
     * 
     * @param sfgkr is the scale factor for Gkr conductance
     */
    void SetScaleFactorGkr(double sfgkr);

     /**
     * Returns the ionic current
     * 
     * @return the total ionic current
     */
    double GetIIonic();

};

#endif // _TENTUSSCHER2006_HPP_
