/*

Copyright (C) University of Oxford, 2008

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
    static const double membrane_C = 75;//nanoF. It must be so because dV/dt is in mV/ms and currents are in nA.
    static const double radius = 0.05;//mm
    static const double length = 2;//mmm
    static const double V_e_ratio=0.1;//dimensionless
    static const double R=8314.472; //J/(kmol*Kelvin)
    static const double T=310; //Kelvin
    static const double F=96485.3415; //C/mol
    
    //conductances and the likes
    static const double g_fna = 3; //microS
    static const double g_fk = 3; //microS
    static const double g_na_b = 0.18; //microS
    static const double g_ca_b = 0.02; //microS
    static const double g_na = 750.0; //microS
    static const double g_k1 = 920; //microS
    static double const g_to = 0.28; //microS/mM
    static const double I_P = 125; //nA
    static const double i_kmax = 180;//nA
    static const double k_naca=0.02; //nA
    static const double P_si=15.0;//nA/mM
    //concentrations
    static const double Nao = 140; //mM
    static const double Cao = 2; //mM
    static const double Kb = 4;//mM
    static const double Kmf= 45; //mM
    static const double Km1= 210; //mM
    static const double Kmto= 10; //mM
    static const double KmK= 1; //mM
    static const double KmNa=40; //mM
    static const double Kmf2=0.001; //mM
    static const double KmCa=0.001;//mM
    static const double Km_Ca=0.0005;//mM
    //others
    static const double n_naca=3; //dimensionless
    static const double gamma=0.5; //dimensionless
    static const double d_naca=0.001; //dimensionless
    static const double rCa=2.0;//dimensionless
    static const double tau_up=0.025;//seconds. Must be seconds to cancel with numerator of F and give nA for i_up
    static const double tau_rep=2; //seconds. Must be seconds to cancel with numerator of F and give nA for i_rep
    static const double tau_rel=0.05; //seconds. Must be seconds to cancel with numerator of F and give nA for i_rel
    static const double Ca_up_max=5;//mM
    static const double pf=0.0007;//per millisecond
    static const double Vecs=0.05;//dimensionless
    //other parameters function of the above constants
    double Vi, Vup, Vrel, Ve;
    double Vcell;
    double RToNF;
    /** 
     *  Range-checking on the current values of the state variables. Make sure
     *  all gating variables have are within zero and one, and all concentrations
     *  are positive
     */
    void VerifyStateVariables();
    
public:
    // Constructor
    DiFrancescoNoble1985OdeSystem(AbstractIvpOdeSolver *pSolver,
                         AbstractStimulusFunction *pIntracellularStimulus);
                               
    // Destructor
    ~DiFrancescoNoble1985OdeSystem();
        
    // This method will compute the RHS of the Difrancesco-Noble model
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    
    double GetIIonic();
    
};

#endif // _DIFRANCESCONOBLE_HPP_
