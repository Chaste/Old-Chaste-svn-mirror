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

#ifndef _CML_noble_varghese_kohl_noble_1998_basic_pe_lut_
#define _CML_noble_varghese_kohl_noble_1998_basic_pe_lut_

// Model: noble_varghese_kohl_noble_1998_basic
// Processed by pycml - CellML Tools in Python
//     (translate: , pycml: )
// on Tue Jun 24 18:57:18 2008

#include <cmath>
#include <cassert>
#include "AbstractCardiacCell.hpp"
#include "Exception.hpp"
#include "AbstractStimulusFunction.hpp"

class CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables
{
public:
    static CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables* Instance()
    {
        if (mpInstance == NULL)
        {
            mpInstance = new CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables;
        }
        return mpInstance;
    }

    // Methods to look up values from lookup tables
    // using linear interpolation
    inline double _lookup_0(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][0];
        double y2 = _lookup_table_0[i+1][0];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_1(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][1];
        double y2 = _lookup_table_0[i+1][1];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_2(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][2];
        double y2 = _lookup_table_0[i+1][2];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_3(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][3];
        double y2 = _lookup_table_0[i+1][3];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_4(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][4];
        double y2 = _lookup_table_0[i+1][4];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_5(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][5];
        double y2 = _lookup_table_0[i+1][5];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_6(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][6];
        double y2 = _lookup_table_0[i+1][6];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_7(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][7];
        double y2 = _lookup_table_0[i+1][7];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_8(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][8];
        double y2 = _lookup_table_0[i+1][8];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_9(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][9];
        double y2 = _lookup_table_0[i+1][9];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_10(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][10];
        double y2 = _lookup_table_0[i+1][10];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_11(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][11];
        double y2 = _lookup_table_0[i+1][11];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_12(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][12];
        double y2 = _lookup_table_0[i+1][12];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_13(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][13];
        double y2 = _lookup_table_0[i+1][13];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_14(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][14];
        double y2 = _lookup_table_0[i+1][14];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_15(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][15];
        double y2 = _lookup_table_0[i+1][15];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_16(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][16];
        double y2 = _lookup_table_0[i+1][16];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_17(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][17];
        double y2 = _lookup_table_0[i+1][17];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_18(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][18];
        double y2 = _lookup_table_0[i+1][18];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_19(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][19];
        double y2 = _lookup_table_0[i+1][19];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_20(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][20];
        double y2 = _lookup_table_0[i+1][20];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_21(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][21];
        double y2 = _lookup_table_0[i+1][21];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_22(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][22];
        double y2 = _lookup_table_0[i+1][22];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_23(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][23];
        double y2 = _lookup_table_0[i+1][23];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_24(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][24];
        double y2 = _lookup_table_0[i+1][24];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_25(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][25];
        double y2 = _lookup_table_0[i+1][25];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_26(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][26];
        double y2 = _lookup_table_0[i+1][26];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_27(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][27];
        double y2 = _lookup_table_0[i+1][27];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_28(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][28];
        double y2 = _lookup_table_0[i+1][28];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_29(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][29];
        double y2 = _lookup_table_0[i+1][29];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_30(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][30];
        double y2 = _lookup_table_0[i+1][30];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_31(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][31];
        double y2 = _lookup_table_0[i+1][31];
        return y1 + (y2-y1)*factor;
    }

    
protected:
    CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables(const CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables&);
    CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables& operator= (const CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables&);
    CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables()
    {
        assert(mpInstance == NULL);
        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][0] = 1.0 + exp((var_membrane__V + 9.0) * 0.0446428571429);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][1] = 50.0 / (1.0 + exp((-(var_membrane__V - 5.0)) * 0.111111111111));
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][2] = 0.05 * exp((-(var_membrane__V - 20.0)) * 0.0666666666667);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][3] = 50.0 / (1.0 + exp((-(var_membrane__V - 5.0)) * 0.111111111111));
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][4] = 0.4 * exp(-pow((var_membrane__V + 30.0) * 0.0333333333333, 3.0));
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][5] = 14.0 / (1.0 + exp((-(var_membrane__V - 40.0)) * 0.111111111111));
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][6] = 1.0 * exp((-var_membrane__V) * 0.0222222222222);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][7] = 8000.0 * exp( -0.056 * (var_membrane__V + 66.0));
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][8] = 20.0 * exp( -0.125 * ((var_membrane__V + 75.0) - 0.0));
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][9] = 2000.0 / (1.0 + (320.0 * exp( -0.1 * ((var_membrane__V + 75.0) - 0.0))));
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][10] = 0.004 / (1.0 + exp((-(var_membrane__V + 52.0)) * 0.125));
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][11] = 1.0 - exp(((-(var_membrane__V - 50.0)) * 96485.3415 * 2.0) * 3.87974901066e-07);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][12] = 2.0 * exp(((-(var_membrane__V - 50.0)) * 96485.3415 * 2.0) * 3.87974901066e-07);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][13] = 1.0 - exp(((-(var_membrane__V - 50.0)) * 96485.3415) * 3.87974901066e-07);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][14] = 4.0 * exp(((-(var_membrane__V - 50.0)) * 96485.3415) * 3.87974901066e-07);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][15] = 1.0 - exp(((-(var_membrane__V - 50.0)) * 96485.3415) * 3.87974901066e-07);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][16] = 140.0 * exp(((-(var_membrane__V - 50.0)) * 96485.3415) * 3.87974901066e-07);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][17] = 1.0 - exp(((-(var_membrane__V - 50.0)) * 96485.3415 * 2.0) * 3.87974901066e-07);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][18] = 2.0 * exp(((-(var_membrane__V - 50.0)) * 96485.3415 * 2.0) * 3.87974901066e-07);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][19] = 1.0 - exp(((-(var_membrane__V - 50.0)) * 96485.3415) * 3.87974901066e-07);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][20] = 4.0 * exp(((-(var_membrane__V - 50.0)) * 96485.3415) * 3.87974901066e-07);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][21] = 1.0 - exp(((-(var_membrane__V - 50.0)) * 96485.3415) * 3.87974901066e-07);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][22] = 140.0 * exp(((-(var_membrane__V - 50.0)) * 96485.3415) * 3.87974901066e-07);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][23] = 12.0 / (1.0 + exp(( -1.0 * (var_membrane__V + 34.0)) * 0.25));
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][24] = 0.033 * exp((-var_membrane__V) * 0.0588235294118);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][25] = 33.0 / (1.0 + exp( -0.125 * (var_membrane__V + 10.0)));
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][26] = 1.0 / (1.0 + exp((-(var_membrane__V + 4.0)) * 0.2));
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][27] = exp((0.5 * 1.0 * var_membrane__V * 96485.3415) * 3.87974901066e-07);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][28] = exp(( -0.5 * 1.0 * var_membrane__V * 96485.3415) * 3.87974901066e-07);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][29] = exp((0.5 * 1.0 * var_membrane__V * 96485.3415) * 3.87974901066e-07);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][30] = exp(( -0.5 * 1.0 * var_membrane__V * 96485.3415) * 3.87974901066e-07);
        }

        for (int i=0 ; i<15001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][31] = 0.0 * exp(0.08 * (var_membrane__V - 40.0));
        }

    }

private:
    /** The single instance of the class */
    static CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables *mpInstance;
    // Lookup tables
    double _lookup_table_0[15001][32];
    
};

CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables* CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::mpInstance = NULL;

class CML_noble_varghese_kohl_noble_1998_basic_pe_lut : public AbstractCardiacCell
{
    // 
    // Settable parameters and readable variables
    // 
    double var_membrane__i_K1;
    double var_membrane__i_to;
    double var_membrane__i_Kr;
    double var_membrane__i_Ks;
    double var_membrane__i_Ca_L_K_cyt;
    double var_membrane__i_Ca_L_K_ds;
    double var_membrane__i_NaK;
    double var_membrane__i_Na;
    double var_membrane__i_b_Na;
    double var_membrane__i_p_Na;
    double var_membrane__i_Ca_L_Na_cyt;
    double var_membrane__i_Ca_L_Na_ds;
    double var_membrane__i_NaCa_cyt;
    double var_membrane__i_NaCa_ds;
    double var_membrane__i_Ca_L_Ca_cyt;
    double var_membrane__i_Ca_L_Ca_ds;
    double var_membrane__i_b_Ca;
    double var_membrane__i_Stim;
public:
    double Get_membrane__i_K1()
    {
        return var_membrane__i_K1;
    }

    double Get_membrane__i_to()
    {
        return var_membrane__i_to;
    }

    double Get_membrane__i_Kr()
    {
        return var_membrane__i_Kr;
    }

    double Get_membrane__i_Ks()
    {
        return var_membrane__i_Ks;
    }

    double Get_membrane__i_Ca_L_K_cyt()
    {
        return var_membrane__i_Ca_L_K_cyt;
    }

    double Get_membrane__i_Ca_L_K_ds()
    {
        return var_membrane__i_Ca_L_K_ds;
    }

    double Get_membrane__i_NaK()
    {
        return var_membrane__i_NaK;
    }

    double Get_membrane__i_Na()
    {
        return var_membrane__i_Na;
    }

    double Get_membrane__i_b_Na()
    {
        return var_membrane__i_b_Na;
    }

    double Get_membrane__i_p_Na()
    {
        return var_membrane__i_p_Na;
    }

    double Get_membrane__i_Ca_L_Na_cyt()
    {
        return var_membrane__i_Ca_L_Na_cyt;
    }

    double Get_membrane__i_Ca_L_Na_ds()
    {
        return var_membrane__i_Ca_L_Na_ds;
    }

    double Get_membrane__i_NaCa_cyt()
    {
        return var_membrane__i_NaCa_cyt;
    }

    double Get_membrane__i_NaCa_ds()
    {
        return var_membrane__i_NaCa_ds;
    }

    double Get_membrane__i_Ca_L_Ca_cyt()
    {
        return var_membrane__i_Ca_L_Ca_cyt;
    }

    double Get_membrane__i_Ca_L_Ca_ds()
    {
        return var_membrane__i_Ca_L_Ca_ds;
    }

    double Get_membrane__i_b_Ca()
    {
        return var_membrane__i_b_Ca;
    }

    double Get_membrane__i_Stim()
    {
        return var_membrane__i_Stim;
    }

    CML_noble_varghese_kohl_noble_1998_basic_pe_lut(AbstractIvpOdeSolver *pSolver, double dt,
                                                    AbstractStimulusFunction *pIntracellularStimulus,
                                                    AbstractStimulusFunction *pExtracellularStimulus=NULL)
        : AbstractCardiacCell(pSolver, 22, 0, dt, pIntracellularStimulus, pExtracellularStimulus)
    {
        // Time units: second
        // 
        mVariableNames.push_back("membrane__V");
        mVariableUnits.push_back("millivolt");
        mInitialConditions.push_back(-92.849333);

        mVariableNames.push_back("rapid_delayed_rectifier_potassium_current_xr1_gate__xr1");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(1.03e-5);

        mVariableNames.push_back("rapid_delayed_rectifier_potassium_current_xr2_gate__xr2");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(2e-7);

        mVariableNames.push_back("slow_delayed_rectifier_potassium_current_xs_gate__xs");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.001302);

        mVariableNames.push_back("fast_sodium_current_m_gate__m");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.0016203);

        mVariableNames.push_back("fast_sodium_current_h_gate__h");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.9944036);

        mVariableNames.push_back("L_type_Ca_channel_d_gate__d");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0);

        mVariableNames.push_back("L_type_Ca_channel_f_gate__f");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(1);

        mVariableNames.push_back("L_type_Ca_channel_f2_gate__f2");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.9349197);

        mVariableNames.push_back("L_type_Ca_channel_f2ds_gate__f2ds");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.9651958);

        mVariableNames.push_back("transient_outward_current_s_gate__s");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.9948645);

        mVariableNames.push_back("transient_outward_current_r_gate__r");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0);

        mVariableNames.push_back("calcium_release__ActFrac");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.0042614);

        mVariableNames.push_back("calcium_release__ProdFrac");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.4068154);

        mVariableNames.push_back("intracellular_sodium_concentration__Na_i");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(7.3321223);

        mVariableNames.push_back("intracellular_potassium_concentration__K_i");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(136.5644281);

        mVariableNames.push_back("intracellular_calcium_concentration__Ca_i");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(1.4e-5);

        mVariableNames.push_back("intracellular_calcium_concentration__Ca_ds");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(1.88e-5);

        mVariableNames.push_back("intracellular_calcium_concentration__Ca_up");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(0.4531889);

        mVariableNames.push_back("intracellular_calcium_concentration__Ca_rel");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(0.4481927);

        mVariableNames.push_back("intracellular_calcium_concentration__Ca_Calmod");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(0.0005555);

        mVariableNames.push_back("intracellular_calcium_concentration__Ca_Trop");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(0.0003542);

        Init();

    }

    ~CML_noble_varghese_kohl_noble_1998_basic_pe_lut(void)
    {
    }

    // Lookup table indices
    unsigned _table_index_0;
    double _factor_0;
    
    void VerifyGatingVariables() {}

    double GetIIonic()
    {
        std::vector<double>& rY = rGetStateVariables();
        double var_membrane__V = rY[0];
        // Units: millivolt; Initial value: -92.849333
        double var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1 = rY[1];
        // Units: dimensionless; Initial value: 1.03e-5
        double var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2 = rY[2];
        // Units: dimensionless; Initial value: 2e-7
        double var_slow_delayed_rectifier_potassium_current_xs_gate__xs = rY[3];
        // Units: dimensionless; Initial value: 0.001302
        double var_fast_sodium_current_m_gate__m = rY[4];
        // Units: dimensionless; Initial value: 0.0016203
        double var_fast_sodium_current_h_gate__h = rY[5];
        // Units: dimensionless; Initial value: 0.9944036
        double var_L_type_Ca_channel_d_gate__d = rY[6];
        // Units: dimensionless; Initial value: 0
        double var_L_type_Ca_channel_f_gate__f = rY[7];
        // Units: dimensionless; Initial value: 1
        double var_L_type_Ca_channel_f2_gate__f2 = rY[8];
        // Units: dimensionless; Initial value: 0.9349197
        double var_L_type_Ca_channel_f2ds_gate__f2ds = rY[9];
        // Units: dimensionless; Initial value: 0.9651958
        double var_transient_outward_current_s_gate__s = rY[10];
        // Units: dimensionless; Initial value: 0.9948645
        double var_transient_outward_current_r_gate__r = rY[11];
        // Units: dimensionless; Initial value: 0
        double var_intracellular_sodium_concentration__Na_i = rY[14];
        // Units: millimolar; Initial value: 7.3321223
        double var_intracellular_potassium_concentration__K_i = rY[15];
        // Units: millimolar; Initial value: 136.5644281
        double var_intracellular_calcium_concentration__Ca_i = rY[16];
        // Units: millimolar; Initial value: 1.4e-5
        double var_intracellular_calcium_concentration__Ca_ds = rY[17];
        // Units: millimolar; Initial value: 1.88e-5
        
        // Lookup table indexing
        if (var_membrane__V>49.9999 || var_membrane__V<-100.0001)
            EXCEPTION(DumpState("V outside lookup table range", rY));
        double _offset_0 = var_membrane__V - -100.0001;
        double _offset_0_over_table_step = _offset_0 * 100.0;
        unsigned _table_index_0 = (unsigned) floor(_offset_0_over_table_step);
        double _factor_0 = _offset_0_over_table_step - _table_index_0;
        
        double var_reversal_potentials__E_K = 26.7137606597 * log(4.0 / var_intracellular_potassium_concentration__K_i);
        double var_time_independent_potassium_current__i_K1 = (0.142857142857 * (var_membrane__V - var_reversal_potentials__E_K)) / (1.0 + exp((((var_membrane__V - var_reversal_potentials__E_K) - 10.0) * 96485.3415 * 1.25) * 3.87974901066e-07));
        var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
        double var_transient_outward_current__i_to = 0.005 * (0.0 + (var_transient_outward_current_s_gate__s * 1.0)) * var_transient_outward_current_r_gate__r * (var_membrane__V - var_reversal_potentials__E_K);
        var_membrane__i_to = var_transient_outward_current__i_to;
        double var_rapid_delayed_rectifier_potassium_current__i_Kr = ((((0.0021 * var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1) + (0.0013 * var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2)) * 1.0) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_0(_table_index_0, _factor_0)) * (var_membrane__V - var_reversal_potentials__E_K);
        var_membrane__i_Kr = var_rapid_delayed_rectifier_potassium_current__i_Kr;
        double var_slow_delayed_rectifier_potassium_current__i_Ks = 0.0026 * pow(var_slow_delayed_rectifier_potassium_current_xs_gate__xs, 2.0) * (var_membrane__V - (26.7137606597 * log(8.2 / (var_intracellular_potassium_concentration__K_i + (0.03 * var_intracellular_sodium_concentration__Na_i)))));
        var_membrane__i_Ks = var_slow_delayed_rectifier_potassium_current__i_Ks;
        double var_L_type_Ca_channel__i_Ca_L_K_cyt = (((0.0 * 0.002 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2_gate__f2 * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_13(_table_index_0, _factor_0)) * ((var_intracellular_potassium_concentration__K_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_14(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_K_cyt = var_L_type_Ca_channel__i_Ca_L_K_cyt;
        double var_L_type_Ca_channel__i_Ca_L_K_ds = (((1.0 * 0.002 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2ds_gate__f2ds * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_19(_table_index_0, _factor_0)) * ((var_intracellular_potassium_concentration__K_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_20(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_K_ds = var_L_type_Ca_channel__i_Ca_L_K_ds;
        double var_sodium_potassium_pump__i_NaK = (0.56 * var_intracellular_sodium_concentration__Na_i) / (40.0 + var_intracellular_sodium_concentration__Na_i);
        var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK;
        double var_fast_sodium_current__i_Na = 2.5 * pow(var_fast_sodium_current_m_gate__m, 3.0) * var_fast_sodium_current_h_gate__h * (var_membrane__V - (26.7137606597 * log(140.48 / (var_intracellular_sodium_concentration__Na_i + (0.12 * var_intracellular_potassium_concentration__K_i)))));
        var_membrane__i_Na = var_fast_sodium_current__i_Na;
        double var_reversal_potentials__E_Na = 26.7137606597 * log(140.0 / var_intracellular_sodium_concentration__Na_i);
        double var_sodium_background_current__i_b_Na = 0.0006 * (var_membrane__V - var_reversal_potentials__E_Na);
        var_membrane__i_b_Na = var_sodium_background_current__i_b_Na;
        double var_persistent_sodium_current__i_p_Na = CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_10(_table_index_0, _factor_0) * (var_membrane__V - var_reversal_potentials__E_Na);
        var_membrane__i_p_Na = var_persistent_sodium_current__i_p_Na;
        double var_L_type_Ca_channel__i_Ca_L_Na_cyt = (((0.0 * 0.01 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2_gate__f2 * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_15(_table_index_0, _factor_0)) * ((var_intracellular_sodium_concentration__Na_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_16(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_Na_cyt = var_L_type_Ca_channel__i_Ca_L_Na_cyt;
        double var_L_type_Ca_channel__i_Ca_L_Na_ds = (((1.0 * 0.01 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2ds_gate__f2ds * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_21(_table_index_0, _factor_0)) * ((var_intracellular_sodium_concentration__Na_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_22(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_Na_ds = var_L_type_Ca_channel__i_Ca_L_Na_ds;
        double var_sodium_calcium_exchanger__i_NaCa_cyt = (0.999 * 0.0005 * ((CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_27(_table_index_0, _factor_0) * pow(var_intracellular_sodium_concentration__Na_i, 3.0) * 2.0) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_28(_table_index_0, _factor_0) * 2744000.0 * var_intracellular_calcium_concentration__Ca_i))) / ((1.0 + (0.0 * ((var_intracellular_calcium_concentration__Ca_i * 2744000.0) + (2.0 * pow(var_intracellular_sodium_concentration__Na_i, 3.0))))) * (1.0 + (var_intracellular_calcium_concentration__Ca_i * 144.927536232)));
        var_membrane__i_NaCa_cyt = var_sodium_calcium_exchanger__i_NaCa_cyt;
        double var_sodium_calcium_exchanger__i_NaCa_ds = (0.001 * 0.0005 * ((CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_29(_table_index_0, _factor_0) * pow(var_intracellular_sodium_concentration__Na_i, 3.0) * 2.0) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_30(_table_index_0, _factor_0) * 2744000.0 * var_intracellular_calcium_concentration__Ca_ds))) / ((1.0 + (0.0 * ((var_intracellular_calcium_concentration__Ca_ds * 2744000.0) + (2.0 * pow(var_intracellular_sodium_concentration__Na_i, 3.0))))) * (1.0 + (var_intracellular_calcium_concentration__Ca_ds * 144.927536232)));
        var_membrane__i_NaCa_ds = var_sodium_calcium_exchanger__i_NaCa_ds;
        double var_L_type_Ca_channel__i_Ca_L_Ca_cyt = (((0.0 * 4.0 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2_gate__f2 * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_11(_table_index_0, _factor_0)) * ((var_intracellular_calcium_concentration__Ca_i * 42.2409057568) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_12(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_Ca_cyt = var_L_type_Ca_channel__i_Ca_L_Ca_cyt;
        double var_L_type_Ca_channel__i_Ca_L_Ca_ds = (((1.0 * 4.0 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2ds_gate__f2ds * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_17(_table_index_0, _factor_0)) * ((var_intracellular_calcium_concentration__Ca_i * 42.2409057568) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_18(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_Ca_ds = var_L_type_Ca_channel__i_Ca_L_Ca_ds;
        double var_calcium_background_current__i_b_Ca = 0.00025 * (var_membrane__V - (13.3568803298 * log(2.0 / var_intracellular_calcium_concentration__Ca_i)));
        var_membrane__i_b_Ca = var_calcium_background_current__i_b_Ca;
        
        return var_membrane__i_K1+var_membrane__i_to+var_membrane__i_Kr+var_membrane__i_Ks+var_membrane__i_Ca_L_K_cyt+var_membrane__i_Ca_L_K_ds+var_membrane__i_NaK+var_membrane__i_Na+var_membrane__i_b_Na+var_membrane__i_p_Na+var_membrane__i_Ca_L_Na_cyt+var_membrane__i_Ca_L_Na_ds+var_membrane__i_NaCa_cyt+var_membrane__i_NaCa_ds+var_membrane__i_Ca_L_Ca_cyt+var_membrane__i_Ca_L_Ca_ds+var_membrane__i_b_Ca;
    }

    void EvaluateYDerivatives (
            double var_environment__time,
            const std::vector<double> &rY,
            std::vector<double> &rDY)
    {
        // Inputs:
        // Time units: second
        double var_membrane__V = rY[0];
        // Units: millivolt; Initial value: -92.849333
        double var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1 = rY[1];
        // Units: dimensionless; Initial value: 1.03e-5
        double var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2 = rY[2];
        // Units: dimensionless; Initial value: 2e-7
        double var_slow_delayed_rectifier_potassium_current_xs_gate__xs = rY[3];
        // Units: dimensionless; Initial value: 0.001302
        double var_fast_sodium_current_m_gate__m = rY[4];
        // Units: dimensionless; Initial value: 0.0016203
        double var_fast_sodium_current_h_gate__h = rY[5];
        // Units: dimensionless; Initial value: 0.9944036
        double var_L_type_Ca_channel_d_gate__d = rY[6];
        // Units: dimensionless; Initial value: 0
        double var_L_type_Ca_channel_f_gate__f = rY[7];
        // Units: dimensionless; Initial value: 1
        double var_L_type_Ca_channel_f2_gate__f2 = rY[8];
        // Units: dimensionless; Initial value: 0.9349197
        double var_L_type_Ca_channel_f2ds_gate__f2ds = rY[9];
        // Units: dimensionless; Initial value: 0.9651958
        double var_transient_outward_current_s_gate__s = rY[10];
        // Units: dimensionless; Initial value: 0.9948645
        double var_transient_outward_current_r_gate__r = rY[11];
        // Units: dimensionless; Initial value: 0
        double var_calcium_release__ActFrac = rY[12];
        // Units: dimensionless; Initial value: 0.0042614
        double var_calcium_release__ProdFrac = rY[13];
        // Units: dimensionless; Initial value: 0.4068154
        double var_intracellular_sodium_concentration__Na_i = rY[14];
        // Units: millimolar; Initial value: 7.3321223
        double var_intracellular_potassium_concentration__K_i = rY[15];
        // Units: millimolar; Initial value: 136.5644281
        double var_intracellular_calcium_concentration__Ca_i = rY[16];
        // Units: millimolar; Initial value: 1.4e-5
        double var_intracellular_calcium_concentration__Ca_ds = rY[17];
        // Units: millimolar; Initial value: 1.88e-5
        double var_intracellular_calcium_concentration__Ca_up = rY[18];
        // Units: millimolar; Initial value: 0.4531889
        double var_intracellular_calcium_concentration__Ca_rel = rY[19];
        // Units: millimolar; Initial value: 0.4481927
        double var_intracellular_calcium_concentration__Ca_Calmod = rY[20];
        // Units: millimolar; Initial value: 0.0005555
        double var_intracellular_calcium_concentration__Ca_Trop = rY[21];
        // Units: millimolar; Initial value: 0.0003542
        
        
        // Lookup table indexing
        if (var_membrane__V>49.9999 || var_membrane__V<-100.0001)
            EXCEPTION(DumpState("V outside lookup table range", rY));
        double _offset_0 = var_membrane__V - -100.0001;
        double _offset_0_over_table_step = _offset_0 * 100.0;
        unsigned _table_index_0 = (unsigned) floor(_offset_0_over_table_step);
        double _factor_0 = _offset_0_over_table_step - _table_index_0;
        
        // Mathematics
        double var_reversal_potentials__E_K = 26.7137606597 * log(4.0 / var_intracellular_potassium_concentration__K_i);
        double var_time_independent_potassium_current__i_K1 = (0.142857142857 * (var_membrane__V - var_reversal_potentials__E_K)) / (1.0 + exp((((var_membrane__V - var_reversal_potentials__E_K) - 10.0) * 96485.3415 * 1.25) * 3.87974901066e-07));
        var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
        double var_transient_outward_current__i_to = 0.005 * (0.0 + (var_transient_outward_current_s_gate__s * 1.0)) * var_transient_outward_current_r_gate__r * (var_membrane__V - var_reversal_potentials__E_K);
        var_membrane__i_to = var_transient_outward_current__i_to;
        double var_rapid_delayed_rectifier_potassium_current__i_Kr = ((((0.0021 * var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1) + (0.0013 * var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2)) * 1.0) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_0(_table_index_0, _factor_0)) * (var_membrane__V - var_reversal_potentials__E_K);
        var_membrane__i_Kr = var_rapid_delayed_rectifier_potassium_current__i_Kr;
        double var_slow_delayed_rectifier_potassium_current__i_Ks = 0.0026 * pow(var_slow_delayed_rectifier_potassium_current_xs_gate__xs, 2.0) * (var_membrane__V - (26.7137606597 * log(8.2 / (var_intracellular_potassium_concentration__K_i + (0.03 * var_intracellular_sodium_concentration__Na_i)))));
        var_membrane__i_Ks = var_slow_delayed_rectifier_potassium_current__i_Ks;
        double var_L_type_Ca_channel__i_Ca_L_K_cyt = (((0.0 * 0.002 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2_gate__f2 * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_13(_table_index_0, _factor_0)) * ((var_intracellular_potassium_concentration__K_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_14(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_K_cyt = var_L_type_Ca_channel__i_Ca_L_K_cyt;
        double var_L_type_Ca_channel__i_Ca_L_K_ds = (((1.0 * 0.002 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2ds_gate__f2ds * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_19(_table_index_0, _factor_0)) * ((var_intracellular_potassium_concentration__K_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_20(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_K_ds = var_L_type_Ca_channel__i_Ca_L_K_ds;
        double var_sodium_potassium_pump__i_NaK = (0.56 * var_intracellular_sodium_concentration__Na_i) / (40.0 + var_intracellular_sodium_concentration__Na_i);
        var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK;
        double var_fast_sodium_current__i_Na = 2.5 * pow(var_fast_sodium_current_m_gate__m, 3.0) * var_fast_sodium_current_h_gate__h * (var_membrane__V - (26.7137606597 * log(140.48 / (var_intracellular_sodium_concentration__Na_i + (0.12 * var_intracellular_potassium_concentration__K_i)))));
        var_membrane__i_Na = var_fast_sodium_current__i_Na;
        double var_reversal_potentials__E_Na = 26.7137606597 * log(140.0 / var_intracellular_sodium_concentration__Na_i);
        double var_sodium_background_current__i_b_Na = 0.0006 * (var_membrane__V - var_reversal_potentials__E_Na);
        var_membrane__i_b_Na = var_sodium_background_current__i_b_Na;
        double var_persistent_sodium_current__i_p_Na = CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_10(_table_index_0, _factor_0) * (var_membrane__V - var_reversal_potentials__E_Na);
        var_membrane__i_p_Na = var_persistent_sodium_current__i_p_Na;
        double var_L_type_Ca_channel__i_Ca_L_Na_cyt = (((0.0 * 0.01 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2_gate__f2 * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_15(_table_index_0, _factor_0)) * ((var_intracellular_sodium_concentration__Na_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_16(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_Na_cyt = var_L_type_Ca_channel__i_Ca_L_Na_cyt;
        double var_L_type_Ca_channel__i_Ca_L_Na_ds = (((1.0 * 0.01 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2ds_gate__f2ds * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_21(_table_index_0, _factor_0)) * ((var_intracellular_sodium_concentration__Na_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_22(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_Na_ds = var_L_type_Ca_channel__i_Ca_L_Na_ds;
        double var_sodium_calcium_exchanger__i_NaCa_cyt = (0.999 * 0.0005 * ((CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_27(_table_index_0, _factor_0) * pow(var_intracellular_sodium_concentration__Na_i, 3.0) * 2.0) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_28(_table_index_0, _factor_0) * 2744000.0 * var_intracellular_calcium_concentration__Ca_i))) / ((1.0 + (0.0 * ((var_intracellular_calcium_concentration__Ca_i * 2744000.0) + (2.0 * pow(var_intracellular_sodium_concentration__Na_i, 3.0))))) * (1.0 + (var_intracellular_calcium_concentration__Ca_i * 144.927536232)));
        var_membrane__i_NaCa_cyt = var_sodium_calcium_exchanger__i_NaCa_cyt;
        double var_sodium_calcium_exchanger__i_NaCa_ds = (0.001 * 0.0005 * ((CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_29(_table_index_0, _factor_0) * pow(var_intracellular_sodium_concentration__Na_i, 3.0) * 2.0) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_30(_table_index_0, _factor_0) * 2744000.0 * var_intracellular_calcium_concentration__Ca_ds))) / ((1.0 + (0.0 * ((var_intracellular_calcium_concentration__Ca_ds * 2744000.0) + (2.0 * pow(var_intracellular_sodium_concentration__Na_i, 3.0))))) * (1.0 + (var_intracellular_calcium_concentration__Ca_ds * 144.927536232)));
        var_membrane__i_NaCa_ds = var_sodium_calcium_exchanger__i_NaCa_ds;
        double var_L_type_Ca_channel__i_Ca_L_Ca_cyt = (((0.0 * 4.0 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2_gate__f2 * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_11(_table_index_0, _factor_0)) * ((var_intracellular_calcium_concentration__Ca_i * 42.2409057568) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_12(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_Ca_cyt = var_L_type_Ca_channel__i_Ca_L_Ca_cyt;
        double var_L_type_Ca_channel__i_Ca_L_Ca_ds = (((1.0 * 4.0 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2ds_gate__f2ds * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_17(_table_index_0, _factor_0)) * ((var_intracellular_calcium_concentration__Ca_i * 42.2409057568) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_18(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_Ca_ds = var_L_type_Ca_channel__i_Ca_L_Ca_ds;
        double var_calcium_background_current__i_b_Ca = 0.00025 * (var_membrane__V - (13.3568803298 * log(2.0 / var_intracellular_calcium_concentration__Ca_i)));
        var_membrane__i_b_Ca = var_calcium_background_current__i_b_Ca;
        double var_membrane__i_Stim = GetStimulus(var_environment__time);
        double var_fast_sodium_current_m_gate__E0_m = var_membrane__V + 41.0;
        double var_L_type_Ca_channel_d_gate__E0_d = (var_membrane__V + 24.0) - 5.0;
        double var_L_type_Ca_channel_f_gate__E0_f = var_membrane__V + 34.0;
        double var_sarcoplasmic_reticulum_calcium_pump__K_2 = var_intracellular_calcium_concentration__Ca_i + (var_intracellular_calcium_concentration__Ca_up * 0.00024) + 0.00012 + 0.0003;
        double var_sarcoplasmic_reticulum_calcium_pump__i_up = ((var_intracellular_calcium_concentration__Ca_i / var_sarcoplasmic_reticulum_calcium_pump__K_2) * 0.4) - (((var_intracellular_calcium_concentration__Ca_up * 0.00024) / var_sarcoplasmic_reticulum_calcium_pump__K_2) * 0.03);
        double var_calcium_translocation__i_trans = 50.0 * (var_intracellular_calcium_concentration__Ca_up - var_intracellular_calcium_concentration__Ca_rel);
        double var_calcium_release__i_rel = ((pow(var_calcium_release__ActFrac / (var_calcium_release__ActFrac + 0.25), 2.0) * 250.0) + 0.05) * var_intracellular_calcium_concentration__Ca_rel;
        double var_calcium_release__CaiReg = var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005);
        double var_calcium_release__RegBindSite = var_calcium_release__CaiReg + ((1.0 - var_calcium_release__CaiReg) * (var_intracellular_calcium_concentration__Ca_ds / (var_intracellular_calcium_concentration__Ca_ds + 0.01)));
        double var_calcium_release__InactRate = 60.0 + (500.0 * pow(var_calcium_release__RegBindSite, 2.0));
        double var_calcium_release__SpeedRel = (var_membrane__V <  -50.0) ? 5.0 : 1.0;
        double d_dt_membrane__V =  -10526.3157895 * (var_membrane__i_Stim + var_membrane__i_K1 + var_membrane__i_to + var_membrane__i_Kr + var_membrane__i_Ks + var_membrane__i_NaK + var_membrane__i_Na + var_membrane__i_b_Na + var_membrane__i_p_Na + var_membrane__i_Ca_L_Na_cyt + var_membrane__i_Ca_L_Na_ds + var_membrane__i_NaCa_cyt + var_membrane__i_NaCa_ds + var_membrane__i_Ca_L_Ca_cyt + var_membrane__i_Ca_L_Ca_ds + var_membrane__i_Ca_L_K_cyt + var_membrane__i_Ca_L_K_ds + var_membrane__i_b_Ca);
        double d_dt_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1 = (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_1(_table_index_0, _factor_0) * (1.0 - var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1)) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_2(_table_index_0, _factor_0) * var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1);
        double d_dt_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2 = (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_3(_table_index_0, _factor_0) * (1.0 - var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2)) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_4(_table_index_0, _factor_0) * var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2);
        double d_dt_slow_delayed_rectifier_potassium_current_xs_gate__xs = (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_5(_table_index_0, _factor_0) * (1.0 - var_slow_delayed_rectifier_potassium_current_xs_gate__xs)) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_6(_table_index_0, _factor_0) * var_slow_delayed_rectifier_potassium_current_xs_gate__xs);
        double d_dt_fast_sodium_current_m_gate__m = (((fabs(var_fast_sodium_current_m_gate__E0_m) < 1e-05) ? 2000.0 : ((200.0 * var_fast_sodium_current_m_gate__E0_m) / (1.0 - exp( -0.1 * var_fast_sodium_current_m_gate__E0_m)))) * (1.0 - var_fast_sodium_current_m_gate__m)) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_7(_table_index_0, _factor_0) * var_fast_sodium_current_m_gate__m);
        double d_dt_fast_sodium_current_h_gate__h = (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_8(_table_index_0, _factor_0) * (1.0 - var_fast_sodium_current_h_gate__h)) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_9(_table_index_0, _factor_0) * var_fast_sodium_current_h_gate__h);
        double d_dt_L_type_Ca_channel_d_gate__d = 3.0 * ((((fabs(var_L_type_Ca_channel_d_gate__E0_d) < 0.0001) ? 120.0 : ((30.0 * var_L_type_Ca_channel_d_gate__E0_d) / (1.0 - exp((-var_L_type_Ca_channel_d_gate__E0_d) * 0.25)))) * (1.0 - var_L_type_Ca_channel_d_gate__d)) - (((fabs(var_L_type_Ca_channel_d_gate__E0_d) < 0.0001) ? 120.0 : ((12.0 * var_L_type_Ca_channel_d_gate__E0_d) / (exp(var_L_type_Ca_channel_d_gate__E0_d * 0.1) - 1.0))) * var_L_type_Ca_channel_d_gate__d));
        double d_dt_L_type_Ca_channel_f_gate__f = 0.3 * ((((fabs(var_L_type_Ca_channel_f_gate__E0_f) < 0.0001) ? 25.0 : ((6.25 * var_L_type_Ca_channel_f_gate__E0_f) / (exp(var_L_type_Ca_channel_f_gate__E0_f * 0.25) - 1.0))) * (1.0 - var_L_type_Ca_channel_f_gate__f)) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_23(_table_index_0, _factor_0) * var_L_type_Ca_channel_f_gate__f));
        double d_dt_L_type_Ca_channel_f2_gate__f2 = 1.0 - (1.0 * ((var_intracellular_calcium_concentration__Ca_i / (100000.0 + var_intracellular_calcium_concentration__Ca_i)) + var_L_type_Ca_channel_f2_gate__f2));
        double d_dt_L_type_Ca_channel_f2ds_gate__f2ds = 20.0 * (1.0 - ((var_intracellular_calcium_concentration__Ca_ds / (0.001 + var_intracellular_calcium_concentration__Ca_ds)) + var_L_type_Ca_channel_f2ds_gate__f2ds));
        double d_dt_transient_outward_current_s_gate__s = (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_24(_table_index_0, _factor_0) * (1.0 - var_transient_outward_current_s_gate__s)) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_25(_table_index_0, _factor_0) * var_transient_outward_current_s_gate__s);
        double d_dt_transient_outward_current_r_gate__r = 333.0 * (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_26(_table_index_0, _factor_0) - var_transient_outward_current_r_gate__r);
        double d_dt_calcium_release__ActFrac = (((1.0 - var_calcium_release__ActFrac) - var_calcium_release__ProdFrac) * var_calcium_release__SpeedRel * (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_LookupTables::Instance()->_lookup_31(_table_index_0, _factor_0) + (500.0 * pow(var_calcium_release__RegBindSite, 2.0)))) - (var_calcium_release__ActFrac * var_calcium_release__SpeedRel * var_calcium_release__InactRate);
        double d_dt_calcium_release__ProdFrac = (var_calcium_release__ActFrac * var_calcium_release__SpeedRel * var_calcium_release__InactRate) - (var_calcium_release__SpeedRel * 1.0 * var_calcium_release__ProdFrac);
        double d_dt_intracellular_sodium_concentration__Na_i =  -0.631827460168 * (var_fast_sodium_current__i_Na + var_persistent_sodium_current__i_p_Na + var_sodium_background_current__i_b_Na + (3.0 * var_sodium_potassium_pump__i_NaK) + (3.0 * var_sodium_calcium_exchanger__i_NaCa_cyt) + var_L_type_Ca_channel__i_Ca_L_Na_cyt + var_L_type_Ca_channel__i_Ca_L_Na_ds);
        double d_dt_intracellular_potassium_concentration__K_i =  -0.631827460168 * ((var_time_independent_potassium_current__i_K1 + var_rapid_delayed_rectifier_potassium_current__i_Kr + var_slow_delayed_rectifier_potassium_current__i_Ks + var_L_type_Ca_channel__i_Ca_L_K_cyt + var_L_type_Ca_channel__i_Ca_L_K_ds + var_transient_outward_current__i_to) - (2.0 * var_sodium_potassium_pump__i_NaK));
        double d_dt_intracellular_calcium_concentration__Ca_Trop = (100000.0 * var_intracellular_calcium_concentration__Ca_i * (0.05 - var_intracellular_calcium_concentration__Ca_Trop)) - (200.0 * var_intracellular_calcium_concentration__Ca_Trop);
        double d_dt_intracellular_calcium_concentration__Ca_Calmod = (100000.0 * var_intracellular_calcium_concentration__Ca_i * (0.02 - var_intracellular_calcium_concentration__Ca_Calmod)) - (50.0 * var_intracellular_calcium_concentration__Ca_Calmod);
        double d_dt_intracellular_calcium_concentration__Ca_i = (((( -0.315913730084 * (((var_L_type_Ca_channel__i_Ca_L_Ca_cyt + var_calcium_background_current__i_b_Ca) - (2.0 * var_sodium_calcium_exchanger__i_NaCa_cyt)) - (2.0 * var_sodium_calcium_exchanger__i_NaCa_ds))) + (var_intracellular_calcium_concentration__Ca_ds * 0.1 * 10.0) + ((var_calcium_release__i_rel * 0.1) * 2.04081632653)) - d_dt_intracellular_calcium_concentration__Ca_Calmod) - d_dt_intracellular_calcium_concentration__Ca_Trop) - var_sarcoplasmic_reticulum_calcium_pump__i_up;
        double d_dt_intracellular_calcium_concentration__Ca_ds = (( -1.0 * var_L_type_Ca_channel__i_Ca_L_Ca_ds) * 3.15913730084) - (var_intracellular_calcium_concentration__Ca_ds * 10.0);
        double d_dt_intracellular_calcium_concentration__Ca_up = (49.0 * var_sarcoplasmic_reticulum_calcium_pump__i_up) - var_calcium_translocation__i_trans;
        double d_dt_intracellular_calcium_concentration__Ca_rel = (0.1 * var_calcium_translocation__i_trans) - var_calcium_release__i_rel;
        
        rDY[0] = d_dt_membrane__V;
        rDY[1] = d_dt_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1;
        rDY[2] = d_dt_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2;
        rDY[3] = d_dt_slow_delayed_rectifier_potassium_current_xs_gate__xs;
        rDY[4] = d_dt_fast_sodium_current_m_gate__m;
        rDY[5] = d_dt_fast_sodium_current_h_gate__h;
        rDY[6] = d_dt_L_type_Ca_channel_d_gate__d;
        rDY[7] = d_dt_L_type_Ca_channel_f_gate__f;
        rDY[8] = d_dt_L_type_Ca_channel_f2_gate__f2;
        rDY[9] = d_dt_L_type_Ca_channel_f2ds_gate__f2ds;
        rDY[10] = d_dt_transient_outward_current_s_gate__s;
        rDY[11] = d_dt_transient_outward_current_r_gate__r;
        rDY[12] = d_dt_calcium_release__ActFrac;
        rDY[13] = d_dt_calcium_release__ProdFrac;
        rDY[14] = d_dt_intracellular_sodium_concentration__Na_i;
        rDY[15] = d_dt_intracellular_potassium_concentration__K_i;
        rDY[16] = d_dt_intracellular_calcium_concentration__Ca_i;
        rDY[17] = d_dt_intracellular_calcium_concentration__Ca_ds;
        rDY[18] = d_dt_intracellular_calcium_concentration__Ca_up;
        rDY[19] = d_dt_intracellular_calcium_concentration__Ca_rel;
        rDY[20] = d_dt_intracellular_calcium_concentration__Ca_Calmod;
        rDY[21] = d_dt_intracellular_calcium_concentration__Ca_Trop;
    }

};

#endif
