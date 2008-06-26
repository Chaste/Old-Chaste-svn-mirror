#ifndef _CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_
#define _CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_

// Model: noble_varghese_kohl_noble_1998_basic
// Processed by pycml - CellML Tools in Python
//     (translate: , pycml: )
// on Tue Jun 24 18:55:53 2008

#include <cmath>
#include <cassert>
#include "AbstractBackwardEulerCardiacCell.hpp"
#include "CardiacNewtonSolver.hpp"
#include "Exception.hpp"
#include "AbstractStimulusFunction.hpp"

class CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables
{
public:
    static CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables* Instance()
    {
        if (mpInstance == NULL)
        {
            mpInstance = new CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables;
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
    CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables(const CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables&);
    CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables& operator= (const CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables&);
    CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables()
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
    static CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables *mpInstance;
    // Lookup tables
    double _lookup_table_0[15001][32];
    
};

CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables* CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::mpInstance = NULL;

class CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be : public AbstractBackwardEulerCardiacCell<12>
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

    CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be(AbstractStimulusFunction *pIntracellularStimulus,
                                                       AbstractStimulusFunction *pExtracellularStimulus=NULL)
        : AbstractBackwardEulerCardiacCell<12>(22, 0, pIntracellularStimulus, pExtracellularStimulus)
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

    ~CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be(void)
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
        _table_index_0 = (unsigned) floor(_offset_0_over_table_step);
        _factor_0 = _offset_0_over_table_step - _table_index_0;
        
        double var_reversal_potentials__E_K = 26.7137606597 * log(4.0 / var_intracellular_potassium_concentration__K_i);
        double var_time_independent_potassium_current__i_K1 = (0.142857142857 * (var_membrane__V - var_reversal_potentials__E_K)) / (1.0 + exp((((var_membrane__V - var_reversal_potentials__E_K) - 10.0) * 96485.3415 * 1.25) * 3.87974901066e-07));
        var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
        double var_transient_outward_current__i_to = 0.005 * (0.0 + (var_transient_outward_current_s_gate__s * 1.0)) * var_transient_outward_current_r_gate__r * (var_membrane__V - var_reversal_potentials__E_K);
        var_membrane__i_to = var_transient_outward_current__i_to;
        double var_rapid_delayed_rectifier_potassium_current__i_Kr = ((((0.0021 * var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1) + (0.0013 * var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2)) * 1.0) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_0(_table_index_0, _factor_0)) * (var_membrane__V - var_reversal_potentials__E_K);
        var_membrane__i_Kr = var_rapid_delayed_rectifier_potassium_current__i_Kr;
        double var_slow_delayed_rectifier_potassium_current__i_Ks = 0.0026 * pow(var_slow_delayed_rectifier_potassium_current_xs_gate__xs, 2.0) * (var_membrane__V - (26.7137606597 * log(8.2 / (var_intracellular_potassium_concentration__K_i + (0.03 * var_intracellular_sodium_concentration__Na_i)))));
        var_membrane__i_Ks = var_slow_delayed_rectifier_potassium_current__i_Ks;
        double var_L_type_Ca_channel__i_Ca_L_K_cyt = (((0.0 * 0.002 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2_gate__f2 * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_13(_table_index_0, _factor_0)) * ((var_intracellular_potassium_concentration__K_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_14(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_K_cyt = var_L_type_Ca_channel__i_Ca_L_K_cyt;
        double var_L_type_Ca_channel__i_Ca_L_K_ds = (((1.0 * 0.002 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2ds_gate__f2ds * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_19(_table_index_0, _factor_0)) * ((var_intracellular_potassium_concentration__K_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_20(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_K_ds = var_L_type_Ca_channel__i_Ca_L_K_ds;
        double var_sodium_potassium_pump__i_NaK = (0.56 * var_intracellular_sodium_concentration__Na_i) / (40.0 + var_intracellular_sodium_concentration__Na_i);
        var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK;
        double var_fast_sodium_current__i_Na = 2.5 * pow(var_fast_sodium_current_m_gate__m, 3.0) * var_fast_sodium_current_h_gate__h * (var_membrane__V - (26.7137606597 * log(140.48 / (var_intracellular_sodium_concentration__Na_i + (0.12 * var_intracellular_potassium_concentration__K_i)))));
        var_membrane__i_Na = var_fast_sodium_current__i_Na;
        double var_reversal_potentials__E_Na = 26.7137606597 * log(140.0 / var_intracellular_sodium_concentration__Na_i);
        double var_sodium_background_current__i_b_Na = 0.0006 * (var_membrane__V - var_reversal_potentials__E_Na);
        var_membrane__i_b_Na = var_sodium_background_current__i_b_Na;
        double var_persistent_sodium_current__i_p_Na = CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_10(_table_index_0, _factor_0) * (var_membrane__V - var_reversal_potentials__E_Na);
        var_membrane__i_p_Na = var_persistent_sodium_current__i_p_Na;
        double var_L_type_Ca_channel__i_Ca_L_Na_cyt = (((0.0 * 0.01 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2_gate__f2 * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_15(_table_index_0, _factor_0)) * ((var_intracellular_sodium_concentration__Na_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_16(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_Na_cyt = var_L_type_Ca_channel__i_Ca_L_Na_cyt;
        double var_L_type_Ca_channel__i_Ca_L_Na_ds = (((1.0 * 0.01 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2ds_gate__f2ds * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_21(_table_index_0, _factor_0)) * ((var_intracellular_sodium_concentration__Na_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_22(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_Na_ds = var_L_type_Ca_channel__i_Ca_L_Na_ds;
        double var_sodium_calcium_exchanger__i_NaCa_cyt = (0.999 * 0.0005 * ((CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_27(_table_index_0, _factor_0) * pow(var_intracellular_sodium_concentration__Na_i, 3.0) * 2.0) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_28(_table_index_0, _factor_0) * 2744000.0 * var_intracellular_calcium_concentration__Ca_i))) / ((1.0 + (0.0 * ((var_intracellular_calcium_concentration__Ca_i * 2744000.0) + (2.0 * pow(var_intracellular_sodium_concentration__Na_i, 3.0))))) * (1.0 + (var_intracellular_calcium_concentration__Ca_i * 144.927536232)));
        var_membrane__i_NaCa_cyt = var_sodium_calcium_exchanger__i_NaCa_cyt;
        double var_sodium_calcium_exchanger__i_NaCa_ds = (0.001 * 0.0005 * ((CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_29(_table_index_0, _factor_0) * pow(var_intracellular_sodium_concentration__Na_i, 3.0) * 2.0) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_30(_table_index_0, _factor_0) * 2744000.0 * var_intracellular_calcium_concentration__Ca_ds))) / ((1.0 + (0.0 * ((var_intracellular_calcium_concentration__Ca_ds * 2744000.0) + (2.0 * pow(var_intracellular_sodium_concentration__Na_i, 3.0))))) * (1.0 + (var_intracellular_calcium_concentration__Ca_ds * 144.927536232)));
        var_membrane__i_NaCa_ds = var_sodium_calcium_exchanger__i_NaCa_ds;
        double var_L_type_Ca_channel__i_Ca_L_Ca_cyt = (((0.0 * 4.0 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2_gate__f2 * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_11(_table_index_0, _factor_0)) * ((var_intracellular_calcium_concentration__Ca_i * 42.2409057568) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_12(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_Ca_cyt = var_L_type_Ca_channel__i_Ca_L_Ca_cyt;
        double var_L_type_Ca_channel__i_Ca_L_Ca_ds = (((1.0 * 4.0 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2ds_gate__f2ds * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_17(_table_index_0, _factor_0)) * ((var_intracellular_calcium_concentration__Ca_i * 42.2409057568) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_18(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_Ca_ds = var_L_type_Ca_channel__i_Ca_L_Ca_ds;
        double var_calcium_background_current__i_b_Ca = 0.00025 * (var_membrane__V - (13.3568803298 * log(2.0 / var_intracellular_calcium_concentration__Ca_i)));
        var_membrane__i_b_Ca = var_calcium_background_current__i_b_Ca;
        
        return var_membrane__i_K1+var_membrane__i_to+var_membrane__i_Kr+var_membrane__i_Ks+var_membrane__i_Ca_L_K_cyt+var_membrane__i_Ca_L_K_ds+var_membrane__i_NaK+var_membrane__i_Na+var_membrane__i_b_Na+var_membrane__i_p_Na+var_membrane__i_Ca_L_Na_cyt+var_membrane__i_Ca_L_Na_ds+var_membrane__i_NaCa_cyt+var_membrane__i_NaCa_ds+var_membrane__i_Ca_L_Ca_cyt+var_membrane__i_Ca_L_Ca_ds+var_membrane__i_b_Ca;
    }

    void ComputeResidual(const double rCurrentGuess[12], double rResidual[12])
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
        double var_transient_outward_current_s_gate__s = rY[10];
        // Units: dimensionless; Initial value: 0.9948645
        double var_transient_outward_current_r_gate__r = rY[11];
        // Units: dimensionless; Initial value: 0
        
        double var_L_type_Ca_channel_f2_gate__f2 = rCurrentGuess[0];
        double var_L_type_Ca_channel_f2ds_gate__f2ds = rCurrentGuess[1];
        double var_calcium_release__ActFrac = rCurrentGuess[2];
        double var_calcium_release__ProdFrac = rCurrentGuess[3];
        double var_intracellular_calcium_concentration__Ca_Calmod = rCurrentGuess[4];
        double var_intracellular_calcium_concentration__Ca_Trop = rCurrentGuess[5];
        double var_intracellular_calcium_concentration__Ca_ds = rCurrentGuess[6];
        double var_intracellular_calcium_concentration__Ca_i = rCurrentGuess[7];
        double var_intracellular_calcium_concentration__Ca_rel = rCurrentGuess[8];
        double var_intracellular_calcium_concentration__Ca_up = rCurrentGuess[9];
        double var_intracellular_potassium_concentration__K_i = rCurrentGuess[10];
        double var_intracellular_sodium_concentration__Na_i = rCurrentGuess[11];
        
        double var_reversal_potentials__E_K = 26.7137606597 * log(4.0 / var_intracellular_potassium_concentration__K_i);
        double var_time_independent_potassium_current__i_K1 = (0.142857142857 * (var_membrane__V - var_reversal_potentials__E_K)) / (1.0 + exp((((var_membrane__V - var_reversal_potentials__E_K) - 10.0) * 96485.3415 * 1.25) * 3.87974901066e-07));
        double var_transient_outward_current__i_to = 0.005 * (0.0 + (var_transient_outward_current_s_gate__s * 1.0)) * var_transient_outward_current_r_gate__r * (var_membrane__V - var_reversal_potentials__E_K);
        double var_rapid_delayed_rectifier_potassium_current__i_Kr = ((((0.0021 * var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1) + (0.0013 * var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2)) * 1.0) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_0(_table_index_0, _factor_0)) * (var_membrane__V - var_reversal_potentials__E_K);
        double var_slow_delayed_rectifier_potassium_current__i_Ks = 0.0026 * pow(var_slow_delayed_rectifier_potassium_current_xs_gate__xs, 2.0) * (var_membrane__V - (26.7137606597 * log(8.2 / (var_intracellular_potassium_concentration__K_i + (0.03 * var_intracellular_sodium_concentration__Na_i)))));
        double var_L_type_Ca_channel__i_Ca_L_K_cyt = (((0.0 * 0.002 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2_gate__f2 * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_13(_table_index_0, _factor_0)) * ((var_intracellular_potassium_concentration__K_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_14(_table_index_0, _factor_0));
        double var_L_type_Ca_channel__i_Ca_L_K_ds = (((1.0 * 0.002 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2ds_gate__f2ds * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_19(_table_index_0, _factor_0)) * ((var_intracellular_potassium_concentration__K_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_20(_table_index_0, _factor_0));
        double var_sodium_potassium_pump__i_NaK = (0.56 * var_intracellular_sodium_concentration__Na_i) / (40.0 + var_intracellular_sodium_concentration__Na_i);
        double var_fast_sodium_current__i_Na = 2.5 * pow(var_fast_sodium_current_m_gate__m, 3.0) * var_fast_sodium_current_h_gate__h * (var_membrane__V - (26.7137606597 * log(140.48 / (var_intracellular_sodium_concentration__Na_i + (0.12 * var_intracellular_potassium_concentration__K_i)))));
        double var_reversal_potentials__E_Na = 26.7137606597 * log(140.0 / var_intracellular_sodium_concentration__Na_i);
        double var_sodium_background_current__i_b_Na = 0.0006 * (var_membrane__V - var_reversal_potentials__E_Na);
        double var_persistent_sodium_current__i_p_Na = CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_10(_table_index_0, _factor_0) * (var_membrane__V - var_reversal_potentials__E_Na);
        double var_L_type_Ca_channel__i_Ca_L_Na_cyt = (((0.0 * 0.01 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2_gate__f2 * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_15(_table_index_0, _factor_0)) * ((var_intracellular_sodium_concentration__Na_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_16(_table_index_0, _factor_0));
        double var_L_type_Ca_channel__i_Ca_L_Na_ds = (((1.0 * 0.01 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2ds_gate__f2ds * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_21(_table_index_0, _factor_0)) * ((var_intracellular_sodium_concentration__Na_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_22(_table_index_0, _factor_0));
        double var_sodium_calcium_exchanger__i_NaCa_cyt = (0.999 * 0.0005 * ((CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_27(_table_index_0, _factor_0) * pow(var_intracellular_sodium_concentration__Na_i, 3.0) * 2.0) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_28(_table_index_0, _factor_0) * 2744000.0 * var_intracellular_calcium_concentration__Ca_i))) / ((1.0 + (0.0 * ((var_intracellular_calcium_concentration__Ca_i * 2744000.0) + (2.0 * pow(var_intracellular_sodium_concentration__Na_i, 3.0))))) * (1.0 + (var_intracellular_calcium_concentration__Ca_i * 144.927536232)));
        double var_sodium_calcium_exchanger__i_NaCa_ds = (0.001 * 0.0005 * ((CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_29(_table_index_0, _factor_0) * pow(var_intracellular_sodium_concentration__Na_i, 3.0) * 2.0) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_30(_table_index_0, _factor_0) * 2744000.0 * var_intracellular_calcium_concentration__Ca_ds))) / ((1.0 + (0.0 * ((var_intracellular_calcium_concentration__Ca_ds * 2744000.0) + (2.0 * pow(var_intracellular_sodium_concentration__Na_i, 3.0))))) * (1.0 + (var_intracellular_calcium_concentration__Ca_ds * 144.927536232)));
        double var_L_type_Ca_channel__i_Ca_L_Ca_cyt = (((0.0 * 4.0 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2_gate__f2 * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_11(_table_index_0, _factor_0)) * ((var_intracellular_calcium_concentration__Ca_i * 42.2409057568) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_12(_table_index_0, _factor_0));
        double var_L_type_Ca_channel__i_Ca_L_Ca_ds = (((1.0 * 4.0 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2ds_gate__f2ds * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_17(_table_index_0, _factor_0)) * ((var_intracellular_calcium_concentration__Ca_i * 42.2409057568) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_18(_table_index_0, _factor_0));
        double var_calcium_background_current__i_b_Ca = 0.00025 * (var_membrane__V - (13.3568803298 * log(2.0 / var_intracellular_calcium_concentration__Ca_i)));
        double var_sarcoplasmic_reticulum_calcium_pump__K_2 = var_intracellular_calcium_concentration__Ca_i + (var_intracellular_calcium_concentration__Ca_up * 0.00024) + 0.00012 + 0.0003;
        double var_sarcoplasmic_reticulum_calcium_pump__i_up = ((var_intracellular_calcium_concentration__Ca_i / var_sarcoplasmic_reticulum_calcium_pump__K_2) * 0.4) - (((var_intracellular_calcium_concentration__Ca_up * 0.00024) / var_sarcoplasmic_reticulum_calcium_pump__K_2) * 0.03);
        double var_calcium_translocation__i_trans = 50.0 * (var_intracellular_calcium_concentration__Ca_up - var_intracellular_calcium_concentration__Ca_rel);
        double var_calcium_release__i_rel = ((pow(var_calcium_release__ActFrac / (var_calcium_release__ActFrac + 0.25), 2.0) * 250.0) + 0.05) * var_intracellular_calcium_concentration__Ca_rel;
        double var_calcium_release__CaiReg = var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005);
        double var_calcium_release__RegBindSite = var_calcium_release__CaiReg + ((1.0 - var_calcium_release__CaiReg) * (var_intracellular_calcium_concentration__Ca_ds / (var_intracellular_calcium_concentration__Ca_ds + 0.01)));
        double var_calcium_release__InactRate = 60.0 + (500.0 * pow(var_calcium_release__RegBindSite, 2.0));
        double var_calcium_release__SpeedRel = (var_membrane__V <  -50.0) ? 5.0 : 1.0;
        double d_dt_L_type_Ca_channel_f2_gate__f2 = 1.0 - (1.0 * ((var_intracellular_calcium_concentration__Ca_i / (100000.0 + var_intracellular_calcium_concentration__Ca_i)) + var_L_type_Ca_channel_f2_gate__f2));
        double d_dt_L_type_Ca_channel_f2ds_gate__f2ds = 20.0 * (1.0 - ((var_intracellular_calcium_concentration__Ca_ds / (0.001 + var_intracellular_calcium_concentration__Ca_ds)) + var_L_type_Ca_channel_f2ds_gate__f2ds));
        double d_dt_calcium_release__ActFrac = (((1.0 - var_calcium_release__ActFrac) - var_calcium_release__ProdFrac) * var_calcium_release__SpeedRel * (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_31(_table_index_0, _factor_0) + (500.0 * pow(var_calcium_release__RegBindSite, 2.0)))) - (var_calcium_release__ActFrac * var_calcium_release__SpeedRel * var_calcium_release__InactRate);
        double d_dt_calcium_release__ProdFrac = (var_calcium_release__ActFrac * var_calcium_release__SpeedRel * var_calcium_release__InactRate) - (var_calcium_release__SpeedRel * 1.0 * var_calcium_release__ProdFrac);
        double d_dt_intracellular_sodium_concentration__Na_i =  -0.631827460168 * (var_fast_sodium_current__i_Na + var_persistent_sodium_current__i_p_Na + var_sodium_background_current__i_b_Na + (3.0 * var_sodium_potassium_pump__i_NaK) + (3.0 * var_sodium_calcium_exchanger__i_NaCa_cyt) + var_L_type_Ca_channel__i_Ca_L_Na_cyt + var_L_type_Ca_channel__i_Ca_L_Na_ds);
        double d_dt_intracellular_potassium_concentration__K_i =  -0.631827460168 * ((var_time_independent_potassium_current__i_K1 + var_rapid_delayed_rectifier_potassium_current__i_Kr + var_slow_delayed_rectifier_potassium_current__i_Ks + var_L_type_Ca_channel__i_Ca_L_K_cyt + var_L_type_Ca_channel__i_Ca_L_K_ds + var_transient_outward_current__i_to) - (2.0 * var_sodium_potassium_pump__i_NaK));
        double d_dt_intracellular_calcium_concentration__Ca_Trop = (100000.0 * var_intracellular_calcium_concentration__Ca_i * (0.05 - var_intracellular_calcium_concentration__Ca_Trop)) - (200.0 * var_intracellular_calcium_concentration__Ca_Trop);
        double d_dt_intracellular_calcium_concentration__Ca_Calmod = (100000.0 * var_intracellular_calcium_concentration__Ca_i * (0.02 - var_intracellular_calcium_concentration__Ca_Calmod)) - (50.0 * var_intracellular_calcium_concentration__Ca_Calmod);
        double d_dt_intracellular_calcium_concentration__Ca_i = (((( -0.315913730084 * (((var_L_type_Ca_channel__i_Ca_L_Ca_cyt + var_calcium_background_current__i_b_Ca) - (2.0 * var_sodium_calcium_exchanger__i_NaCa_cyt)) - (2.0 * var_sodium_calcium_exchanger__i_NaCa_ds))) + (var_intracellular_calcium_concentration__Ca_ds * 0.1 * 10.0) + ((var_calcium_release__i_rel * 0.1) * 2.04081632653)) - d_dt_intracellular_calcium_concentration__Ca_Calmod) - d_dt_intracellular_calcium_concentration__Ca_Trop) - var_sarcoplasmic_reticulum_calcium_pump__i_up;
        double d_dt_intracellular_calcium_concentration__Ca_ds = (( -1.0 * var_L_type_Ca_channel__i_Ca_L_Ca_ds) * 3.15913730084) - (var_intracellular_calcium_concentration__Ca_ds * 10.0);
        double d_dt_intracellular_calcium_concentration__Ca_up = (49.0 * var_sarcoplasmic_reticulum_calcium_pump__i_up) - var_calcium_translocation__i_trans;
        double d_dt_intracellular_calcium_concentration__Ca_rel = (0.1 * var_calcium_translocation__i_trans) - var_calcium_release__i_rel;
        
        rResidual[0] = rCurrentGuess[0] - rY[8] - mDt*0.001*d_dt_L_type_Ca_channel_f2_gate__f2;
        rResidual[1] = rCurrentGuess[1] - rY[9] - mDt*0.001*d_dt_L_type_Ca_channel_f2ds_gate__f2ds;
        rResidual[2] = rCurrentGuess[2] - rY[12] - mDt*0.001*d_dt_calcium_release__ActFrac;
        rResidual[3] = rCurrentGuess[3] - rY[13] - mDt*0.001*d_dt_calcium_release__ProdFrac;
        rResidual[11] = rCurrentGuess[11] - rY[14] - mDt*0.001*d_dt_intracellular_sodium_concentration__Na_i;
        rResidual[10] = rCurrentGuess[10] - rY[15] - mDt*0.001*d_dt_intracellular_potassium_concentration__K_i;
        rResidual[7] = rCurrentGuess[7] - rY[16] - mDt*0.001*d_dt_intracellular_calcium_concentration__Ca_i;
        rResidual[6] = rCurrentGuess[6] - rY[17] - mDt*0.001*d_dt_intracellular_calcium_concentration__Ca_ds;
        rResidual[9] = rCurrentGuess[9] - rY[18] - mDt*0.001*d_dt_intracellular_calcium_concentration__Ca_up;
        rResidual[8] = rCurrentGuess[8] - rY[19] - mDt*0.001*d_dt_intracellular_calcium_concentration__Ca_rel;
        rResidual[4] = rCurrentGuess[4] - rY[20] - mDt*0.001*d_dt_intracellular_calcium_concentration__Ca_Calmod;
        rResidual[5] = rCurrentGuess[5] - rY[21] - mDt*0.001*d_dt_intracellular_calcium_concentration__Ca_Trop;
    }

    void ComputeJacobian(const double rCurrentGuess[12], double rJacobian[12][12])
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
        double var_transient_outward_current_s_gate__s = rY[10];
        // Units: dimensionless; Initial value: 0.9948645
        double var_transient_outward_current_r_gate__r = rY[11];
        // Units: dimensionless; Initial value: 0
        
        double var_L_type_Ca_channel_f2_gate__f2 = rCurrentGuess[0];
        double var_L_type_Ca_channel_f2ds_gate__f2ds = rCurrentGuess[1];
        double var_calcium_release__ActFrac = rCurrentGuess[2];
        double var_calcium_release__ProdFrac = rCurrentGuess[3];
        double var_intracellular_calcium_concentration__Ca_Calmod = rCurrentGuess[4];
        double var_intracellular_calcium_concentration__Ca_Trop = rCurrentGuess[5];
        double var_intracellular_calcium_concentration__Ca_ds = rCurrentGuess[6];
        double var_intracellular_calcium_concentration__Ca_i = rCurrentGuess[7];
        double var_intracellular_calcium_concentration__Ca_rel = rCurrentGuess[8];
        double var_intracellular_calcium_concentration__Ca_up = rCurrentGuess[9];
        double var_intracellular_potassium_concentration__K_i = rCurrentGuess[10];
        double var_intracellular_sodium_concentration__Na_i = rCurrentGuess[11];
        
        
        rJacobian[0][0] = 1.0 + mDt;
        rJacobian[0][1] = 0.0;
        rJacobian[0][2] = 0.0;
        rJacobian[0][3] = 0.0;
        rJacobian[0][4] = 0.0;
        rJacobian[0][5] = 0.0;
        rJacobian[0][6] = 0.0;
        rJacobian[0][7] = (-mDt) * (((-1.0) / (100000.0 + var_intracellular_calcium_concentration__Ca_i)) + (var_intracellular_calcium_concentration__Ca_i / pow(100000.0 + var_intracellular_calcium_concentration__Ca_i, 2.0)));
        rJacobian[0][8] = 0.0;
        rJacobian[0][9] = 0.0;
        rJacobian[0][10] = 0.0;
        rJacobian[0][11] = 0.0;
        rJacobian[1][0] = 0.0;
        rJacobian[1][1] = 1.0 + (20.0 * mDt);
        rJacobian[1][2] = 0.0;
        rJacobian[1][3] = 0.0;
        rJacobian[1][4] = 0.0;
        rJacobian[1][5] = 0.0;
        rJacobian[1][6] = (-mDt) * (((-20.0) / (0.001 + var_intracellular_calcium_concentration__Ca_ds)) + ((20.0 * var_intracellular_calcium_concentration__Ca_ds) / pow(0.001 + var_intracellular_calcium_concentration__Ca_ds, 2.0)));
        rJacobian[1][7] = 0.0;
        rJacobian[1][8] = 0.0;
        rJacobian[1][9] = 0.0;
        rJacobian[1][10] = 0.0;
        rJacobian[1][11] = 0.0;
        rJacobian[2][0] = 0.0;
        rJacobian[2][1] = 0.0;
        rJacobian[2][2] = 1.0 - (mDt * ((((-500.0) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * pow((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + 0.01)), 2.0)) - (((var_membrane__V < (-50.0)) ? 5.0 : 1.0) * (60.0 + (500.0 * pow((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + 0.01)), 2.0))))));
        rJacobian[2][3] = ((500.0 * mDt) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * pow((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + 0.01)), 2.0);
        rJacobian[2][4] = 0.0;
        rJacobian[2][5] = 0.0;
        rJacobian[2][6] = (-mDt) * (((((1000.0 * ((1.0 - var_calcium_release__ActFrac) - var_calcium_release__ProdFrac)) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * ((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + 0.01)))) * (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005))) / (var_intracellular_calcium_concentration__Ca_ds + 0.01)) - (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005))) * var_intracellular_calcium_concentration__Ca_ds) / pow(var_intracellular_calcium_concentration__Ca_ds + 0.01, 2.0)))) - ((((1000.0 * var_calcium_release__ActFrac) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * ((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + 0.01)))) * (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005))) / (var_intracellular_calcium_concentration__Ca_ds + 0.01)) - (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005))) * var_intracellular_calcium_concentration__Ca_ds) / pow(var_intracellular_calcium_concentration__Ca_ds + 0.01, 2.0)))));
        rJacobian[2][7] = (-mDt) * (((((1000.0 * ((1.0 - var_calcium_release__ActFrac) - var_calcium_release__ProdFrac)) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * ((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + 0.01)))) * (((1.0 / (var_intracellular_calcium_concentration__Ca_i + 0.0005)) - (var_intracellular_calcium_concentration__Ca_i / pow(var_intracellular_calcium_concentration__Ca_i + 0.0005, 2.0))) + (((((-1.0) / (var_intracellular_calcium_concentration__Ca_i + 0.0005)) + (var_intracellular_calcium_concentration__Ca_i / pow(var_intracellular_calcium_concentration__Ca_i + 0.0005, 2.0))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + 0.01)))) - ((((1000.0 * var_calcium_release__ActFrac) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * ((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + 0.01)))) * (((1.0 / (var_intracellular_calcium_concentration__Ca_i + 0.0005)) - (var_intracellular_calcium_concentration__Ca_i / pow(var_intracellular_calcium_concentration__Ca_i + 0.0005, 2.0))) + (((((-1.0) / (var_intracellular_calcium_concentration__Ca_i + 0.0005)) + (var_intracellular_calcium_concentration__Ca_i / pow(var_intracellular_calcium_concentration__Ca_i + 0.0005, 2.0))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + 0.01)))));
        rJacobian[2][8] = 0.0;
        rJacobian[2][9] = 0.0;
        rJacobian[2][10] = 0.0;
        rJacobian[2][11] = 0.0;
        rJacobian[3][0] = 0.0;
        rJacobian[3][1] = 0.0;
        rJacobian[3][2] = ((-mDt) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * (60.0 + (500.0 * pow((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + 0.01)), 2.0)));
        rJacobian[3][3] = 1.0 + (mDt * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0));
        rJacobian[3][4] = 0.0;
        rJacobian[3][5] = 0.0;
        rJacobian[3][6] = (((((-1000.0) * mDt) * var_calcium_release__ActFrac) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * ((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + 0.01)))) * (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005))) / (var_intracellular_calcium_concentration__Ca_ds + 0.01)) - (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005))) * var_intracellular_calcium_concentration__Ca_ds) / pow(var_intracellular_calcium_concentration__Ca_ds + 0.01, 2.0)));
        rJacobian[3][7] = (((((-1000.0) * mDt) * var_calcium_release__ActFrac) * ((var_membrane__V < (-50.0)) ? 5.0 : 1.0)) * ((var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005)) + (((1.0 - (var_intracellular_calcium_concentration__Ca_i / (var_intracellular_calcium_concentration__Ca_i + 0.0005))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + 0.01)))) * (((1.0 / (var_intracellular_calcium_concentration__Ca_i + 0.0005)) - (var_intracellular_calcium_concentration__Ca_i / pow(var_intracellular_calcium_concentration__Ca_i + 0.0005, 2.0))) + (((((-1.0) / (var_intracellular_calcium_concentration__Ca_i + 0.0005)) + (var_intracellular_calcium_concentration__Ca_i / pow(var_intracellular_calcium_concentration__Ca_i + 0.0005, 2.0))) * var_intracellular_calcium_concentration__Ca_ds) / (var_intracellular_calcium_concentration__Ca_ds + 0.01)));
        rJacobian[3][8] = 0.0;
        rJacobian[3][9] = 0.0;
        rJacobian[3][10] = 0.0;
        rJacobian[3][11] = 0.0;
        rJacobian[4][0] = 0.0;
        rJacobian[4][1] = 0.0;
        rJacobian[4][2] = 0.0;
        rJacobian[4][3] = 0.0;
        rJacobian[4][4] = 1.0 - (mDt * (((-100000.0) * var_intracellular_calcium_concentration__Ca_i) - 50.0));
        rJacobian[4][5] = 0.0;
        rJacobian[4][6] = 0.0;
        rJacobian[4][7] = (-mDt) * (2000.0 - (100000.0 * var_intracellular_calcium_concentration__Ca_Calmod));
        rJacobian[4][8] = 0.0;
        rJacobian[4][9] = 0.0;
        rJacobian[4][10] = 0.0;
        rJacobian[4][11] = 0.0;
        rJacobian[5][0] = 0.0;
        rJacobian[5][1] = 0.0;
        rJacobian[5][2] = 0.0;
        rJacobian[5][3] = 0.0;
        rJacobian[5][4] = 0.0;
        rJacobian[5][5] = 1.0 - (mDt * (((-100000.0) * var_intracellular_calcium_concentration__Ca_i) - 200.0));
        rJacobian[5][6] = 0.0;
        rJacobian[5][7] = (-mDt) * (5000.0 - (100000.0 * var_intracellular_calcium_concentration__Ca_Trop));
        rJacobian[5][8] = 0.0;
        rJacobian[5][9] = 0.0;
        rJacobian[5][10] = 0.0;
        rJacobian[5][11] = 0.0;
        rJacobian[6][0] = 0.0;
        rJacobian[6][1] = (((((0.04730352036 * mDt) * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * (var_membrane__V - 50.0)) / (1.0 - exp(((-0.07486778168) * var_membrane__V) + 3.743389084))) * ((42.24090583 * var_intracellular_calcium_concentration__Ca_i) - (2.0 * exp(((-0.07486778168) * var_membrane__V) + 3.743389084)));
        rJacobian[6][2] = 0.0;
        rJacobian[6][3] = 0.0;
        rJacobian[6][4] = 0.0;
        rJacobian[6][5] = 0.0;
        rJacobian[6][6] = 1.0 + (10.0 * mDt);
        rJacobian[6][7] = (((((1.998143549 * mDt) * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * var_L_type_Ca_channel_f2ds_gate__f2ds) * (var_membrane__V - 50.0)) / (1.0 - exp(((-0.07486778168) * var_membrane__V) + 3.743389084));
        rJacobian[6][8] = 0.0;
        rJacobian[6][9] = 0.0;
        rJacobian[6][10] = 0.0;
        rJacobian[6][11] = 0.0;
        rJacobian[7][0] = 0.0;
        rJacobian[7][1] = 0.0;
        rJacobian[7][2] = (((-0.2040816327) * mDt) * (((500.0 * var_calcium_release__ActFrac) / pow(var_calcium_release__ActFrac + 0.25, 2.0)) - ((500.0 * pow(var_calcium_release__ActFrac, 2.0)) / pow(var_calcium_release__ActFrac + 0.25, 3.0)))) * var_intracellular_calcium_concentration__Ca_rel;
        rJacobian[7][3] = 0.0;
        rJacobian[7][4] = (-mDt) * ((100000.0 * var_intracellular_calcium_concentration__Ca_i) + 50.0);
        rJacobian[7][5] = (-mDt) * ((100000.0 * var_intracellular_calcium_concentration__Ca_i) + 200.0);
        rJacobian[7][6] = (-mDt) * (((((-0.8668672751) * exp((-0.01871694541) * var_membrane__V)) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_ds))) - ((4.578459854e-05 * (((2.0 * exp(0.01871694541 * var_membrane__V)) * pow(var_intracellular_sodium_concentration__Na_i, 3.0)) - ((2744000.0 * exp((-0.01871694541) * var_membrane__V)) * var_intracellular_calcium_concentration__Ca_ds))) / pow(1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_ds), 2.0))) + 1.0);
        rJacobian[7][7] = 1.0 - (mDt * ((((((((((-0.001054905472) / var_intracellular_calcium_concentration__Ca_i) - ((866.0004077 * exp((-0.01871694541) * var_membrane__V)) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_i)))) - ((0.04573881393 * (((2.0 * exp(0.01871694541 * var_membrane__V)) * pow(var_intracellular_sodium_concentration__Na_i, 3.0)) - ((2744000.0 * exp((-0.01871694541) * var_membrane__V)) * var_intracellular_calcium_concentration__Ca_i))) / pow(1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_i), 2.0))) - 7000.0) + (100000.0 * var_intracellular_calcium_concentration__Ca_Calmod)) + (100000.0 * var_intracellular_calcium_concentration__Ca_Trop)) - (0.4 / ((var_intracellular_calcium_concentration__Ca_i + (0.00024 * var_intracellular_calcium_concentration__Ca_up)) + 0.00042))) + ((0.4 * var_intracellular_calcium_concentration__Ca_i) / pow((var_intracellular_calcium_concentration__Ca_i + (0.00024 * var_intracellular_calcium_concentration__Ca_up)) + 0.00042, 2.0))) - ((7.2e-06 * var_intracellular_calcium_concentration__Ca_up) / pow((var_intracellular_calcium_concentration__Ca_i + (0.00024 * var_intracellular_calcium_concentration__Ca_up)) + 0.00042, 2.0))));
        rJacobian[7][8] = (-mDt) * (((51.02040818 * pow(var_calcium_release__ActFrac, 2.0)) / pow(var_calcium_release__ActFrac + 0.25, 2.0)) + 0.01020408164);
        rJacobian[7][9] = (-mDt) * ((((9.6e-05 * var_intracellular_calcium_concentration__Ca_i) / pow((var_intracellular_calcium_concentration__Ca_i + (0.00024 * var_intracellular_calcium_concentration__Ca_up)) + 0.00042, 2.0)) + (7.2e-06 / ((var_intracellular_calcium_concentration__Ca_i + (0.00024 * var_intracellular_calcium_concentration__Ca_up)) + 0.00042))) - ((1.728e-09 * var_intracellular_calcium_concentration__Ca_up) / pow((var_intracellular_calcium_concentration__Ca_i + (0.00024 * var_intracellular_calcium_concentration__Ca_up)) + 0.00042, 2.0)));
        rJacobian[7][10] = 0.0;
        rJacobian[7][11] = (-mDt) * ((((0.001893586897 * exp(0.01871694541 * var_membrane__V)) * pow(var_intracellular_sodium_concentration__Na_i, 2.0)) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_i))) + (((1.89548238e-06 * exp(0.01871694541 * var_membrane__V)) * pow(var_intracellular_sodium_concentration__Na_i, 2.0)) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_ds))));
        rJacobian[8][0] = 0.0;
        rJacobian[8][1] = 0.0;
        rJacobian[8][2] = (mDt * (((500.0 * var_calcium_release__ActFrac) / pow(var_calcium_release__ActFrac + 0.25, 2.0)) - ((500.0 * pow(var_calcium_release__ActFrac, 2.0)) / pow(var_calcium_release__ActFrac + 0.25, 3.0)))) * var_intracellular_calcium_concentration__Ca_rel;
        rJacobian[8][3] = 0.0;
        rJacobian[8][4] = 0.0;
        rJacobian[8][5] = 0.0;
        rJacobian[8][6] = 0.0;
        rJacobian[8][7] = 0.0;
        rJacobian[8][8] = 1.0 - (mDt * ((-5.05) - ((250.0 * pow(var_calcium_release__ActFrac, 2.0)) / pow(var_calcium_release__ActFrac + 0.25, 2.0))));
        rJacobian[8][9] = (-5.0) * mDt;
        rJacobian[8][10] = 0.0;
        rJacobian[8][11] = 0.0;
        rJacobian[9][0] = 0.0;
        rJacobian[9][1] = 0.0;
        rJacobian[9][2] = 0.0;
        rJacobian[9][3] = 0.0;
        rJacobian[9][4] = 0.0;
        rJacobian[9][5] = 0.0;
        rJacobian[9][6] = 0.0;
        rJacobian[9][7] = (-mDt) * (((19.6 / ((var_intracellular_calcium_concentration__Ca_i + (0.00024 * var_intracellular_calcium_concentration__Ca_up)) + 0.00042)) - ((19.6 * var_intracellular_calcium_concentration__Ca_i) / pow((var_intracellular_calcium_concentration__Ca_i + (0.00024 * var_intracellular_calcium_concentration__Ca_up)) + 0.00042, 2.0))) + ((0.0003528 * var_intracellular_calcium_concentration__Ca_up) / pow((var_intracellular_calcium_concentration__Ca_i + (0.00024 * var_intracellular_calcium_concentration__Ca_up)) + 0.00042, 2.0)));
        rJacobian[9][8] = (-50.0) * mDt;
        rJacobian[9][9] = 1.0 - (mDt * ((((((-0.004704) * var_intracellular_calcium_concentration__Ca_i) / pow((var_intracellular_calcium_concentration__Ca_i + (0.00024 * var_intracellular_calcium_concentration__Ca_up)) + 0.00042, 2.0)) - (0.0003528 / ((var_intracellular_calcium_concentration__Ca_i + (0.00024 * var_intracellular_calcium_concentration__Ca_up)) + 0.00042))) + ((8.4672e-08 * var_intracellular_calcium_concentration__Ca_up) / pow((var_intracellular_calcium_concentration__Ca_i + (0.00024 * var_intracellular_calcium_concentration__Ca_up)) + 0.00042, 2.0))) - 50.0));
        rJacobian[9][10] = 0.0;
        rJacobian[9][11] = 0.0;
        rJacobian[10][0] = 0.0;
        rJacobian[10][1] = (((((4.730352032e-06 * mDt) * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * (var_membrane__V - 50.0)) / (1.0 - exp(((-0.03743389084) * var_membrane__V) + 1.871694542))) * ((6.499300411 * var_intracellular_potassium_concentration__K_i) - (4.0 * exp(((-0.03743389084) * var_membrane__V) + 1.871694542)));
        rJacobian[10][2] = 0.0;
        rJacobian[10][3] = 0.0;
        rJacobian[10][4] = 0.0;
        rJacobian[10][5] = 0.0;
        rJacobian[10][6] = 0.0;
        rJacobian[10][7] = 0.0;
        rJacobian[10][8] = 0.0;
        rJacobian[10][9] = 0.0;
        rJacobian[10][10] = 1.0 - (mDt * ((((((((-2.411212507) / var_intracellular_potassium_concentration__K_i) / (1.0 + exp(((0.04679236355 * var_membrane__V) - (1.25 * log(4.0 / var_intracellular_potassium_concentration__K_i))) - 0.4679236355))) + ((((0.1128263322 * (var_membrane__V - (26.71376066 * log(4.0 / var_intracellular_potassium_concentration__K_i)))) / pow(1.0 + exp(((0.04679236355 * var_membrane__V) - (1.25 * log(4.0 / var_intracellular_potassium_concentration__K_i))) - 0.4679236355), 2.0)) / var_intracellular_potassium_concentration__K_i) * exp(((0.04679236355 * var_membrane__V) - (1.25 * log(4.0 / var_intracellular_potassium_concentration__K_i))) - 0.4679236355))) - (((16.87848754 * ((0.0021 * var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1) + (0.0013 * var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2))) / (1.0 + exp((0.04464285714 * var_membrane__V) + 0.4017857143))) / var_intracellular_potassium_concentration__K_i)) - ((0.04388406762 * pow(var_slow_delayed_rectifier_potassium_current_xs_gate__xs, 2.0)) / (var_intracellular_potassium_concentration__K_i + (0.03 * var_intracellular_sodium_concentration__Na_i)))) - (((((3.074397891e-05 * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * var_L_type_Ca_channel_f2ds_gate__f2ds) * (var_membrane__V - 50.0)) / (1.0 - exp(((-0.03743389084) * var_membrane__V) + 1.871694542)))) - (((0.08439243772 * var_transient_outward_current_s_gate__s) * var_transient_outward_current_r_gate__r) / var_intracellular_potassium_concentration__K_i)));
        rJacobian[10][11] = (-mDt) * (((((-0.001316522028) * pow(var_slow_delayed_rectifier_potassium_current_xs_gate__xs, 2.0)) / (var_intracellular_potassium_concentration__K_i + (0.03 * var_intracellular_sodium_concentration__Na_i))) + (0.7076467552 / (40.0 + var_intracellular_sodium_concentration__Na_i))) - ((0.7076467552 * var_intracellular_sodium_concentration__Na_i) / pow(40.0 + var_intracellular_sodium_concentration__Na_i, 2.0)));
        rJacobian[11][0] = 0.0;
        rJacobian[11][1] = (((((2.365176017e-05 * mDt) * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * (var_membrane__V - 50.0)) / (1.0 - exp(((-0.03743389084) * var_membrane__V) + 1.871694542))) * ((6.499300411 * var_intracellular_sodium_concentration__Na_i) - (140.0 * exp(((-0.03743389084) * var_membrane__V) + 1.871694542)));
        rJacobian[11][2] = 0.0;
        rJacobian[11][3] = 0.0;
        rJacobian[11][4] = 0.0;
        rJacobian[11][5] = 0.0;
        rJacobian[11][6] = 0.0;
        rJacobian[11][7] = (-mDt) * (((2598.001224 * exp((-0.01871694541) * var_membrane__V)) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_i))) + ((0.1372164418 * (((2.0 * exp(0.01871694541 * var_membrane__V)) * pow(var_intracellular_sodium_concentration__Na_i, 3.0)) - ((2744000.0 * exp((-0.01871694541) * var_membrane__V)) * var_intracellular_calcium_concentration__Ca_i))) / pow(1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_i), 2.0)));
        rJacobian[11][8] = 0.0;
        rJacobian[11][9] = 0.0;
        rJacobian[11][10] = (((5.063546263 * mDt) * pow(var_fast_sodium_current_m_gate__m, 3.0)) * var_fast_sodium_current_h_gate__h) / (var_intracellular_sodium_concentration__Na_i + (0.12 * var_intracellular_potassium_concentration__K_i));
        rJacobian[11][11] = 1.0 - (mDt * ((((((((((-42.19621886) * pow(var_fast_sodium_current_m_gate__m, 3.0)) * var_fast_sodium_current_h_gate__h) / (var_intracellular_sodium_concentration__Na_i + (0.12 * var_intracellular_potassium_concentration__K_i))) - ((0.06751395018 / (1.0 + exp((((-1.0) / 8.0) * var_membrane__V) - (13.0 / 2.0)))) / var_intracellular_sodium_concentration__Na_i)) - (0.01012709253 / var_intracellular_sodium_concentration__Na_i)) - (1.061470133 / (40.0 + var_intracellular_sodium_concentration__Na_i))) + ((1.061470133 * var_intracellular_sodium_concentration__Na_i) / pow(40.0 + var_intracellular_sodium_concentration__Na_i, 2.0))) - (((0.005680760693 * exp(0.01871694541 * var_membrane__V)) * pow(var_intracellular_sodium_concentration__Na_i, 2.0)) / (1.0 + (144.9275362 * var_intracellular_calcium_concentration__Ca_i)))) - (((((0.0001537198946 * var_L_type_Ca_channel_d_gate__d) * var_L_type_Ca_channel_f_gate__f) * var_L_type_Ca_channel_f2ds_gate__f2ds) * (var_membrane__V - 50.0)) / (1.0 - exp(((-0.03743389084) * var_membrane__V) + 1.871694542)))));
    }

protected:
    void UpdateTransmembranePotential(double var_environment__time)
    {
        // Time units: second
        var_environment__time *= 0.001;
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
        _table_index_0 = (unsigned) floor(_offset_0_over_table_step);
        _factor_0 = _offset_0_over_table_step - _table_index_0;
        
        double var_reversal_potentials__E_K = 26.7137606597 * log(4.0 / var_intracellular_potassium_concentration__K_i);
        double var_time_independent_potassium_current__i_K1 = (0.142857142857 * (var_membrane__V - var_reversal_potentials__E_K)) / (1.0 + exp((((var_membrane__V - var_reversal_potentials__E_K) - 10.0) * 96485.3415 * 1.25) * 3.87974901066e-07));
        var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
        double var_transient_outward_current__i_to = 0.005 * (0.0 + (var_transient_outward_current_s_gate__s * 1.0)) * var_transient_outward_current_r_gate__r * (var_membrane__V - var_reversal_potentials__E_K);
        var_membrane__i_to = var_transient_outward_current__i_to;
        double var_rapid_delayed_rectifier_potassium_current__i_Kr = ((((0.0021 * var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1) + (0.0013 * var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2)) * 1.0) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_0(_table_index_0, _factor_0)) * (var_membrane__V - var_reversal_potentials__E_K);
        var_membrane__i_Kr = var_rapid_delayed_rectifier_potassium_current__i_Kr;
        double var_slow_delayed_rectifier_potassium_current__i_Ks = 0.0026 * pow(var_slow_delayed_rectifier_potassium_current_xs_gate__xs, 2.0) * (var_membrane__V - (26.7137606597 * log(8.2 / (var_intracellular_potassium_concentration__K_i + (0.03 * var_intracellular_sodium_concentration__Na_i)))));
        var_membrane__i_Ks = var_slow_delayed_rectifier_potassium_current__i_Ks;
        double var_L_type_Ca_channel__i_Ca_L_K_cyt = (((0.0 * 0.002 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2_gate__f2 * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_13(_table_index_0, _factor_0)) * ((var_intracellular_potassium_concentration__K_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_14(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_K_cyt = var_L_type_Ca_channel__i_Ca_L_K_cyt;
        double var_L_type_Ca_channel__i_Ca_L_K_ds = (((1.0 * 0.002 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2ds_gate__f2ds * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_19(_table_index_0, _factor_0)) * ((var_intracellular_potassium_concentration__K_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_20(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_K_ds = var_L_type_Ca_channel__i_Ca_L_K_ds;
        double var_sodium_potassium_pump__i_NaK = (0.56 * var_intracellular_sodium_concentration__Na_i) / (40.0 + var_intracellular_sodium_concentration__Na_i);
        var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK;
        double var_fast_sodium_current__i_Na = 2.5 * pow(var_fast_sodium_current_m_gate__m, 3.0) * var_fast_sodium_current_h_gate__h * (var_membrane__V - (26.7137606597 * log(140.48 / (var_intracellular_sodium_concentration__Na_i + (0.12 * var_intracellular_potassium_concentration__K_i)))));
        var_membrane__i_Na = var_fast_sodium_current__i_Na;
        double var_reversal_potentials__E_Na = 26.7137606597 * log(140.0 / var_intracellular_sodium_concentration__Na_i);
        double var_sodium_background_current__i_b_Na = 0.0006 * (var_membrane__V - var_reversal_potentials__E_Na);
        var_membrane__i_b_Na = var_sodium_background_current__i_b_Na;
        double var_persistent_sodium_current__i_p_Na = CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_10(_table_index_0, _factor_0) * (var_membrane__V - var_reversal_potentials__E_Na);
        var_membrane__i_p_Na = var_persistent_sodium_current__i_p_Na;
        double var_L_type_Ca_channel__i_Ca_L_Na_cyt = (((0.0 * 0.01 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2_gate__f2 * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_15(_table_index_0, _factor_0)) * ((var_intracellular_sodium_concentration__Na_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_16(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_Na_cyt = var_L_type_Ca_channel__i_Ca_L_Na_cyt;
        double var_L_type_Ca_channel__i_Ca_L_Na_ds = (((1.0 * 0.01 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2ds_gate__f2ds * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_21(_table_index_0, _factor_0)) * ((var_intracellular_sodium_concentration__Na_i * 6.49930040518) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_22(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_Na_ds = var_L_type_Ca_channel__i_Ca_L_Na_ds;
        double var_sodium_calcium_exchanger__i_NaCa_cyt = (0.999 * 0.0005 * ((CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_27(_table_index_0, _factor_0) * pow(var_intracellular_sodium_concentration__Na_i, 3.0) * 2.0) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_28(_table_index_0, _factor_0) * 2744000.0 * var_intracellular_calcium_concentration__Ca_i))) / ((1.0 + (0.0 * ((var_intracellular_calcium_concentration__Ca_i * 2744000.0) + (2.0 * pow(var_intracellular_sodium_concentration__Na_i, 3.0))))) * (1.0 + (var_intracellular_calcium_concentration__Ca_i * 144.927536232)));
        var_membrane__i_NaCa_cyt = var_sodium_calcium_exchanger__i_NaCa_cyt;
        double var_sodium_calcium_exchanger__i_NaCa_ds = (0.001 * 0.0005 * ((CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_29(_table_index_0, _factor_0) * pow(var_intracellular_sodium_concentration__Na_i, 3.0) * 2.0) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_30(_table_index_0, _factor_0) * 2744000.0 * var_intracellular_calcium_concentration__Ca_ds))) / ((1.0 + (0.0 * ((var_intracellular_calcium_concentration__Ca_ds * 2744000.0) + (2.0 * pow(var_intracellular_sodium_concentration__Na_i, 3.0))))) * (1.0 + (var_intracellular_calcium_concentration__Ca_ds * 144.927536232)));
        var_membrane__i_NaCa_ds = var_sodium_calcium_exchanger__i_NaCa_ds;
        double var_L_type_Ca_channel__i_Ca_L_Ca_cyt = (((0.0 * 4.0 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2_gate__f2 * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_11(_table_index_0, _factor_0)) * ((var_intracellular_calcium_concentration__Ca_i * 42.2409057568) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_12(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_Ca_cyt = var_L_type_Ca_channel__i_Ca_L_Ca_cyt;
        double var_L_type_Ca_channel__i_Ca_L_Ca_ds = (((1.0 * 4.0 * 0.1 * var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f2ds_gate__f2ds * (var_membrane__V - 50.0) * 96485.3415) * 3.87974901066e-07) / CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_17(_table_index_0, _factor_0)) * ((var_intracellular_calcium_concentration__Ca_i * 42.2409057568) - CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_18(_table_index_0, _factor_0));
        var_membrane__i_Ca_L_Ca_ds = var_L_type_Ca_channel__i_Ca_L_Ca_ds;
        double var_calcium_background_current__i_b_Ca = 0.00025 * (var_membrane__V - (13.3568803298 * log(2.0 / var_intracellular_calcium_concentration__Ca_i)));
        var_membrane__i_b_Ca = var_calcium_background_current__i_b_Ca;
        double var_membrane__i_Stim = GetStimulus((1.0/0.001)*var_environment__time);
        double d_dt_membrane__V =  -10526.3157895 * (var_membrane__i_Stim + var_membrane__i_K1 + var_membrane__i_to + var_membrane__i_Kr + var_membrane__i_Ks + var_membrane__i_NaK + var_membrane__i_Na + var_membrane__i_b_Na + var_membrane__i_p_Na + var_membrane__i_Ca_L_Na_cyt + var_membrane__i_Ca_L_Na_ds + var_membrane__i_NaCa_cyt + var_membrane__i_NaCa_ds + var_membrane__i_Ca_L_Ca_cyt + var_membrane__i_Ca_L_Ca_ds + var_membrane__i_Ca_L_K_cyt + var_membrane__i_Ca_L_K_ds + var_membrane__i_b_Ca);
        
        rY[0] += mDt * 0.001*d_dt_membrane__V;
    }

    void ComputeOneStepExceptVoltage(double var_environment__time)
    {
        // Time units: second
        var_environment__time *= 0.001;
        std::vector<double>& rY = rGetStateVariables();
        double var_membrane__V = rY[0];
        // Units: millivolt; Initial value: -92.849333
        
        // Lookup table indexing
        if (var_membrane__V>49.9999 || var_membrane__V<-100.0001)
            EXCEPTION(DumpState("V outside lookup table range", rY));
        double _offset_0 = var_membrane__V - -100.0001;
        double _offset_0_over_table_step = _offset_0 * 100.0;
        _table_index_0 = (unsigned) floor(_offset_0_over_table_step);
        _factor_0 = _offset_0_over_table_step - _table_index_0;
        
        double var_fast_sodium_current_m_gate__E0_m = var_membrane__V + 41.0;
        double var_L_type_Ca_channel_d_gate__E0_d = (var_membrane__V + 24.0) - 5.0;
        double var_L_type_Ca_channel_f_gate__E0_f = var_membrane__V + 34.0;
        
        const double _g_0 = 3.0 * (((fabs(var_L_type_Ca_channel_d_gate__E0_d) < 0.0001) ? 120.0 : ((30.0 * var_L_type_Ca_channel_d_gate__E0_d) / (1.0 - exp((-var_L_type_Ca_channel_d_gate__E0_d) * 0.25)))) * 1.0);
        const double _h_0 = 3.0 * ((((fabs(var_L_type_Ca_channel_d_gate__E0_d) < 0.0001) ? 120.0 : ((30.0 * var_L_type_Ca_channel_d_gate__E0_d) / (1.0 - exp((-var_L_type_Ca_channel_d_gate__E0_d) * 0.25)))) * (-1.0)) - (((fabs(var_L_type_Ca_channel_d_gate__E0_d) < 0.0001) ? 120.0 : ((12.0 * var_L_type_Ca_channel_d_gate__E0_d) / (exp(var_L_type_Ca_channel_d_gate__E0_d * 0.1) - 1.0))) * 1.0));
        const double _g_1 = 0.3 * (((fabs(var_L_type_Ca_channel_f_gate__E0_f) < 0.0001) ? 25.0 : ((6.25 * var_L_type_Ca_channel_f_gate__E0_f) / (exp(var_L_type_Ca_channel_f_gate__E0_f * 0.25) - 1.0))) * 1.0);
        const double _h_1 = 0.3 * ((((fabs(var_L_type_Ca_channel_f_gate__E0_f) < 0.0001) ? 25.0 : ((6.25 * var_L_type_Ca_channel_f_gate__E0_f) / (exp(var_L_type_Ca_channel_f_gate__E0_f * 0.25) - 1.0))) * (-1.0)) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_23(_table_index_0, _factor_0) * 1.0));
        const double _g_2 = CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_8(_table_index_0, _factor_0) * 1.0;
        const double _h_2 = (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_8(_table_index_0, _factor_0) * (-1.0)) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_9(_table_index_0, _factor_0) * 1.0);
        const double _g_3 = ((fabs(var_fast_sodium_current_m_gate__E0_m) < 1e-05) ? 2000.0 : ((200.0 * var_fast_sodium_current_m_gate__E0_m) / (1.0 - exp( -0.1 * var_fast_sodium_current_m_gate__E0_m)))) * 1.0;
        const double _h_3 = (((fabs(var_fast_sodium_current_m_gate__E0_m) < 1e-05) ? 2000.0 : ((200.0 * var_fast_sodium_current_m_gate__E0_m) / (1.0 - exp( -0.1 * var_fast_sodium_current_m_gate__E0_m)))) * (-1.0)) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_7(_table_index_0, _factor_0) * 1.0);
        const double _g_4 = CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_1(_table_index_0, _factor_0) * 1.0;
        const double _h_4 = (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_1(_table_index_0, _factor_0) * (-1.0)) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_2(_table_index_0, _factor_0) * 1.0);
        const double _g_5 = CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_3(_table_index_0, _factor_0) * 1.0;
        const double _h_5 = (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_3(_table_index_0, _factor_0) * (-1.0)) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_4(_table_index_0, _factor_0) * 1.0);
        const double _g_6 = CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_5(_table_index_0, _factor_0) * 1.0;
        const double _h_6 = (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_5(_table_index_0, _factor_0) * (-1.0)) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_6(_table_index_0, _factor_0) * 1.0);
        const double _g_7 = 333.0 * CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_26(_table_index_0, _factor_0);
        const double _h_7 = 333.0 * (-1.0);
        const double _g_8 = CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_24(_table_index_0, _factor_0) * 1.0;
        const double _h_8 = (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_24(_table_index_0, _factor_0) * (-1.0)) - (CML_noble_varghese_kohl_noble_1998_basic_pe_lut_be_LookupTables::Instance()->_lookup_25(_table_index_0, _factor_0) * 1.0);
        
        const double dt = mDt*0.001;
        rY[1] = (rY[1] + _g_4*dt) / (1 - _h_4*dt);
        rY[2] = (rY[2] + _g_5*dt) / (1 - _h_5*dt);
        rY[3] = (rY[3] + _g_6*dt) / (1 - _h_6*dt);
        rY[4] = (rY[4] + _g_3*dt) / (1 - _h_3*dt);
        rY[5] = (rY[5] + _g_2*dt) / (1 - _h_2*dt);
        rY[6] = (rY[6] + _g_0*dt) / (1 - _h_0*dt);
        rY[7] = (rY[7] + _g_1*dt) / (1 - _h_1*dt);
        rY[10] = (rY[10] + _g_8*dt) / (1 - _h_8*dt);
        rY[11] = (rY[11] + _g_7*dt) / (1 - _h_7*dt);
        
        double _guess[12] = {rY[8],rY[9],rY[12],rY[13],rY[14],rY[15],rY[16],rY[17],rY[18],rY[19],rY[20],rY[21]};
        CardiacNewtonSolver<12> *_solver = CardiacNewtonSolver<12>::Instance();
        _solver->Solve(*this, _guess);
        rY[8] = _guess[0];
        rY[9] = _guess[1];
        rY[12] = _guess[2];
        rY[13] = _guess[3];
        rY[20] = _guess[4];
        rY[21] = _guess[5];
        rY[17] = _guess[6];
        rY[16] = _guess[7];
        rY[19] = _guess[8];
        rY[18] = _guess[9];
        rY[15] = _guess[10];
        rY[14] = _guess[11];
    }

};

#endif
