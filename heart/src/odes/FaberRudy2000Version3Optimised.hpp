#ifndef _FaberRudy2000Version3Optimised
#define _FaberRudy2000Version3Optimised

// Model: LR_Dynamic_model_2000
// Processed by pycml - CellML Tools in Python
//     (translate: 1004, pycml: 896)
// on Wed Dec 19 11:02:41 2007

#include <cmath>
#include <cassert>
#include "AbstractCardiacCell.hpp"
#include "Exception.hpp"
#include "AbstractStimulusFunction.hpp"

class FaberRudy2000Version3OptimisedLookupTables
{
public:
    static FaberRudy2000Version3OptimisedLookupTables* Instance()
    {
        if (mpInstance == NULL)
        {
            mpInstance = new FaberRudy2000Version3OptimisedLookupTables;
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

    inline double _lookup_32(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][32];
        double y2 = _lookup_table_0[i+1][32];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_33(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][33];
        double y2 = _lookup_table_0[i+1][33];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_34(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][34];
        double y2 = _lookup_table_0[i+1][34];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_35(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][35];
        double y2 = _lookup_table_0[i+1][35];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_36(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][36];
        double y2 = _lookup_table_0[i+1][36];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_37(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][37];
        double y2 = _lookup_table_0[i+1][37];
        return y1 + (y2-y1)*factor;
    }

    inline double _lookup_38(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][38];
        double y2 = _lookup_table_0[i+1][38];
        return y1 + (y2-y1)*factor;
    }


protected:
    FaberRudy2000Version3OptimisedLookupTables(const FaberRudy2000Version3OptimisedLookupTables&);
    FaberRudy2000Version3OptimisedLookupTables& operator= (const FaberRudy2000Version3OptimisedLookupTables&);
    FaberRudy2000Version3OptimisedLookupTables()
    {
        assert(mpInstance == NULL);
        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][0] = 80.0 * exp((-var_membrane__V) * 0.0909090909091);
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][1] = (var_membrane__V <  -40.0) ? (135.0 * exp((80.0 + var_membrane__V) *  -0.147058823529)) : 0.0;
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][2] = (var_membrane__V <  -40.0) ? ((3560.0 * exp(0.079 * var_membrane__V)) + (310000000.0 * exp(0.35 * var_membrane__V))) : (1000.0 / (0.13 * (1.0 + exp((var_membrane__V + 10.66) *  -0.0900900900901))));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][3] = (var_membrane__V <  -40.0) ? ((1000.0 * (-((127140.0 * exp(0.2444 * var_membrane__V)) + (3.474e-05 * exp( -0.04391 * var_membrane__V)))) * (var_membrane__V + 37.78)) / (1.0 + exp(0.311 * (var_membrane__V + 79.23)))) : 0.0;
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][4] = (var_membrane__V <  -40.0) ? ((121.2 * exp( -0.01052 * var_membrane__V)) / (1.0 + exp( -0.1378 * (var_membrane__V + 40.14)))) : ((300.0 * exp( -2.535e-07 * var_membrane__V)) / (1.0 + exp( -0.1 * (var_membrane__V + 32.0))));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][5] = exp((2.0 * var_membrane__V * 96485.0) * 3.87996927064e-07);
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][6] = exp((2.0 * var_membrane__V * 96485.0) * 3.87996927064e-07) - 1.0;
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][7] = exp((1.0 * var_membrane__V * 96485.0) * 3.87996927064e-07);
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][8] = exp((1.0 * var_membrane__V * 96485.0) * 3.87996927064e-07) - 1.0;
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][9] = exp((1.0 * var_membrane__V * 96485.0) * 3.87996927064e-07);
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][10] = exp((1.0 * var_membrane__V * 96485.0) * 3.87996927064e-07) - 1.0;
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][11] = (1.0 / (1.0 + exp((var_membrane__V + 32.0) * 0.125))) + (0.6 / (1.0 + exp((50.0 - var_membrane__V) * 0.05)));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][12] = 0.001 / ((0.0197 * exp(-pow(0.0337 * (var_membrane__V + 10.0), 2.0))) + 0.02);
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][13] = 1.0 / (1.0 + exp((-(var_membrane__V + 14.0)) * 0.0925925925926));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][14] = 0.0037 + (0.0061 / (1.0 + exp((var_membrane__V + 25.0) * 0.222222222222)));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][15] = 1.0 / (1.0 + exp((var_membrane__V + 60.0) * 0.178571428571));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][16] = 1.0 / (1.0 + exp((var_membrane__V + 9.0) * 0.0446428571429));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][17] = 1.0 / (1.0 + exp((-(var_membrane__V + 21.5)) * 0.133333333333));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][18] = 0.001 / (((0.00138 * (var_membrane__V + 14.2)) / (1.0 - exp( -0.123 * (var_membrane__V + 14.2)))) + ((0.00061 * (var_membrane__V + 38.9)) / (exp(0.145 * (var_membrane__V + 38.9)) - 1.0)));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][19] = 1.0 / (1.0 + exp((-(var_membrane__V - 1.5)) * 0.059880239521));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][20] = 0.001 / (((7.19e-05 * (var_membrane__V + 30.0)) / (1.0 - exp( -0.148 * (var_membrane__V + 30.0)))) + ((0.000131 * (var_membrane__V + 30.0)) / (exp(0.0687 * (var_membrane__V + 30.0)) - 1.0)));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][21] = 1.0 / (1.0 + exp((-(var_membrane__V - 1.5)) * 0.059880239521));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][22] = 0.004 / (((7.19e-05 * (var_membrane__V + 30.0)) / (1.0 - exp( -0.148 * (var_membrane__V + 30.0)))) + ((0.000131 * (var_membrane__V + 30.0)) / (exp(0.0687 * (var_membrane__V + 30.0)) - 1.0)));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][23] = 1.0 / (1.0 + exp((7.488 - var_membrane__V) * 0.167224080268));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][24] = 0.8 - (0.65 / (1.0 + exp((var_membrane__V + 125.0) * 0.0666666666667)));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][25] = exp(var_membrane__V * 0.01);
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][26] = (10000.0 * exp((var_membrane__V - 40.0) * 0.04)) / (1.0 + exp((var_membrane__V - 40.0) * 0.04));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][27] = (10000.0 * exp((-(var_membrane__V + 90.0)) * 0.04)) / (1.0 + exp((-(var_membrane__V + 90.0)) * 0.04));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][28] = 15.0 / (1.0 + exp((var_membrane__V + 60.0) * 0.2));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][29] = (100.0 * exp((var_membrane__V + 25.0) * 0.2)) / (1.0 + exp((var_membrane__V + 25.0) * 0.2));
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][30] = 2.25 * (1.0 / (1.0 + (0.1245 * exp(( -0.1 * var_membrane__V * 96485.0) * 3.87996927064e-07)) + (0.0365 * 0.872719796652 * exp(((-var_membrane__V) * 96485.0) * 3.87996927064e-07)))) * 1.0;
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][31] = exp((1.0 * var_membrane__V * 96485.0) * 3.87996927064e-07);
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][32] = exp((1.0 * var_membrane__V * 96485.0) * 3.87996927064e-07) - 1.0;
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][33] = exp((1.0 * var_membrane__V * 96485.0) * 3.87996927064e-07);
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][34] = exp((1.0 * var_membrane__V * 96485.0) * 3.87996927064e-07) - 1.0;
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][35] = exp(( -0.85 * var_membrane__V * 96485.0) * 3.87996927064e-07);
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][36] = exp((var_membrane__V * 96485.0) * 3.87996927064e-07);
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][37] = exp(( -0.85 * var_membrane__V * 96485.0) * 3.87996927064e-07);
        }

        for (int i=0; i<20001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][38] = exp((var_membrane__V * 96485.0) * 3.87996927064e-07);
        }

    }
private:
    /** The single instance of the class */
    static FaberRudy2000Version3OptimisedLookupTables *mpInstance;
    // Lookup tables
    double _lookup_table_0[20001][39];

};

FaberRudy2000Version3OptimisedLookupTables* FaberRudy2000Version3OptimisedLookupTables::mpInstance = NULL;

class FaberRudy2000Version3Optimised : public AbstractCardiacCell
{
public:
    FaberRudy2000Version3Optimised(AbstractIvpOdeSolver *pSolver, double dt,
                                     AbstractStimulusFunction *pIntracellularStimulus,
                                     AbstractStimulusFunction *pExtracellularStimulus=NULL)
        : AbstractCardiacCell(pSolver, 25, 0, dt, pIntracellularStimulus, pExtracellularStimulus)
    {
        // Time units: second

        mVariableNames.push_back("V");
        mVariableUnits.push_back("millivolt");
        mInitialConditions.push_back(-90);

        mVariableNames.push_back("m");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.0008);

        mVariableNames.push_back("h");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.993771);

        mVariableNames.push_back("j");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.995727);

        mVariableNames.push_back("d");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(3.210618e-6);

        mVariableNames.push_back("f");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.999837);

        mVariableNames.push_back("b");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.000970231);

        mVariableNames.push_back("g");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.994305);

        mVariableNames.push_back("xr");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.000124042);

        mVariableNames.push_back("xs1");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.00445683);

        mVariableNames.push_back("xs2");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.00445683);

        mVariableNames.push_back("zdv");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.5);

        mVariableNames.push_back("ydv");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0.5);

        mVariableNames.push_back("CaI");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(6e-5);

        mVariableNames.push_back("Ca_JSR");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(1.8);

        mVariableNames.push_back("Ca_NSR");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(1.8);

        mVariableNames.push_back("APtrack");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0);

        mVariableNames.push_back("APtrack2");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0);

        mVariableNames.push_back("APtrack3");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0);

        mVariableNames.push_back("Cainfluxtrack");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0);

        mVariableNames.push_back("OVRLDtrack");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0);

        mVariableNames.push_back("OVRLDtrack2");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0);

        mVariableNames.push_back("OVRLDtrack3");
        mVariableUnits.push_back("dimensionless");
        mInitialConditions.push_back(0);

        mVariableNames.push_back("Nai");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(9);

        mVariableNames.push_back("Ki");
        mVariableUnits.push_back("millimolar");
        mInitialConditions.push_back(141.2);

        Init();

    }

    ~FaberRudy2000Version3Optimised(void)
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
        // Units: millivolt; Initial value: -90
        double var_fast_sodium_current_m_gate__m = rY[1];
        // Units: dimensionless; Initial value: 0.0008
        double var_fast_sodium_current_h_gate__h = rY[2];
        // Units: dimensionless; Initial value: 0.993771
        double var_fast_sodium_current_j_gate__j = rY[3];
        // Units: dimensionless; Initial value: 0.995727
        double var_L_type_Ca_channel_d_gate__d = rY[4];
        // Units: dimensionless; Initial value: 3.210618e-6
        double var_L_type_Ca_channel_f_gate__f = rY[5];
        // Units: dimensionless; Initial value: 0.999837
        double var_T_type_Ca_channel_b_gate__b = rY[6];
        // Units: dimensionless; Initial value: 0.000970231
        double var_T_type_Ca_channel_g_gate__g = rY[7];
        // Units: dimensionless; Initial value: 0.994305
        double var_rapid_delayed_rectifier_potassium_current_xr_gate__xr = rY[8];
        // Units: dimensionless; Initial value: 0.000124042
        double var_slow_delayed_rectifier_potassium_current_xs1_gate__xs1 = rY[9];
        // Units: dimensionless; Initial value: 0.00445683
        double var_slow_delayed_rectifier_potassium_current_xs2_gate__xs2 = rY[10];
        // Units: dimensionless; Initial value: 0.00445683
        double var_transient_outward_current_zdv_gate__zdv = rY[11];
        // Units: dimensionless; Initial value: 0.5
        double var_transient_outward_current_ydv_gate__ydv = rY[12];
        // Units: dimensionless; Initial value: 0.5
        double var_calcium_dynamics__Cai = rY[13];
        // Units: millimolar; Initial value: 6e-5
        double var_ionic_concentrations__Nai = rY[23];
        // Units: millimolar; Initial value: 9
        double var_ionic_concentrations__Ki = rY[24];
        // Units: millimolar; Initial value: 141.2

        // Lookup table indexing
#define COVERAGE_IGNORE
        if (var_membrane__V>99.9999 || var_membrane__V<-100.0001)
            EXCEPTION(DumpState("V outside lookup table range"));
#undef COVERAGE_IGNORE
        double _offset_0 = var_membrane__V - -100.0001;
        double _offset_0_over_table_step = _offset_0 * 100.0;
        unsigned _table_index_0 = (unsigned) floor(_offset_0_over_table_step);
        double _factor_0 = _offset_0_over_table_step - _table_index_0;

        double var_fast_sodium_current__E_Na = 26.7123387055 * log(132.0 / var_ionic_concentrations__Nai);
        double var_fast_sodium_current__i_Na = 16.0 * pow(var_fast_sodium_current_m_gate__m, 3.0) * var_fast_sodium_current_h_gate__h * var_fast_sodium_current_j_gate__j * (var_membrane__V - var_fast_sodium_current__E_Na);
        double var_membrane__i_Na = var_fast_sodium_current__i_Na;
        double var_L_type_Ca_channel_f_Ca_gate__f_Ca = 1.0 / (1.0 + (var_calcium_dynamics__Cai * 1666.66666667));
        double var_L_type_Ca_channel__i_CaCa = var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f_Ca_gate__f_Ca * ((((0.00054 * 4.0 * var_membrane__V * 9309355225.0) * 3.87996927064e-07) * ((1.0 * var_calcium_dynamics__Cai * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_5(_table_index_0, _factor_0)) - 0.6138)) / FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_6(_table_index_0, _factor_0));
        double var_L_type_Ca_channel__i_CaNa = var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f_Ca_gate__f_Ca * ((((6.75e-07 * 1.0 * var_membrane__V * 9309355225.0) * 3.87996927064e-07) * ((0.75 * var_ionic_concentrations__Nai * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_7(_table_index_0, _factor_0)) - 99.0)) / FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_8(_table_index_0, _factor_0));
        double var_L_type_Ca_channel__i_CaK = var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f_Ca_gate__f_Ca * ((((1.93e-07 * 1.0 * var_membrane__V * 9309355225.0) * 3.87996927064e-07) * ((0.75 * var_ionic_concentrations__Ki * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_9(_table_index_0, _factor_0)) - 3.375)) / FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_10(_table_index_0, _factor_0));
        double var_L_type_Ca_channel__i_Ca_L = var_L_type_Ca_channel__i_CaCa + var_L_type_Ca_channel__i_CaK + var_L_type_Ca_channel__i_CaNa;
        double var_membrane__i_Ca_L = var_L_type_Ca_channel__i_Ca_L;
        double var_calcium_background_current__E_Ca = 13.3561693527 * log(1.8 / var_calcium_dynamics__Cai);
        double var_T_type_Ca_channel__i_Ca_T = 0.05 * var_T_type_Ca_channel_b_gate__b * var_T_type_Ca_channel_b_gate__b * var_T_type_Ca_channel_g_gate__g * (var_membrane__V - var_calcium_background_current__E_Ca);
        double var_membrane__i_Ca_T = var_T_type_Ca_channel__i_Ca_T;
        double var_time_independent_potassium_current__E_K = 26.7123387055 * log(4.5 / var_ionic_concentrations__Ki);
        double var_rapid_delayed_rectifier_potassium_current__i_Kr = 0.0238624460886 * var_rapid_delayed_rectifier_potassium_current_xr_gate__xr * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_16(_table_index_0, _factor_0) * (var_membrane__V - var_time_independent_potassium_current__E_K);
        double var_membrane__i_Kr = var_rapid_delayed_rectifier_potassium_current__i_Kr;
        double var_slow_delayed_rectifier_potassium_current__i_Ks = (0.433 * (1.0 + (0.6 / (1.0 + pow(3.8e-05 / var_calcium_dynamics__Cai, 1.4))))) * var_slow_delayed_rectifier_potassium_current_xs1_gate__xs1 * var_slow_delayed_rectifier_potassium_current_xs2_gate__xs2 * (var_membrane__V - (26.7123387055 * log(6.91956 / (var_ionic_concentrations__Ki + (0.01833 * var_ionic_concentrations__Nai)))));
        double var_membrane__i_Ks = var_slow_delayed_rectifier_potassium_current__i_Ks;
        double var_sodium_activated_potassium_current__i_K_Na = 0.0 * (0.85 / (1.0 + pow(66.0 / var_ionic_concentrations__Nai, 2.8))) * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_24(_table_index_0, _factor_0) * (var_membrane__V - var_time_independent_potassium_current__E_K);
        double var_membrane__i_K_Na = var_sodium_activated_potassium_current__i_K_Na;
        double var_ATP_sensitive_potassium_current__i_K_ATP = 2.75741043608e-08 * (var_membrane__V - var_time_independent_potassium_current__E_K);
        double var_membrane__i_K_ATP = var_ATP_sensitive_potassium_current__i_K_ATP;
        double var_transient_outward_current__i_to = 0.0 * pow(var_transient_outward_current_zdv_gate__zdv, 3.0) * var_transient_outward_current_ydv_gate__ydv * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_25(_table_index_0, _factor_0) * (var_membrane__V - var_time_independent_potassium_current__E_K);
        double var_membrane__i_to = var_transient_outward_current__i_to;
        double var_Na_Ca_exchanger__i_NaCa = (0.00025 * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_35(_table_index_0, _factor_0) * ((FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_36(_table_index_0, _factor_0) * pow(var_ionic_concentrations__Nai, 3.0) * 1.8) - (2299968.0 * var_calcium_dynamics__Cai))) / (1.0 + (0.0001 * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_37(_table_index_0, _factor_0) * ((FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_38(_table_index_0, _factor_0) * pow(var_ionic_concentrations__Nai, 3.0) * 1.8) + (2299968.0 * var_calcium_dynamics__Cai))));
        double var_membrane__i_NaCa = var_Na_Ca_exchanger__i_NaCa;
        double var_time_independent_potassium_current_K1_gate__alpha_K1 = 1020.0 / (1.0 + exp(0.2385 * ((var_membrane__V - var_time_independent_potassium_current__E_K) - 59.215)));
        double var_time_independent_potassium_current__i_K1 = 0.684653196881 * (var_time_independent_potassium_current_K1_gate__alpha_K1 / (var_time_independent_potassium_current_K1_gate__alpha_K1 + ((1000.0 * ((0.49124 * exp(0.08032 * ((var_membrane__V - var_time_independent_potassium_current__E_K) + 5.476))) + exp(0.06175 * ((var_membrane__V - var_time_independent_potassium_current__E_K) - 594.31)))) / (1.0 + exp( -0.5143 * ((var_membrane__V - var_time_independent_potassium_current__E_K) + 4.753)))))) * (var_membrane__V - var_time_independent_potassium_current__E_K);
        double var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
        double var_plateau_potassium_current__i_Kp = 0.00552 * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_23(_table_index_0, _factor_0) * (var_membrane__V - var_time_independent_potassium_current__E_K);
        double var_membrane__i_Kp = var_plateau_potassium_current__i_Kp;
        double var_sarcolemmal_calcium_pump__i_p_Ca = (1.15 * var_calcium_dynamics__Cai) / (0.0005 + var_calcium_dynamics__Cai);
        double var_membrane__i_p_Ca = var_sarcolemmal_calcium_pump__i_p_Ca;
        double var_sodium_background_current__i_Na_b = 0.004 * (var_membrane__V - var_fast_sodium_current__E_Na);
        double var_membrane__i_Na_b = var_sodium_background_current__i_Na_b;
        double var_calcium_background_current__i_Ca_b = 0.003016 * (var_membrane__V - var_calcium_background_current__E_Ca);
        double var_membrane__i_Ca_b = var_calcium_background_current__i_Ca_b;
        double var_sodium_potassium_pump__i_NaK = ((FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_30(_table_index_0, _factor_0) / (1.0 + pow(10.0 / var_ionic_concentrations__Nai, 2.0))) * 4.5) * 0.166666666667;
        double var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK;
        double var_non_specific_calcium_activated_current__i_ns_Na = (((((0.0 * 1.0 * var_membrane__V * 9309355225.0) * 3.87996927064e-07) * ((0.75 * var_ionic_concentrations__Nai * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_31(_table_index_0, _factor_0)) - 99.0)) / FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_32(_table_index_0, _factor_0)) * 1.0) / (1.0 + pow(0.0012 / var_calcium_dynamics__Cai, 3.0));
        double var_non_specific_calcium_activated_current__i_ns_K = (((((0.0 * 1.0 * var_membrane__V * 9309355225.0) * 3.87996927064e-07) * ((0.75 * var_ionic_concentrations__Ki * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_33(_table_index_0, _factor_0)) - 3.375)) / FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_34(_table_index_0, _factor_0)) * 1.0) / (1.0 + pow(0.0012 / var_calcium_dynamics__Cai, 3.0));
        double var_non_specific_calcium_activated_current__i_ns_Ca = var_non_specific_calcium_activated_current__i_ns_Na + var_non_specific_calcium_activated_current__i_ns_K;
        double var_membrane__i_ns_Ca = var_non_specific_calcium_activated_current__i_ns_Ca;

        return var_membrane__i_Na+var_membrane__i_Ca_L+var_membrane__i_Ca_T+var_membrane__i_Kr+var_membrane__i_Ks+var_membrane__i_K_Na+var_membrane__i_K_ATP+var_membrane__i_to+var_membrane__i_NaCa+var_membrane__i_K1+var_membrane__i_Kp+var_membrane__i_p_Ca+var_membrane__i_Na_b+var_membrane__i_Ca_b+var_membrane__i_NaK+var_membrane__i_ns_Ca;
    }

    void EvaluateYDerivatives (
            double var_environment__time,
            const std::vector<double> &rY,
            std::vector<double> &rDY)
    {
        // Inputs:
        // Time units: second
        var_environment__time *= 1e-3;
        double var_membrane__V = rY[0];
        // Units: millivolt; Initial value: -90
        double var_fast_sodium_current_m_gate__m = rY[1];
        // Units: dimensionless; Initial value: 0.0008
        double var_fast_sodium_current_h_gate__h = rY[2];
        // Units: dimensionless; Initial value: 0.993771
        double var_fast_sodium_current_j_gate__j = rY[3];
        // Units: dimensionless; Initial value: 0.995727
        double var_L_type_Ca_channel_d_gate__d = rY[4];
        // Units: dimensionless; Initial value: 3.210618e-6
        double var_L_type_Ca_channel_f_gate__f = rY[5];
        // Units: dimensionless; Initial value: 0.999837
        double var_T_type_Ca_channel_b_gate__b = rY[6];
        // Units: dimensionless; Initial value: 0.000970231
        double var_T_type_Ca_channel_g_gate__g = rY[7];
        // Units: dimensionless; Initial value: 0.994305
        double var_rapid_delayed_rectifier_potassium_current_xr_gate__xr = rY[8];
        // Units: dimensionless; Initial value: 0.000124042
        double var_slow_delayed_rectifier_potassium_current_xs1_gate__xs1 = rY[9];
        // Units: dimensionless; Initial value: 0.00445683
        double var_slow_delayed_rectifier_potassium_current_xs2_gate__xs2 = rY[10];
        // Units: dimensionless; Initial value: 0.00445683
        double var_transient_outward_current_zdv_gate__zdv = rY[11];
        // Units: dimensionless; Initial value: 0.5
        double var_transient_outward_current_ydv_gate__ydv = rY[12];
        // Units: dimensionless; Initial value: 0.5
        double var_calcium_dynamics__Cai = rY[13];
        // Units: millimolar; Initial value: 6e-5
        double var_calcium_dynamics__Ca_JSR = rY[14];
        // Units: millimolar; Initial value: 1.8
        double var_calcium_dynamics__Ca_NSR = rY[15];
        // Units: millimolar; Initial value: 1.8
        double var_calcium_dynamics__APtrack = rY[16];
        // Units: dimensionless; Initial value: 0
        double var_calcium_dynamics__APtrack2 = rY[17];
        // Units: dimensionless; Initial value: 0
        double var_calcium_dynamics__APtrack3 = rY[18];
        // Units: dimensionless; Initial value: 0
        double var_calcium_dynamics__Cainfluxtrack = rY[19];
        // Units: dimensionless; Initial value: 0
        double var_calcium_dynamics__OVRLDtrack = rY[20];
        // Units: dimensionless; Initial value: 0
        double var_calcium_dynamics__OVRLDtrack2 = rY[21];
        // Units: dimensionless; Initial value: 0
        double var_calcium_dynamics__OVRLDtrack3 = rY[22];
        // Units: dimensionless; Initial value: 0
        double var_ionic_concentrations__Nai = rY[23];
        // Units: millimolar; Initial value: 9
        double var_ionic_concentrations__Ki = rY[24];
        // Units: millimolar; Initial value: 141.2


        // Lookup table indexing
#define COVERAGE_IGNORE
        if (var_membrane__V>99.9999 || var_membrane__V<-100.0001)
            EXCEPTION(DumpState("V outside lookup table range"));
#undef COVERAGE_IGNORE
        double _offset_0 = var_membrane__V - -100.0001;
        double _offset_0_over_table_step = _offset_0 * 100.0;
        unsigned _table_index_0 = (unsigned) floor(_offset_0_over_table_step);
        double _factor_0 = _offset_0_over_table_step - _table_index_0;

        // Mathematics
        double var_membrane__I_st = GetStimulus(var_environment__time*1000);
        double var_fast_sodium_current__E_Na = 26.7123387055 * log(132.0 / var_ionic_concentrations__Nai);
        double var_fast_sodium_current__i_Na = 16.0 * pow(var_fast_sodium_current_m_gate__m, 3.0) * var_fast_sodium_current_h_gate__h * var_fast_sodium_current_j_gate__j * (var_membrane__V - var_fast_sodium_current__E_Na);
        double var_membrane__i_Na = var_fast_sodium_current__i_Na;
        double var_L_type_Ca_channel_f_Ca_gate__f_Ca = 1.0 / (1.0 + (var_calcium_dynamics__Cai * 1666.66666667));
        double var_L_type_Ca_channel__i_CaCa = var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f_Ca_gate__f_Ca * ((((0.00054 * 4.0 * var_membrane__V * 9309355225.0) * 3.87996927064e-07) * ((1.0 * var_calcium_dynamics__Cai * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_5(_table_index_0, _factor_0)) - 0.6138)) / FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_6(_table_index_0, _factor_0));
        double var_L_type_Ca_channel__i_CaNa = var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f_Ca_gate__f_Ca * ((((6.75e-07 * 1.0 * var_membrane__V * 9309355225.0) * 3.87996927064e-07) * ((0.75 * var_ionic_concentrations__Nai * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_7(_table_index_0, _factor_0)) - 99.0)) / FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_8(_table_index_0, _factor_0));
        double var_L_type_Ca_channel__i_CaK = var_L_type_Ca_channel_d_gate__d * var_L_type_Ca_channel_f_gate__f * var_L_type_Ca_channel_f_Ca_gate__f_Ca * ((((1.93e-07 * 1.0 * var_membrane__V * 9309355225.0) * 3.87996927064e-07) * ((0.75 * var_ionic_concentrations__Ki * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_9(_table_index_0, _factor_0)) - 3.375)) / FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_10(_table_index_0, _factor_0));
        double var_L_type_Ca_channel__i_Ca_L = var_L_type_Ca_channel__i_CaCa + var_L_type_Ca_channel__i_CaK + var_L_type_Ca_channel__i_CaNa;
        double var_membrane__i_Ca_L = var_L_type_Ca_channel__i_Ca_L;
        double var_calcium_background_current__E_Ca = 13.3561693527 * log(1.8 / var_calcium_dynamics__Cai);
        double var_T_type_Ca_channel__i_Ca_T = 0.05 * var_T_type_Ca_channel_b_gate__b * var_T_type_Ca_channel_b_gate__b * var_T_type_Ca_channel_g_gate__g * (var_membrane__V - var_calcium_background_current__E_Ca);
        double var_membrane__i_Ca_T = var_T_type_Ca_channel__i_Ca_T;
        double var_time_independent_potassium_current__E_K = 26.7123387055 * log(4.5 / var_ionic_concentrations__Ki);
        double var_rapid_delayed_rectifier_potassium_current__i_Kr = 0.0238624460886 * var_rapid_delayed_rectifier_potassium_current_xr_gate__xr * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_16(_table_index_0, _factor_0) * (var_membrane__V - var_time_independent_potassium_current__E_K);
        double var_membrane__i_Kr = var_rapid_delayed_rectifier_potassium_current__i_Kr;
        double var_slow_delayed_rectifier_potassium_current__i_Ks = (0.433 * (1.0 + (0.6 / (1.0 + pow(3.8e-05 / var_calcium_dynamics__Cai, 1.4))))) * var_slow_delayed_rectifier_potassium_current_xs1_gate__xs1 * var_slow_delayed_rectifier_potassium_current_xs2_gate__xs2 * (var_membrane__V - (26.7123387055 * log(6.91956 / (var_ionic_concentrations__Ki + (0.01833 * var_ionic_concentrations__Nai)))));
        double var_membrane__i_Ks = var_slow_delayed_rectifier_potassium_current__i_Ks;
        double var_sodium_activated_potassium_current__i_K_Na = 0.0 * (0.85 / (1.0 + pow(66.0 / var_ionic_concentrations__Nai, 2.8))) * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_24(_table_index_0, _factor_0) * (var_membrane__V - var_time_independent_potassium_current__E_K);
        double var_membrane__i_K_Na = var_sodium_activated_potassium_current__i_K_Na;
        double var_ATP_sensitive_potassium_current__i_K_ATP = 2.75741043608e-08 * (var_membrane__V - var_time_independent_potassium_current__E_K);
        double var_membrane__i_K_ATP = var_ATP_sensitive_potassium_current__i_K_ATP;
        double var_transient_outward_current__i_to = 0.0 * pow(var_transient_outward_current_zdv_gate__zdv, 3.0) * var_transient_outward_current_ydv_gate__ydv * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_25(_table_index_0, _factor_0) * (var_membrane__V - var_time_independent_potassium_current__E_K);
        double var_membrane__i_to = var_transient_outward_current__i_to;
        double var_Na_Ca_exchanger__i_NaCa = (0.00025 * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_35(_table_index_0, _factor_0) * ((FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_36(_table_index_0, _factor_0) * pow(var_ionic_concentrations__Nai, 3.0) * 1.8) - (2299968.0 * var_calcium_dynamics__Cai))) / (1.0 + (0.0001 * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_37(_table_index_0, _factor_0) * ((FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_38(_table_index_0, _factor_0) * pow(var_ionic_concentrations__Nai, 3.0) * 1.8) + (2299968.0 * var_calcium_dynamics__Cai))));
        double var_membrane__i_NaCa = var_Na_Ca_exchanger__i_NaCa;
        double var_time_independent_potassium_current_K1_gate__alpha_K1 = 1020.0 / (1.0 + exp(0.2385 * ((var_membrane__V - var_time_independent_potassium_current__E_K) - 59.215)));
        double var_time_independent_potassium_current__i_K1 = 0.684653196881 * (var_time_independent_potassium_current_K1_gate__alpha_K1 / (var_time_independent_potassium_current_K1_gate__alpha_K1 + ((1000.0 * ((0.49124 * exp(0.08032 * ((var_membrane__V - var_time_independent_potassium_current__E_K) + 5.476))) + exp(0.06175 * ((var_membrane__V - var_time_independent_potassium_current__E_K) - 594.31)))) / (1.0 + exp( -0.5143 * ((var_membrane__V - var_time_independent_potassium_current__E_K) + 4.753)))))) * (var_membrane__V - var_time_independent_potassium_current__E_K);
        double var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
        double var_plateau_potassium_current__i_Kp = 0.00552 * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_23(_table_index_0, _factor_0) * (var_membrane__V - var_time_independent_potassium_current__E_K);
        double var_membrane__i_Kp = var_plateau_potassium_current__i_Kp;
        double var_sarcolemmal_calcium_pump__i_p_Ca = (1.15 * var_calcium_dynamics__Cai) / (0.0005 + var_calcium_dynamics__Cai);
        double var_membrane__i_p_Ca = var_sarcolemmal_calcium_pump__i_p_Ca;
        double var_sodium_background_current__i_Na_b = 0.004 * (var_membrane__V - var_fast_sodium_current__E_Na);
        double var_membrane__i_Na_b = var_sodium_background_current__i_Na_b;
        double var_calcium_background_current__i_Ca_b = 0.003016 * (var_membrane__V - var_calcium_background_current__E_Ca);
        double var_membrane__i_Ca_b = var_calcium_background_current__i_Ca_b;
        double var_sodium_potassium_pump__i_NaK = ((FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_30(_table_index_0, _factor_0) / (1.0 + pow(10.0 / var_ionic_concentrations__Nai, 2.0))) * 4.5) * 0.166666666667;
        double var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK;
        double var_non_specific_calcium_activated_current__i_ns_Na = (((((0.0 * 1.0 * var_membrane__V * 9309355225.0) * 3.87996927064e-07) * ((0.75 * var_ionic_concentrations__Nai * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_31(_table_index_0, _factor_0)) - 99.0)) / FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_32(_table_index_0, _factor_0)) * 1.0) / (1.0 + pow(0.0012 / var_calcium_dynamics__Cai, 3.0));
        double var_non_specific_calcium_activated_current__i_ns_K = (((((0.0 * 1.0 * var_membrane__V * 9309355225.0) * 3.87996927064e-07) * ((0.75 * var_ionic_concentrations__Ki * FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_33(_table_index_0, _factor_0)) - 3.375)) / FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_34(_table_index_0, _factor_0)) * 1.0) / (1.0 + pow(0.0012 / var_calcium_dynamics__Cai, 3.0));
        double var_non_specific_calcium_activated_current__i_ns_Ca = var_non_specific_calcium_activated_current__i_ns_Na + var_non_specific_calcium_activated_current__i_ns_K;
        double var_membrane__i_ns_Ca = var_non_specific_calcium_activated_current__i_ns_Ca;
        double var_fast_sodium_current_m_gate__E0_m = var_membrane__V + 47.13;
        double var_L_type_Ca_channel_d_gate__E0_d = var_membrane__V + 10.0;
        double var_L_type_Ca_channel_d_gate__d_infinity = 1.0 / (1.0 + exp((-var_L_type_Ca_channel_d_gate__E0_d) * 0.160256410256));
        double var_L_type_Ca_channel_d_gate__tau_d = (fabs(var_L_type_Ca_channel_d_gate__E0_d) < 1e-05) ? 0.00457875457875 : ((0.001 * var_L_type_Ca_channel_d_gate__d_infinity * (1.0 - exp((-var_L_type_Ca_channel_d_gate__E0_d) * 0.160256410256))) / (0.035 * var_L_type_Ca_channel_d_gate__E0_d));
        double var_L_type_Ca_channel_f_gate__f_infinity = FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_11(_table_index_0, _factor_0);
        double var_L_type_Ca_channel_f_gate__tau_f = FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_12(_table_index_0, _factor_0);
        double var_transient_outward_current_zdv_gate__alpha_zdv = FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_26(_table_index_0, _factor_0);
        double var_transient_outward_current_zdv_gate__beta_zdv = FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_27(_table_index_0, _factor_0);
        double var_transient_outward_current_ydv_gate__alpha_ydv = FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_28(_table_index_0, _factor_0);
        double var_transient_outward_current_ydv_gate__beta_ydv = FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_29(_table_index_0, _factor_0);
        double var_calcium_dynamics__i_rel = ((var_calcium_dynamics__Cainfluxtrack > 0.00018) ? (((60000.0 * (var_calcium_dynamics__Cainfluxtrack - 0.00018)) / ((0.0008 + var_calcium_dynamics__Cainfluxtrack) - 0.00018)) * (1.0 - var_calcium_dynamics__APtrack2) * var_calcium_dynamics__APtrack2) : ((var_calcium_dynamics__Cainfluxtrack <= 0.00018) && (var_calcium_dynamics__OVRLDtrack2 > 0.0)) ? (4000.0 * (1.0 - var_calcium_dynamics__OVRLDtrack2) * var_calcium_dynamics__OVRLDtrack2) : 0.0) * (var_calcium_dynamics__Ca_JSR - var_calcium_dynamics__Cai);
        double var_calcium_dynamics__i_up = (8.75 * var_calcium_dynamics__Cai) / (var_calcium_dynamics__Cai + 0.00092);
        double var_calcium_dynamics__i_leak = 0.583333333333 * var_calcium_dynamics__Ca_NSR;
        double var_calcium_dynamics__i_tr = (var_calcium_dynamics__Ca_NSR - var_calcium_dynamics__Ca_JSR) * 5.55555555556;
        double d_dt_membrane__V =  -1000.0 * (var_membrane__i_Na + var_membrane__i_Ca_L + var_membrane__i_Ca_T + var_membrane__i_Kr + var_membrane__i_Ks + var_membrane__i_K_Na + var_membrane__i_K_ATP + var_membrane__i_to + var_membrane__i_K1 + var_membrane__i_Kp + var_membrane__i_NaCa + var_membrane__i_p_Ca + var_membrane__i_Na_b + var_membrane__i_Ca_b + var_membrane__i_NaK + var_membrane__i_ns_Ca + var_membrane__I_st);
        double d_dt_fast_sodium_current_m_gate__m = (((fabs(var_fast_sodium_current_m_gate__E0_m) >= 1e-05) ? ((320.0 * var_fast_sodium_current_m_gate__E0_m) / (1.0 - exp( -0.1 * var_fast_sodium_current_m_gate__E0_m))) : 3200.0) * (1.0 - var_fast_sodium_current_m_gate__m)) - (FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_0(_table_index_0, _factor_0) * var_fast_sodium_current_m_gate__m);
        double d_dt_fast_sodium_current_h_gate__h = (FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_1(_table_index_0, _factor_0) * (1.0 - var_fast_sodium_current_h_gate__h)) - (FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_2(_table_index_0, _factor_0) * var_fast_sodium_current_h_gate__h);
        double d_dt_fast_sodium_current_j_gate__j = (FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_3(_table_index_0, _factor_0) * (1.0 - var_fast_sodium_current_j_gate__j)) - (FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_4(_table_index_0, _factor_0) * var_fast_sodium_current_j_gate__j);
        double d_dt_L_type_Ca_channel_d_gate__d = ((var_L_type_Ca_channel_d_gate__d_infinity / var_L_type_Ca_channel_d_gate__tau_d) * (1.0 - var_L_type_Ca_channel_d_gate__d)) - (((1.0 - var_L_type_Ca_channel_d_gate__d_infinity) / var_L_type_Ca_channel_d_gate__tau_d) * var_L_type_Ca_channel_d_gate__d);
        double d_dt_L_type_Ca_channel_f_gate__f = ((var_L_type_Ca_channel_f_gate__f_infinity / var_L_type_Ca_channel_f_gate__tau_f) * (1.0 - var_L_type_Ca_channel_f_gate__f)) - (((1.0 - var_L_type_Ca_channel_f_gate__f_infinity) / var_L_type_Ca_channel_f_gate__tau_f) * var_L_type_Ca_channel_f_gate__f);
        double d_dt_T_type_Ca_channel_b_gate__b = (FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_13(_table_index_0, _factor_0) - var_T_type_Ca_channel_b_gate__b) / FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_14(_table_index_0, _factor_0);
        double d_dt_T_type_Ca_channel_g_gate__g = (FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_15(_table_index_0, _factor_0) - var_T_type_Ca_channel_g_gate__g) / ((var_membrane__V <= 0.0) ? (( -0.000875 * var_membrane__V) + 0.012) : 0.012);
        double d_dt_rapid_delayed_rectifier_potassium_current_xr_gate__xr = (FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_17(_table_index_0, _factor_0) - var_rapid_delayed_rectifier_potassium_current_xr_gate__xr) / FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_18(_table_index_0, _factor_0);
        double d_dt_slow_delayed_rectifier_potassium_current_xs1_gate__xs1 = (FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_19(_table_index_0, _factor_0) - var_slow_delayed_rectifier_potassium_current_xs1_gate__xs1) / FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_20(_table_index_0, _factor_0);
        double d_dt_slow_delayed_rectifier_potassium_current_xs2_gate__xs2 = (FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_21(_table_index_0, _factor_0) - var_slow_delayed_rectifier_potassium_current_xs2_gate__xs2) / FaberRudy2000Version3OptimisedLookupTables::Instance()->_lookup_22(_table_index_0, _factor_0);
        double d_dt_transient_outward_current_zdv_gate__zdv = ((var_transient_outward_current_zdv_gate__alpha_zdv / (var_transient_outward_current_zdv_gate__alpha_zdv + var_transient_outward_current_zdv_gate__beta_zdv)) - var_transient_outward_current_zdv_gate__zdv) / (1.0 / (var_transient_outward_current_zdv_gate__alpha_zdv + var_transient_outward_current_zdv_gate__beta_zdv));
        double d_dt_transient_outward_current_ydv_gate__ydv = ((var_transient_outward_current_ydv_gate__alpha_ydv / (var_transient_outward_current_ydv_gate__alpha_ydv + var_transient_outward_current_ydv_gate__beta_ydv)) - var_transient_outward_current_ydv_gate__ydv) / (1.0 / (var_transient_outward_current_ydv_gate__alpha_ydv + var_transient_outward_current_ydv_gate__beta_ydv));
        double d_dt_calcium_dynamics__APtrack = (( -1000.0 * (var_membrane__i_Na + var_membrane__i_Ca_L + var_membrane__i_Ca_T + var_membrane__i_Kr + var_membrane__i_Ks + var_membrane__i_K_Na + var_membrane__i_K_ATP + var_membrane__i_to + var_membrane__i_K1 + var_membrane__i_Kp + var_membrane__i_NaCa + var_membrane__i_p_Ca + var_membrane__i_Na_b + var_membrane__i_Ca_b + var_membrane__i_NaK + var_membrane__i_ns_Ca + var_membrane__I_st)) > 150000.0) ? ((100000.0 * (1.0 - var_calcium_dynamics__APtrack)) - (500.0 * var_calcium_dynamics__APtrack)) : ( -500.0 * var_calcium_dynamics__APtrack);
        double d_dt_calcium_dynamics__APtrack2 = ((var_calcium_dynamics__APtrack < 0.2) && (var_calcium_dynamics__APtrack > 0.18)) ? ((100000.0 * (1.0 - var_calcium_dynamics__APtrack2)) - (500.0 * var_calcium_dynamics__APtrack2)) : ( -500.0 * var_calcium_dynamics__APtrack2);
        double d_dt_calcium_dynamics__APtrack3 = ((var_calcium_dynamics__APtrack < 0.2) && (var_calcium_dynamics__APtrack > 0.18)) ? ((100000.0 * (1.0 - var_calcium_dynamics__APtrack3)) - (500.0 * var_calcium_dynamics__APtrack3)) : ( -10.0 * var_calcium_dynamics__APtrack3);
        double d_dt_calcium_dynamics__Cainfluxtrack = (var_calcium_dynamics__APtrack > 0.2) ? (( -1.434e-07 * (((var_L_type_Ca_channel__i_CaCa + var_T_type_Ca_channel__i_Ca_T) - var_Na_Ca_exchanger__i_NaCa) + var_sarcolemmal_calcium_pump__i_p_Ca + var_calcium_background_current__i_Ca_b)) * 200477.689034) : ((var_calcium_dynamics__APtrack2 > 0.01) && (var_calcium_dynamics__APtrack <= 0.2)) ? 0.0 : ( -500.0 * var_calcium_dynamics__Cainfluxtrack);
        double d_dt_calcium_dynamics__OVRLDtrack = (((1.0 / (1.0 + (0.8 / var_calcium_dynamics__Ca_JSR))) > 0.7) && (var_calcium_dynamics__OVRLDtrack3 < 0.37) && (var_calcium_dynamics__APtrack3 < 0.37)) ? (50000.0 * (1.0 - var_calcium_dynamics__OVRLDtrack)) : ( -500.0 * var_calcium_dynamics__OVRLDtrack);
        double d_dt_calcium_dynamics__OVRLDtrack2 = ((var_calcium_dynamics__OVRLDtrack > 0.98) && (var_calcium_dynamics__OVRLDtrack2 < 0.98)) ? (50000.0 * (1.0 - var_calcium_dynamics__OVRLDtrack2)) : ( -500.0 * var_calcium_dynamics__OVRLDtrack2);
        double d_dt_calcium_dynamics__OVRLDtrack3 = ((var_calcium_dynamics__OVRLDtrack > 0.98) && (var_calcium_dynamics__OVRLDtrack3 < 0.98)) ? (50000.0 * (1.0 - var_calcium_dynamics__OVRLDtrack3)) : ( -10.0 * var_calcium_dynamics__OVRLDtrack3);
        double d_dt_calcium_dynamics__Ca_JSR = (1.0 / (1.0 + (8.0 / pow(0.8 + var_calcium_dynamics__Ca_JSR, 2.0)))) * (var_calcium_dynamics__i_tr - var_calcium_dynamics__i_rel);
        double d_dt_calcium_dynamics__Ca_NSR = ((((-var_calcium_dynamics__i_tr) * 1.8246370132e-13) * 476568879782.0) - var_calcium_dynamics__i_leak) + var_calcium_dynamics__i_up;
        double d_dt_calcium_dynamics__Cai = (1.0 / (1.0 + (0.000119 / pow(0.00238 + var_calcium_dynamics__Cai, 2.0)) + (3.5e-05 / pow(0.0005 + var_calcium_dynamics__Cai, 2.0)))) * ((( -1.434e-07 * (((var_L_type_Ca_channel__i_CaCa + var_T_type_Ca_channel__i_Ca_T) - (2.0 * var_Na_Ca_exchanger__i_NaCa)) + var_sarcolemmal_calcium_pump__i_p_Ca + var_calcium_background_current__i_Ca_b)) * 200477.689034) + ((var_calcium_dynamics__i_rel * 1.8246370132e-13) * 38686179652.9) + (((var_calcium_dynamics__i_leak - var_calcium_dynamics__i_up) * 2.09833256519e-12) * 38686179652.9));
        double d_dt_ionic_concentrations__Nai = ((-(var_fast_sodium_current__i_Na + var_L_type_Ca_channel__i_CaNa + var_sodium_background_current__i_Na_b + var_non_specific_calcium_activated_current__i_ns_Na + (var_Na_Ca_exchanger__i_NaCa * 3.0) + (var_sodium_potassium_pump__i_NaK * 3.0))) * 1.434e-07) * 400955.378068;
        double d_dt_ionic_concentrations__Ki = ((-(var_L_type_Ca_channel__i_CaK + var_rapid_delayed_rectifier_potassium_current__i_Kr + var_slow_delayed_rectifier_potassium_current__i_Ks + var_time_independent_potassium_current__i_K1 + var_plateau_potassium_current__i_Kp + var_sodium_activated_potassium_current__i_K_Na + var_ATP_sensitive_potassium_current__i_K_ATP + var_transient_outward_current__i_to + var_non_specific_calcium_activated_current__i_ns_K + ((-var_sodium_potassium_pump__i_NaK) * 2.0))) * 1.434e-07) * 400955.378068;

        rDY[0] = d_dt_membrane__V*1e-3;
        rDY[1] = d_dt_fast_sodium_current_m_gate__m*1e-3;
        rDY[2] = d_dt_fast_sodium_current_h_gate__h*1e-3;
        rDY[3] = d_dt_fast_sodium_current_j_gate__j*1e-3;
        rDY[4] = d_dt_L_type_Ca_channel_d_gate__d*1e-3;
        rDY[5] = d_dt_L_type_Ca_channel_f_gate__f*1e-3;
        rDY[6] = d_dt_T_type_Ca_channel_b_gate__b*1e-3;
        rDY[7] = d_dt_T_type_Ca_channel_g_gate__g*1e-3;
        rDY[8] = d_dt_rapid_delayed_rectifier_potassium_current_xr_gate__xr*1e-3;
        rDY[9] = d_dt_slow_delayed_rectifier_potassium_current_xs1_gate__xs1*1e-3;
        rDY[10] = d_dt_slow_delayed_rectifier_potassium_current_xs2_gate__xs2*1e-3;
        rDY[11] = d_dt_transient_outward_current_zdv_gate__zdv*1e-3;
        rDY[12] = d_dt_transient_outward_current_ydv_gate__ydv*1e-3;
        rDY[13] = d_dt_calcium_dynamics__Cai*1e-3;
        rDY[14] = d_dt_calcium_dynamics__Ca_JSR*1e-3;
        rDY[15] = d_dt_calcium_dynamics__Ca_NSR*1e-3;
        rDY[16] = d_dt_calcium_dynamics__APtrack*1e-3;
        rDY[17] = d_dt_calcium_dynamics__APtrack2*1e-3;
        rDY[18] = d_dt_calcium_dynamics__APtrack3*1e-3;
        rDY[19] = d_dt_calcium_dynamics__Cainfluxtrack*1e-3;
        rDY[20] = d_dt_calcium_dynamics__OVRLDtrack*1e-3;
        rDY[21] = d_dt_calcium_dynamics__OVRLDtrack2*1e-3;
        rDY[22] = d_dt_calcium_dynamics__OVRLDtrack3*1e-3;
        rDY[23] = d_dt_ionic_concentrations__Nai*1e-3;
        rDY[24] = d_dt_ionic_concentrations__Ki*1e-3;
    }

};

#endif
