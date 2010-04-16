#ifdef CHASTE_CVODE
//! @file
//! 
//! This source file was generated from CellML.
//! 
//! Model: luo_rudy_1991
//! 
//! Processed by pycml - CellML Tools in Python
//!     (translate: 8617, pycml: 8196)
//! on Fri Apr 16 13:33:46 2010
//! 
//! <autogenerated>

#include "luo_rudy_1991CvodeOpt.hpp"
#include <cmath>
#include <cassert>
#include "Exception.hpp"
#include "OdeSystemInformation.hpp"

class Cellluo_rudy_1991FromCellMLCvodeOpt_LookupTables
{
public:
    static Cellluo_rudy_1991FromCellMLCvodeOpt_LookupTables* Instance()
    {
        if (mpInstance == NULL)
        {
            mpInstance = new Cellluo_rudy_1991FromCellMLCvodeOpt_LookupTables;
        }
        return mpInstance;
    }
    
    // Row lookup methods
    // using linear interpolation
    double* _lookup_0_row(unsigned i, double factor)
    {
        for (unsigned j=0; j<16; j++)
        {
            double y1 = _lookup_table_0[i][j];
            double y2 = _lookup_table_0[i+1][j];
            _lookup_table_0_row[j] = y1 + (y2-y1)*factor;
        }
        return _lookup_table_0_row;
    }
    
    
protected:
    Cellluo_rudy_1991FromCellMLCvodeOpt_LookupTables(const Cellluo_rudy_1991FromCellMLCvodeOpt_LookupTables&);
    Cellluo_rudy_1991FromCellMLCvodeOpt_LookupTables& operator= (const Cellluo_rudy_1991FromCellMLCvodeOpt_LookupTables&);
    Cellluo_rudy_1991FromCellMLCvodeOpt_LookupTables()
    {
        assert(mpInstance == NULL);
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][0] = (0.32 * (var_membrane__V + 47.13)) / (1.0 - exp( -0.1 * (var_membrane__V + 47.13)));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][1] = 0.08 * exp((-var_membrane__V) * 0.0909090909091);
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][2] = (var_membrane__V <  -40.0) ? (0.135 * exp((80.0 + var_membrane__V) *  -0.147058823529)) : 0.0;
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][3] = (var_membrane__V <  -40.0) ? ((3.56 * exp(0.079 * var_membrane__V)) + (310000.0 * exp(0.35 * var_membrane__V))) : (1.0 / (0.13 * (1.0 + exp((var_membrane__V + 10.66) *  -0.0900900900901))));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][4] = (var_membrane__V <  -40.0) ? (((( -127140.0 * exp(0.2444 * var_membrane__V)) - (3.474e-05 * exp( -0.04391 * var_membrane__V))) * (var_membrane__V + 37.78)) / (1.0 + exp(0.311 * (var_membrane__V + 79.23)))) : 0.0;
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][5] = (var_membrane__V <  -40.0) ? ((0.1212 * exp( -0.01052 * var_membrane__V)) / (1.0 + exp( -0.1378 * (var_membrane__V + 40.14)))) : ((0.3 * exp( -2.535e-07 * var_membrane__V)) / (1.0 + exp( -0.1 * (var_membrane__V + 32.0))));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][6] = (0.095 * exp( -0.01 * (var_membrane__V - 5.0))) / (1.0 + exp( -0.072 * (var_membrane__V - 5.0)));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][7] = (0.07 * exp( -0.017 * (var_membrane__V + 44.0))) / (1.0 + exp(0.05 * (var_membrane__V + 44.0)));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][8] = (0.012 * exp( -0.008 * (var_membrane__V + 28.0))) / (1.0 + exp(0.15 * (var_membrane__V + 28.0)));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][9] = (0.0065 * exp( -0.02 * (var_membrane__V + 30.0))) / (1.0 + exp( -0.2 * (var_membrane__V + 30.0)));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][10] = (var_membrane__V >  -100.0) ? ((2.837 * (exp(0.04 * (var_membrane__V + 77.0)) - 1.0)) / ((var_membrane__V + 77.0) * exp(0.04 * (var_membrane__V + 35.0)))) : 1.0;
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][11] = (0.0005 * exp(0.083 * (var_membrane__V + 50.0))) / (1.0 + exp(0.057 * (var_membrane__V + 50.0)));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][12] = (0.0013 * exp( -0.06 * (var_membrane__V + 20.0))) / (1.0 + exp( -0.04 * (var_membrane__V + 20.0)));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][13] = ((0.49124 * exp(0.08032 * ((var_membrane__V + 5.476) -  -87.8929017138))) + (1.0 * exp(0.06175 * (var_membrane__V - 506.417098286)))) / (1.0 + exp( -0.5143 * ((var_membrane__V -  -87.8929017138) + 4.753)));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][14] = 1.02 / (1.0 + exp(0.2385 * ((var_membrane__V -  -87.8929017138) - 59.215)));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][15] = 0.0183 * (1.0 / (1.0 + exp((7.488 - var_membrane__V) * 0.167224080268))) * (var_membrane__V -  -87.8929017138);
        }
        
    }
    
private:
    /** The single instance of the class */
    static Cellluo_rudy_1991FromCellMLCvodeOpt_LookupTables *mpInstance;

    // Row lookup methods memory
    double _lookup_table_0_row[16];
    
    // Lookup tables
    double _lookup_table_0[16001][16];
    
};

Cellluo_rudy_1991FromCellMLCvodeOpt_LookupTables* Cellluo_rudy_1991FromCellMLCvodeOpt_LookupTables::mpInstance = NULL;

    double Cellluo_rudy_1991FromCellMLCvodeOpt::Get_membrane__I_stim()
    {
        return var_membrane__I_stim;
    }
    
    double Cellluo_rudy_1991FromCellMLCvodeOpt::Get_membrane__i_Na()
    {
        return var_membrane__i_Na;
    }
    
    double Cellluo_rudy_1991FromCellMLCvodeOpt::Get_membrane__i_si()
    {
        return var_membrane__i_si;
    }
    
    double Cellluo_rudy_1991FromCellMLCvodeOpt::Get_membrane__i_K()
    {
        return var_membrane__i_K;
    }
    
    double Cellluo_rudy_1991FromCellMLCvodeOpt::Get_membrane__i_K1()
    {
        return var_membrane__i_K1;
    }
    
    double Cellluo_rudy_1991FromCellMLCvodeOpt::Get_membrane__i_Kp()
    {
        return var_membrane__i_Kp;
    }
    
    double Cellluo_rudy_1991FromCellMLCvodeOpt::Get_membrane__i_b()
    {
        return var_membrane__i_b;
    }
    
    Cellluo_rudy_1991FromCellMLCvodeOpt::Cellluo_rudy_1991FromCellMLCvodeOpt(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver /* unused; should be empty */, boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
        : AbstractCvodeCell(
                pOdeSolver,
                8,
                0,
                pIntracellularStimulus)
    {
        // Time units: millisecond
        // 
        mpSystemInfo = OdeSystemInformation<Cellluo_rudy_1991FromCellMLCvodeOpt>::Instance();
        Init();

        NV_Ith_S(mParameters, 0) = 23;// var_fast_sodium_current__g_Na
    }
    
    Cellluo_rudy_1991FromCellMLCvodeOpt::~Cellluo_rudy_1991FromCellMLCvodeOpt()
    {
    }
    
    void Cellluo_rudy_1991FromCellMLCvodeOpt::VerifyStateVariables()
    {}

    double Cellluo_rudy_1991FromCellMLCvodeOpt::GetIIonic()
    {
        N_Vector rY = rGetStateVariables();
        double var_membrane__V = NV_Ith_S(rY, 0);
        // Units: millivolt; Initial value: -83.853
        double var_fast_sodium_current_m_gate__m = NV_Ith_S(rY, 1);
        // Units: dimensionless; Initial value: 0.00187018
        double var_fast_sodium_current_h_gate__h = NV_Ith_S(rY, 2);
        // Units: dimensionless; Initial value: 0.9804713
        double var_fast_sodium_current_j_gate__j = NV_Ith_S(rY, 3);
        // Units: dimensionless; Initial value: 0.98767124
        double var_slow_inward_current_d_gate__d = NV_Ith_S(rY, 4);
        // Units: dimensionless; Initial value: 0.00316354
        double var_slow_inward_current_f_gate__f = NV_Ith_S(rY, 5);
        // Units: dimensionless; Initial value: 0.99427859
        double var_time_dependent_potassium_current_X_gate__X = NV_Ith_S(rY, 6);
        // Units: dimensionless; Initial value: 0.16647703
        double var_intracellular_calcium_concentration__Cai = NV_Ith_S(rY, 7);
        // Units: millimolar; Initial value: 0.0002
        
        // Lookup table indexing
        if (var_membrane__V>59.9999 || var_membrane__V<-100.0001)
        {
#define COVERAGE_IGNORE
            EXCEPTION(DumpState("V outside lookup table range", rY));
#undef COVERAGE_IGNORE
        }
        
        double _offset_0 = var_membrane__V - -100.0001;
        double _offset_0_over_table_step = _offset_0 * 100.0;
        unsigned _table_index_0 = (unsigned)(_offset_0_over_table_step);
        double _factor_0 = _offset_0_over_table_step - _table_index_0;
        double* _lt_0_row = Cellluo_rudy_1991FromCellMLCvodeOpt_LookupTables::Instance()->_lookup_0_row(_table_index_0, _factor_0);
        
        const double var_fast_sodium_current__g_Na = NV_Ith_S(mParameters, 0);
        var_membrane__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current_m_gate__m, 3.0) * var_fast_sodium_current_h_gate__h * var_fast_sodium_current_j_gate__j * (var_membrane__V - 54.7944639351);
        const double var_slow_inward_current__i_si = 0.09 * var_slow_inward_current_d_gate__d * var_slow_inward_current_f_gate__f * (var_membrane__V - (7.7 - (13.0287 * log(var_intracellular_calcium_concentration__Cai * 1.0))));
        var_membrane__i_si = var_slow_inward_current__i_si;
        var_membrane__i_K = 0.282 * var_time_dependent_potassium_current_X_gate__X * _lt_0_row[10] * (var_membrane__V -  -77.5675843853);
        const double var_time_independent_potassium_current_K1_gate__alpha_K1 = _lt_0_row[14];
        var_membrane__i_K1 = 0.6047 * (var_time_independent_potassium_current_K1_gate__alpha_K1 / (var_time_independent_potassium_current_K1_gate__alpha_K1 + _lt_0_row[13])) * (var_membrane__V -  -87.8929017138);
        var_membrane__i_Kp = _lt_0_row[15];
        var_membrane__i_b = 0.03921 * (var_membrane__V -  -59.87);
        
        return (var_membrane__i_Na+var_membrane__i_si+var_membrane__i_K+var_membrane__i_K1+var_membrane__i_Kp+var_membrane__i_b);
    }
    
    void Cellluo_rudy_1991FromCellMLCvodeOpt::EvaluateRhs(double var_environment__time, const N_Vector rY, N_Vector rDY)
    {
        // Inputs:
        // Time units: millisecond
        var_environment__time *= 1.0;
        double var_membrane__V = NV_Ith_S(rY, 0);
        // Units: millivolt; Initial value: -83.853
        double var_fast_sodium_current_m_gate__m = NV_Ith_S(rY, 1);
        // Units: dimensionless; Initial value: 0.00187018
        double var_fast_sodium_current_h_gate__h = NV_Ith_S(rY, 2);
        // Units: dimensionless; Initial value: 0.9804713
        double var_fast_sodium_current_j_gate__j = NV_Ith_S(rY, 3);
        // Units: dimensionless; Initial value: 0.98767124
        double var_slow_inward_current_d_gate__d = NV_Ith_S(rY, 4);
        // Units: dimensionless; Initial value: 0.00316354
        double var_slow_inward_current_f_gate__f = NV_Ith_S(rY, 5);
        // Units: dimensionless; Initial value: 0.99427859
        double var_time_dependent_potassium_current_X_gate__X = NV_Ith_S(rY, 6);
        // Units: dimensionless; Initial value: 0.16647703
        double var_intracellular_calcium_concentration__Cai = NV_Ith_S(rY, 7);
        // Units: millimolar; Initial value: 0.0002
        
        
        // Lookup table indexing
        if (var_membrane__V>59.9999 || var_membrane__V<-100.0001)
        {
#define COVERAGE_IGNORE
            EXCEPTION(DumpState("V outside lookup table range", rY));
#undef COVERAGE_IGNORE
        }
        
        double _offset_0 = var_membrane__V - -100.0001;
        double _offset_0_over_table_step = _offset_0 * 100.0;
        unsigned _table_index_0 = (unsigned)(_offset_0_over_table_step);
        double _factor_0 = _offset_0_over_table_step - _table_index_0;
        double* _lt_0_row = Cellluo_rudy_1991FromCellMLCvodeOpt_LookupTables::Instance()->_lookup_0_row(_table_index_0, _factor_0);
        
        // Mathematics
        var_membrane__I_stim = GetStimulus((1.0/1)*var_environment__time);
        const double var_fast_sodium_current__g_Na = NV_Ith_S(mParameters, 0);
        var_membrane__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current_m_gate__m, 3.0) * var_fast_sodium_current_h_gate__h * var_fast_sodium_current_j_gate__j * (var_membrane__V - 54.7944639351);
        const double var_slow_inward_current__i_si = 0.09 * var_slow_inward_current_d_gate__d * var_slow_inward_current_f_gate__f * (var_membrane__V - (7.7 - (13.0287 * log(var_intracellular_calcium_concentration__Cai * 1.0))));
        var_membrane__i_si = var_slow_inward_current__i_si;
        var_membrane__i_K = 0.282 * var_time_dependent_potassium_current_X_gate__X * _lt_0_row[10] * (var_membrane__V -  -77.5675843853);
        const double var_time_independent_potassium_current_K1_gate__alpha_K1 = _lt_0_row[14];
        var_membrane__i_K1 = 0.6047 * (var_time_independent_potassium_current_K1_gate__alpha_K1 / (var_time_independent_potassium_current_K1_gate__alpha_K1 + _lt_0_row[13])) * (var_membrane__V -  -87.8929017138);
        var_membrane__i_Kp = _lt_0_row[15];
        var_membrane__i_b = 0.03921 * (var_membrane__V -  -59.87);
        
        double d_dt_membrane__V;
        if (mSetVoltageDerivativeToZero)
        {
            d_dt_membrane__V = 0.0;
        }
        else
        {
            d_dt_membrane__V =  -1.0 * (var_membrane__I_stim + var_membrane__i_Na + var_membrane__i_si + var_membrane__i_K + var_membrane__i_K1 + var_membrane__i_Kp + var_membrane__i_b);
        }
        
        const double d_dt_fast_sodium_current_m_gate__m = (_lt_0_row[0] * (1.0 - var_fast_sodium_current_m_gate__m)) - (_lt_0_row[1] * var_fast_sodium_current_m_gate__m);
        const double d_dt_fast_sodium_current_h_gate__h = (_lt_0_row[2] * (1.0 - var_fast_sodium_current_h_gate__h)) - (_lt_0_row[3] * var_fast_sodium_current_h_gate__h);
        const double d_dt_fast_sodium_current_j_gate__j = (_lt_0_row[4] * (1.0 - var_fast_sodium_current_j_gate__j)) - (_lt_0_row[5] * var_fast_sodium_current_j_gate__j);
        const double d_dt_slow_inward_current_d_gate__d = (_lt_0_row[6] * (1.0 - var_slow_inward_current_d_gate__d)) - (_lt_0_row[7] * var_slow_inward_current_d_gate__d);
        const double d_dt_slow_inward_current_f_gate__f = (_lt_0_row[8] * (1.0 - var_slow_inward_current_f_gate__f)) - (_lt_0_row[9] * var_slow_inward_current_f_gate__f);
        const double d_dt_time_dependent_potassium_current_X_gate__X = (_lt_0_row[11] * (1.0 - var_time_dependent_potassium_current_X_gate__X)) - (_lt_0_row[12] * var_time_dependent_potassium_current_X_gate__X);
        const double d_dt_intracellular_calcium_concentration__Cai = ( -0.0001 * var_slow_inward_current__i_si) + (0.07 * (0.0001 - var_intracellular_calcium_concentration__Cai));
        
        NV_Ith_S(rDY, 0) = 1.0*d_dt_membrane__V;
        NV_Ith_S(rDY, 1) = 1.0*d_dt_fast_sodium_current_m_gate__m;
        NV_Ith_S(rDY, 2) = 1.0*d_dt_fast_sodium_current_h_gate__h;
        NV_Ith_S(rDY, 3) = 1.0*d_dt_fast_sodium_current_j_gate__j;
        NV_Ith_S(rDY, 4) = 1.0*d_dt_slow_inward_current_d_gate__d;
        NV_Ith_S(rDY, 5) = 1.0*d_dt_slow_inward_current_f_gate__f;
        NV_Ith_S(rDY, 6) = 1.0*d_dt_time_dependent_potassium_current_X_gate__X;
        NV_Ith_S(rDY, 7) = 1.0*d_dt_intracellular_calcium_concentration__Cai;
    }
    
template<>
void OdeSystemInformation<Cellluo_rudy_1991FromCellMLCvodeOpt>::Initialise(void)
{
    // Time units: millisecond
    // 
    this->mVariableNames.push_back("membrane__V");
    this->mVariableUnits.push_back("millivolt");
    this->mInitialConditions.push_back(-83.853);

    this->mVariableNames.push_back("fast_sodium_current_m_gate__m");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.00187018);

    this->mVariableNames.push_back("fast_sodium_current_h_gate__h");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.9804713);

    this->mVariableNames.push_back("fast_sodium_current_j_gate__j");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.98767124);

    this->mVariableNames.push_back("slow_inward_current_d_gate__d");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.00316354);

    this->mVariableNames.push_back("slow_inward_current_f_gate__f");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.99427859);

    this->mVariableNames.push_back("time_dependent_potassium_current_X_gate__X");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.16647703);

    this->mVariableNames.push_back("intracellular_calcium_concentration__Cai");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.0002);

    this->mParameterNames.push_back("fast_sodium_current__g_Na");
    this->mParameterUnits.push_back("milliS_per_cm2");
    
    this->mInitialised = true;
}


#endif // CHASTE_CVODE
