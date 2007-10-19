#ifndef NHSCELLULARMECHANICSODESYSTEM_HPP_
#define NHSCELLULARMECHANICSODESYSTEM_HPP_

#include "AbstractOdeSystem.hpp"


/**
 *  NHS (Niederer, Hunter, Smith) model of active tension in cardiac cells.
 * 
 *  A system of ODEs which determines the active potential, given the intracellular
 *  calcium concentration, the stretch (lambda) of the cell, and the stretch rate 
 *  (dlambda_dt) of the cell.
 * 
 *  The state variables are, in order: Calcium_troponin, z, Q1, Q2, Q3
 * 
 *  Reference: S.A. Niederer, N.P. Smith, P.J. Hunter, "New developments in a strongly
 *  coupled cardiac electro-mechanical model" Europace 7, S118-S127
 * 
 *  THE ACTIVE TENSION IS RETURNED IN KILOPASCALS!!
 */
class NhsCellularMechanicsOdeSystem  : public AbstractOdeSystem
{
friend class TestCellularMechanicsOdeSystems;

protected :
    /*< The stretch. To be specified by the caller */
    double mLambda;
    /*< The stretch rate. To be specified by the caller */
    double mDLambdaDt;
    /*< The intracellular calcium concentration. To be specified by the caller */
    double mCalciumI;


    /*< A parameter only dependent on constants and lambda, so updated whenever lambda is updated */
    double mCalciumTrop50; 

    /*< A constant determined from the other constrants. Set up in the constructor */ 
    double mK1;
    /*< A constant determined from the other constrants. Set up in the constructor */ 
    double mK2;

    // Parameters  
    
    /** FILL IN. (mMols)^-1 (ms)^-1 */
    static const double mKon = 100;

    /** FILL IN. (ms)^-1 */
    static const double mKrefoff = 0.2;

    /** FILL IN. Dimensionless */
    static const double mGamma = 2;

    /** FILL IN. mMols */
    static const double mCalciumTroponinMax = 0.07;

    /** FILL IN. (ms)^-1 */
    static const double mAlphaR1 = 0.002;

    /** FILL IN. (ms)^-1 */
    static const double mAlphaR2 = 0.00175;

    /** FILL IN. Dimensionless */
    static const double mKZ = 0.15;

    /** FILL IN. Dimensionless */
    static const double mNr = 3;

    /** FILL IN. Dimensionless */
    static const double mBeta1 = -4;

    /** FILL IN. (ms)^-1 */
    static const double mAlpha0 = 0.008;

    /** FILL IN. Dimensionless */
    static const double mN = 3;

    /** FILL IN. Dimensionless */
    static const double mZp = 0.85;

    /** FILL IN. mMols */
    static const double mCalcium50ref = 0.00105;

    /** FILL IN. kPa */
    static const double mTref = 56.2;

    /** FILL IN. Dimensionless */
    static const double mBeta0 = 4.9;

    /** FILL IN. Dimensionless */
    static const double mA = 0.35;  

    /** FILL IN. Dimensionless */
    static const double mA1 = -29;

    /** FILL IN. Dimensionless */
    static const double mA2 = 138;

    /** FILL IN. Dimensionless */
    static const double mA3 = 129;

    /** FILL IN. (ms)^-1 */
    static const double mAlpha1 = 0.03;

    /** FILL IN. (ms)^-1 */
    static const double mAlpha2 = 0.130;

    /** FILL IN. (ms)^-1 */
    static const double mAlpha3 = 0.625;

    /**  
     *  Compute the calcium_trop50 concentration. This is a function of constants and 
     *  lambda, so only needs to be called in the constructor or when lambda is set
     */
    void CalculateCalciumTrop50();

    /**
     *  Calculate T0. This is a function of constants, lambda and z
     */
    double CalculateT0(double z);
     
public :
    /** 
     *  Constructor. Initialises all state variables to zero, lambda to 1, dlambda_dt
     *  to 0 and intracellular calcium concentration to 0
     */ 
    NhsCellularMechanicsOdeSystem();
    
    /**
     *  Set the current stretch and the stretch rate of the cell/fibre
     */
    void SetLambdaAndDerivative(double lambda, double dlambdaDt);    

    /** 
     *  Set the current intracellular calcium concentration
     */
    void SetIntracellularCalciumConcentration(double calciumI);

    /**
     *  Get the current Calcium Troponin (one of the state variables) value. This
     *  may be needed if the cell model has Calcium troponin and might need overwriting
     */
    double GetCalciumTroponinValue();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    
    /**
     *  Get the active tension, which is a function of the constants and current state variables
     */
    double GetActiveTension();
    
    /**
     *  Get the current stretch rate
     */
    double GetLambda();
};
#endif /*NHSCELLULARMECHANICSODESYSTEM_HPP_*/
