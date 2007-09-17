#ifndef NHSCELLULARMECHANICSODESYSTEM_HPP_
#define NHSCELLULARMECHANICSODESYSTEM_HPP_

#include "AbstractOdeSystem.hpp"

///\todo: doxy

class NHSCellularMechanicsOdeSystem  : public AbstractOdeSystem
{
friend class TestCellularMechanicsOdeSystems;

protected :
    double mLambda1;
    double mDLambda1Dt;
    double mCalciumI;
    double mCalciumTrop50; // only dependent on constants and lambda, so updated whenever lambda is updated

    double mK1;
    double mK2;

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
    
    void CalculateCalciumTrop50();

    double CalculateT0(double z);
     
public : 
    NHSCellularMechanicsOdeSystem();
    
    void SetLambda1AndDerivative(double lambda1, double dlambda1Dt);    

    void SetIntracellularCalciumConcentration(double calciumI);

    double GetCalciumTroponinValue();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    
    double GetActiveTension();
};
#endif /*NHSCELLULARMECHANICSODESYSTEM_HPP_*/
