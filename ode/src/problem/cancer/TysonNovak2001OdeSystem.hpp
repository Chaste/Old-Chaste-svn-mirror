#ifndef _TYSONNOVAK2001ODESYSTEM_HPP_
#define _TYSONNOVAK2001ODESYSTEM_HPP_

#include <vector>
#include <cmath>
#include "AbstractOdeSystemWithAnalyticJacobian.hpp"

/**
 * Represents the Tyson & Novak 2001 system of ODEs.
 */
class TysonNovak2001OdeSystem: public AbstractOdeSystemWithAnalyticJacobian
{
private:

    // Constants for the Tyson Novak Model
    double mK1;
    double mK2d;
    double mK2dd;
    double mK2ddd;
    double mCycB_threshold;
    double mK3d;
    double mK3dd;
    double mK4d;
    double mK4;
    double mJ3;
    double mJ4;
    double mK5d;
    double mK5dd;
    double mK6;
    double mJ5;
    double mN;
    double mK7;
    double mK8;
    double mJ7;
    double mJ8;
    double mMad;
    double mK9;
    double mK10;
    double mK11;
    double mK12d;
    double mK12dd;
    double mK12ddd;
    double mKeq;
    double mK13;
    double mK14;
    double mK15d;
    double mK15dd;
    double mK16d;
    double mK16dd;
    double mJ15;
    double mJ16;
    double mMu;
    double mMstar;
    
public:
    // Constructor
    TysonNovak2001OdeSystem();
    
    // Destructor
    ~TysonNovak2001OdeSystem();
    
    void Init(); //Sets up the parameter values
    
    // Compute the RHS of the T&N system of ODEs
    std::vector<double> EvaluateYDerivatives(double time, const std::vector<double> &rY);
    
    bool CalculateStoppingEvent(double time, const std::vector<double> &rY)
    {
         return (fabs(rY[0]-0.1) < 1.0e-2 && EvaluateYDerivatives(time, rY)[0] < 0.0);
    }
    
    PetscErrorCode AnalyticJacobian(Vec solutionGuess, Mat *pJacobian, double time, double timeStep);
};

#endif //_TYSONNOVAK2001ODESYSTEM_HPP_
