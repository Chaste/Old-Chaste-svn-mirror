#ifndef KSPCONVERGENCETESTER_HPP_
#define KSPCONVERGENCETESTER_HPP_

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM>
class KspConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM>
{
public:
    void SetInitialConvergenceParameters()
    {
        this->SetKspRelativeTolerance(1e-2);
    }
    void UpdateConvergenceParameters()
    {
        this->SetKspRelativeTolerance(this->GetKspRelativeTolerance()*0.1);
    
    }
    bool GiveUpConvergence()
    {
        return this->GetKspRelativeTolerance()<1e-9;
    }
    double Abscissa()
    {
        return this->GetKspRelativeTolerance();
    }
};

#endif /*KSPCONVERGENCETESTER_HPP_*/
