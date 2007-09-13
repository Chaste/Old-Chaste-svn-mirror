#ifndef KSPCONVERGENCETESTER_HPP_
#define KSPCONVERGENCETESTER_HPP_

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM>
class KspConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM>
{
public:
    void SetInitialConvergenceParameters()
    {
        this->KspRtol=1e-2;
    }
    void UpdateConvergenceParameters()
    {
        this->KspRtol *= 0.1;
    
    }
    bool GiveUpConvergence()
    {
        return this->KspRtol<1e-9;
    }
    double Abscissa()
    {
        return this->KspRtol;
    }
};

#endif /*KSPCONVERGENCETESTER_HPP_*/
