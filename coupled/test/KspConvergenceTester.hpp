#ifndef KSPCONVERGENCETESTER_HPP_
#define KSPCONVERGENCETESTER_HPP_

template<class CELL, class CARDIAC_PROBLEM, unsigned DIM>
class KspConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM>
{
public:
    void SetInitialConvergenceParameters()
    {
        this->mKspRtol=1e-1;
    }
    void UpdateConvergenceParameters()
    {
        this->mKspRtol *= 0.1;
    
    }
    bool GiveUpConvergence()
    {
        return this->mKspRtol<1e-9;
    }
    double Abscissa()
    {
        return this->mKspRtol;
    }
};

#endif /*KSPCONVERGENCETESTER_HPP_*/
