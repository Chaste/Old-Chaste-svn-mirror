#ifndef CHECKMONOLR91VARS_HPP_
#define CHECKMONOLR91VARS_HPP_

#include "MonodomainProblem.hpp"

template<unsigned SPACE_DIM>
void CheckMonoLr91Vars(MonodomainProblem<SPACE_DIM>& problem)
{

    DistributedVector voltage(problem.GetVoltage());
    for (DistributedVector::Iterator index = DistributedVector::Begin();
         index != DistributedVector::End();
         ++index)
        {
        // assuming LR model has Ena = 54.4 and Ek = -77
        double Ena   =  54.4;
        double Ek    = -77.0;
        
        TS_ASSERT_LESS_THAN_EQUALS( voltage[index] , Ena +  30);
        TS_ASSERT_LESS_THAN_EQUALS(-voltage[index] + (Ek-30), 0);
        
        std::vector<double> odeVars = problem.GetMonodomainPde()->GetCardiacCell(index.Global)->rGetStateVariables();
        for (int j=0; j<8; j++)
        {
            // if not voltage or calcium ion conc, test whether between 0 and 1
            if ((j!=4) && (j!=3))
            {
                TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);
                TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);
            }
        }
    }
};

#endif /*CHECKMONOLR91VARS_HPP_*/
