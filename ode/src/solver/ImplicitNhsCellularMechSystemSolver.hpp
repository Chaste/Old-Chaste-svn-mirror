#ifndef IMPLICITNHSCELLULARMECHSYSTEMSOLVER_HPP_
#define IMPLICITNHSCELLULARMECHSYSTEMSOLVER_HPP_

#include "NHSCellularMechanicsOdeSystem.hpp"
#include "EulerIvpOdeSolver.hpp"

class ImplicitNhsCellularMechSystemSolver
{
private:
    NHSCellularMechanicsOdeSystem& mrSystem;    
    std::vector<double> mStoredStateVariables;
    
public :  
    ImplicitNhsCellularMechSystemSolver(NHSCellularMechanicsOdeSystem& rSystem)
        : mrSystem(rSystem)
    {
        mStoredStateVariables.resize(rSystem.GetNumberOfStateVariables());
    }
    
    void SolveDoNotUpdate(double startTime, 
                          double endTime, 
                          double timestep)
    {
        assert(startTime < endTime);
        
        // store the initial state vars
        std::vector<double> old_state_vars = mrSystem.rGetStateVariables();
        
//todo: proper implicit solve
        // solve
        EulerIvpOdeSolver solver;
        solver.SolveAndUpdateStateVariable(&mrSystem, startTime, endTime, timestep);


        // store the new state vars
        mStoredStateVariables = mrSystem.rGetStateVariables();
        
        // rewrite with the old state vars
        for(unsigned i=0; i<mrSystem.rGetStateVariables().size(); i++)
        {
            mrSystem.rGetStateVariables()[i] = old_state_vars[i];
        }        
    }
    
    void UpdateSystem()
    {
        for(unsigned i=0; i<mrSystem.rGetStateVariables().size(); i++)
        {
            mrSystem.rGetStateVariables()[i] = mStoredStateVariables[i];
        }      
    }
};

#endif /*IMPLICITNHSCELLULARMECHSYSTEMSOLVER_HPP_*/
