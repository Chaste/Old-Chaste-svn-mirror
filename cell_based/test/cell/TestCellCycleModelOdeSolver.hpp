/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef TESTCELLCYCLEMODELODESOLVER_HPP_
#define TESTCELLCYCLEMODELODESOLVER_HPP_

#include <cxxtest/TestSuite.h>

#include "CellCycleModelOdeSolver.hpp"
#include "TysonNovakCellCycleModel.hpp"

#include "RungeKutta4IvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "CvodeAdaptor.hpp"
#include "OdeSystemInformation.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/**
 * Simple ODE system for use in the test suite. Defines the
 * IVP dy/dt = 1, y(0) = 0.
 */
class SimpleOde : public AbstractOdeSystem
{
public:
    SimpleOde() : AbstractOdeSystem(1) // 1 here is the number of variables
    {
        mpSystemInfo = OdeSystemInformation<SimpleOde>::Instance();
        SetStateVariables(GetInitialConditions());
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        rDY[0] = 1.0;
    }
};

template<>
void OdeSystemInformation<SimpleOde>::Initialise()
{
    this->mVariableNames.push_back("Variable 1");
    this->mVariableUnits.push_back("Units 1");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}

/**
 * Simple ODE system with stopping event for use in the test suite.
 * Defines the IVP dx/dt = y, dy/dt = -x, x(0) = 1, y(0) = 0.
 * Solutions to this system form circles about the origin. The
 * stopping event is x(0) < 0, which should first occur at time
 * t = pi/2.
 */
class OdeSecondOrderWithEvents : public AbstractOdeSystem
{
public :
    OdeSecondOrderWithEvents() : AbstractOdeSystem(2)
    {
        mpSystemInfo = OdeSystemInformation<OdeSecondOrderWithEvents>::Instance();
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        rDY[0] =  rY[1];
        rDY[1] = -rY[0];
    }

    bool CalculateStoppingEvent(double time, const std::vector<double>& rY)
    {
        return (rY[0] < 0);
    }
};

template<>
void OdeSystemInformation<OdeSecondOrderWithEvents>::Initialise()
{
    this->mVariableNames.push_back("Variable 1");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("Variable 2");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}

class TestCellCycleModelOdeSolver : public AbstractCellBasedTestSuite
{
public:

    void TestMethods() throw(Exception)
    {
        // Check we can create an instance
        boost::shared_ptr<CellCycleModelOdeSolver<TysonNovakCellCycleModel,RungeKutta4IvpOdeSolver> > p_solver
            = CellCycleModelOdeSolver<TysonNovakCellCycleModel,RungeKutta4IvpOdeSolver>::Instance();
        TS_ASSERT(p_solver.get() != NULL);

        // Check singleton-ness
        boost::shared_ptr<CellCycleModelOdeSolver<TysonNovakCellCycleModel,RungeKutta4IvpOdeSolver> > p_solver2
            = CellCycleModelOdeSolver<TysonNovakCellCycleModel,RungeKutta4IvpOdeSolver>::Instance();
        TS_ASSERT_EQUALS(p_solver, p_solver2);

        p_solver->Initialise();

        // Check the solver can be called for a simple ODE system
        SimpleOde ode;
        double last_time = 0.0;
        double current_time = 5;
        double dt = 1e-5;

        ode.SetStateVariables(ode.GetInitialConditions());
        p_solver->SolveAndUpdateStateVariable(&ode, last_time, current_time, dt);

        // No stopping event is specified, so check the solver did not stop
        TS_ASSERT_EQUALS(p_solver->StoppingEventOccurred(), false);

        // Check the solver can be called for another ODE system, this time with a stopping event
        OdeSecondOrderWithEvents ode_with_events;

        ode_with_events.SetStateVariables(ode_with_events.GetInitialConditions());
        p_solver->SolveAndUpdateStateVariable(&ode_with_events, last_time, current_time, dt);

        // Check the solver stopped at the correct time
        TS_ASSERT_EQUALS(p_solver->StoppingEventOccurred(), true);
        TS_ASSERT_DELTA(p_solver->GetStoppingTime(), M_PI_2, 1e-4);
    }

    void TestWithBackwardEulerIvpOdeSolver() throw(Exception)
    {
        // Check we can create an instance
        boost::shared_ptr<CellCycleModelOdeSolver<TysonNovakCellCycleModel, BackwardEulerIvpOdeSolver> > p_solver
            = CellCycleModelOdeSolver<TysonNovakCellCycleModel, BackwardEulerIvpOdeSolver>::Instance();
        TS_ASSERT(p_solver.get() != NULL);

        // Check singleton-ness
        boost::shared_ptr<CellCycleModelOdeSolver<TysonNovakCellCycleModel, BackwardEulerIvpOdeSolver> > p_solver2
            = CellCycleModelOdeSolver<TysonNovakCellCycleModel, BackwardEulerIvpOdeSolver>::Instance();
        TS_ASSERT_EQUALS(p_solver, p_solver2);

        TS_ASSERT_THROWS_THIS(p_solver->Initialise(), "SetSizeOfOdeSystem() must be called before calling Initialise()");
        p_solver->SetSizeOfOdeSystem(1);
        p_solver->Initialise();
       
        // Check the solver can be called for a simple ODE system
        SimpleOde ode;
        double last_time = 0.0;
        double current_time = 5;
        double dt = 1e-5;

        ode.SetStateVariables(ode.GetInitialConditions());
        p_solver->SolveAndUpdateStateVariable(&ode, last_time, current_time, dt);

        // No stopping event is specified, so check the solver did not stop
        TS_ASSERT_EQUALS(p_solver->StoppingEventOccurred(), false);

        p_solver->Reset();
        p_solver->SetSizeOfOdeSystem(2);
        p_solver->Initialise();

        // Check the solver can be called for another ODE system, this time with a stopping event
        OdeSecondOrderWithEvents ode_with_events;

        ode_with_events.SetStateVariables(ode_with_events.GetInitialConditions());
        p_solver->SolveAndUpdateStateVariable(&ode_with_events, last_time, current_time, dt);

        // Check the solver stopped at the correct time
        TS_ASSERT_EQUALS(p_solver->StoppingEventOccurred(), true);

        ///\todo The following line currently fails - work out why (see #1427)
//        TS_ASSERT_DELTA(p_solver->GetStoppingTime(), M_PI_2, 1e-4);
    }

    void TestWithCvodeAdaptor() throw(Exception)
    {
#ifdef CHASTE_CVODE
        // Check we can create an instance
        boost::shared_ptr<CellCycleModelOdeSolver<TysonNovakCellCycleModel,CvodeAdaptor> > p_solver = CellCycleModelOdeSolver<TysonNovakCellCycleModel, CvodeAdaptor>::Instance();
        TS_ASSERT(p_solver.get() != NULL);

        // Check singleton-ness
        boost::shared_ptr<CellCycleModelOdeSolver<TysonNovakCellCycleModel,CvodeAdaptor> > p_solver2 = CellCycleModelOdeSolver<TysonNovakCellCycleModel, CvodeAdaptor>::Instance();
        TS_ASSERT_EQUALS(p_solver, p_solver2);

        p_solver->Initialise();

        TS_ASSERT_THROWS_NOTHING(p_solver->CheckForStoppingEvents());
        TS_ASSERT_THROWS_NOTHING(p_solver->SetMaxSteps(1000));
        TS_ASSERT_THROWS_NOTHING(p_solver->SetTolerances(1e-5, 1e-5));
#endif // CHASTE_CVODE
    }

};


#endif /*TESTCELLCYCLEMODELODESOLVER_HPP_*/
