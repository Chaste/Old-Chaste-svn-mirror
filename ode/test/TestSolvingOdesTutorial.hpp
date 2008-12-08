/*

Copyright (C) University of Oxford, 2008

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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */
#ifndef TESTSOLVINGODESTUTORIAL_HPP_
#define TESTSOLVINGODESTUTORIAL_HPP_
/*
 * = Introduction =
 *
 * In this tutorial we show how Chaste can be solve an ODE.
 *
 * EMPTYLINE
 *
 * The following header files need to be included.
 * First we include the header needed to define this class as a test suite. */
#include <cxxtest/TestSuite.h>
/*
 * We will use a simple forward euler solver to solve the ODE, so the following
 * needs to be included
 */
#include "EulerIvpOdeSolver.hpp"
/* All the ODE solvers take in a concrete ODE system class, which is user-defined
 * and must inherit from the following class, which defines an ODE interface.
 */
#include "AbstractOdeSystem.hpp"
/* In order to convenient define useful information about the ODE system, such
 * as the names and units of variables, and suggested initial conditions, we
 * need the following header.
 */
#include "OdeSystemInformation.hpp"
/*
 * EMPTYLINE
 *
 * = Defining the ODE classes =
 *
 * Let us solve the ODE dy/dt = y^2^+t^2^, with y(0) = 1. To do so, we have to define
 * our own ODE class, inheriting from {{{AbstractOdeSystem}}}, which implements that
 * {{{EvaluateYDerivatives()}}} method.
 */
class MyOde : public AbstractOdeSystem
{
public:
/* The constructor does very little.
 * It calls the base constructor, passing the number of state variables in the
 * ODE system (here, 1, i.e. y is a 1d vector).
 * It also sets the object to use to retrieve system information (see later).
 */
    MyOde() : AbstractOdeSystem(1)
    {
        mpSystemInfo = OdeSystemInformation<MyOde>::Instance();
    }

/* The ODE solvers will repeatedly call a method called EvaluateYDerivatives(), which needs
 * to be implemented in this concrete class. This takes in the time, a {{{std::vector}}} of
 * y values (in this, of size 1), and a reference to a {{{std::vector}}} in which the
 * derivative(s) should be filled in by the method..
 */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY,
                              std::vector<double> &rDY)
    {
        /*..so we set {{{rDY[0]}}} to be y^2^ + t^2^. */
        rDY[0] = rY[0]*rY[0] + time*time;
    }
};

/* The following ''template specialisation'' defines the information for this
 * ODE system.  Note that we use the ODE system class as a template parameter.
 */
template<>
void OdeSystemInformation<MyOde>::Initialise(void)
{
    this->mVariableNames.push_back("y");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);
    
    this->mInitialised = true;
}

/* That would be all that is needed for this class to solve the ODE. However, rather
 * than solving up to a fixed time, suppose we wanted to solve until some function
 * of y (and t) reached a certain value, e.g. let's say we wanted to solve the ODE until
 * y reached 2.5. To do this, we have to define a stopping event, by implementing
 * the method {{{CalculateStoppingEvent()}}} in {{{AbstractOdeSystem}}}. For this, let us
 * define a new class, inheriting from the above class (i.e. representing the same ODE)
 * but with a stopping event defined.
 */
class MyOdeWithStoppingEvent : public MyOde
{
public:
    /* All we have to do is implement the following function. This is defined in
     * the base class ({{{AbstractOdeSystem}}}), where it always returns false, and here we override it
     * to return true if y>=2.5
     */
    bool CalculateStoppingEvent(double time, const std::vector<double> &rY)
    {
        return (rY[0]>=2.5);
    }
};

/* (Ignore this class until solving with state variables is discussed)
 *
 * Another class which sets up a state variable. Note that this is done in the
 * constructor, and the {{{EvaluateYDerivatives}}} is identical to before */
class MyOdeUsingStateVariables : public AbstractOdeSystem
{
public:
    MyOdeUsingStateVariables() : AbstractOdeSystem(1)
    {
        mpSystemInfo = OdeSystemInformation<MyOdeUsingStateVariables>::Instance();
        mStateVariables.push_back(1.0);
    }

    void EvaluateYDerivatives(double time, const std::vector<double> &rY,
                              std::vector<double> &rDY)
    {
        rDY[0] = rY[0]*rY[0] + time*time;
    }
};

/* Again we need to define the ODE system information.
 */
template<>
void OdeSystemInformation<MyOdeUsingStateVariables>::Initialise(void)
{
    this->mVariableNames.push_back("y");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);
    
    this->mInitialised = true;
}

/* This class is another simple ODE class, just as an example of how a 2d ODE is solved. Here
 * we solve the ODE dy,,1,,/dt = y,,2,,, dy,,2,,/dt = (y,,1,,)^2^ (which represents the second-order ODE d^2^y/dt^2^ = y^2^
 */
class My2dOde : public AbstractOdeSystem
{
public:
    My2dOde() : AbstractOdeSystem(2)
    {
        mpSystemInfo = OdeSystemInformation<My2dOde>::Instance();
    }

    void EvaluateYDerivatives(double time, const std::vector<double> &rY,
                              std::vector<double> &rDY)
    {
        rDY[0] = rY[1];
        rDY[1] = rY[0]*rY[0];
    }
};

/* Again we need to define the ODE system information.
 */
template<>
void OdeSystemInformation<My2dOde>::Initialise(void)
{
    this->mVariableNames.push_back("y");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);
    
    this->mVariableNames.push_back("dy/dt");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);
    
    this->mInitialised = true;
}

/*
 * EMPTYLINE
 *
 * = The Tests =
 *
 * EMPTYLINE
 *
 * == Standard ODE Solving ==
 *
 * Now we can define the test, where the ODEs are solved. */
class TestSolvingOdesTutorial: public CxxTest::TestSuite
{
public:
    void TestSolvingOdes() throw(Exception)
    {
        /* First, create an instance of the ODE class to be solved. */
        MyOde my_ode;
        /* Next, create a solver. */
        EulerIvpOdeSolver euler_solver;
        /* We will need to provide an initial condition, which needs to
         * be a {{{std::vector}}}.*/
        std::vector<double> initial_condition;
        initial_condition.push_back(1.0);
        /* Then, just call Solve(), passing in a pointer to the ODE, the
         * initial condition, the start time, end time, the solving timestep,
         * and sampling timestep (how often we want the returned solution).
         * Here we solve from 0 to 1, with a timestep of 0.01 but a sampling
         * timestep of 0.1. The return value is an object of type {{{OdeSolution}}}
         * (which is basically just a list of times and solutions).
         */
        OdeSolution solutions = euler_solver.Solve(&my_ode, initial_condition, 0, 1, 0.01, 0.1);
        /* Let's look at the results, which can be obtained from the {{{OdeSolutions}}}
         * object using the methods {{{rGetTimes()}}} and {{{rGetSolutions()}}}, which
         * return a {{{std::vector}}} and a {{{std::vector}}} of {{{std::vector}}}s
         * respectively. */
        for(unsigned i=0; i<solutions.rGetTimes().size(); i++)
        {
            /* the {{{[0]}}} here is because getting the zeroth component of y (a 1-dimensional vector) */
            std::cout << solutions.rGetTimes()[i] << " " << solutions.rGetSolutions()[i][0] << "\n";
        }

        /* Alternatively, we can print the solution directly to a file, using the {{{WriteToFile}}}
         * method on the {{{OdeSolution}}} class. (To do this, we need to provide the
         * ODE system and the units time is in, since a header line is written (this is because
         * the ODE system can contain the names and units of the variables, although our's
         * does not have that defined)) */
        solutions.WriteToFile("SolvingOdesTutorial", "ode1.txt", &my_ode, "sec");

        /* We can see from the printed out results that y goes above 2.5 somewhere just
         * before 0.6. To solve only up until y=2.5, we can solve the ODE that has the
         * stopping event defined, using the same solver as before. */
        MyOdeWithStoppingEvent my_ode_stopping;

        /* '''Note:''' ''when a {{{std::vector}}} is passed in as an initial condition
         * to a {{{Solve}}} call, it gets updated as the solve takes place''. Therefore, if
         * we want to use the same initial condition again, we have to reset it back to 1.0 */
        initial_condition[0] = 1.0;
        solutions = euler_solver.Solve(&my_ode_stopping, initial_condition, 0, 1, 0.01, 0.1);
        /* We can check with the solver that it stopped because of the stopping event, rather than because
         * it reached to end time. */
        assert(euler_solver.StoppingEventOccured()==true);
        /* Finally, let's print the time of the stopping event (to the nearest dt or so). */
        std::cout << "Stopping event occured at t="<<solutions.rGetTimes().back()<<"\n";
    }

    /*
     * EMPTYLINE
     *
     * == ODE Solving Using the State Variable ==
     *
     * In this second test, we show how to do an alternative version of ODE solving, which
     * does not involve passing in initial conditions and returning a {{{OdeSolution}}}.
     * The {{{AbstractOdeSystem}}} has a variable called the ''state variable'', which can
     * be used to hold the solution, and will be updated if a particular version of Solve
     * is called. This can be useful for embedding ODE models in a bigger system, since
     * the ODE models will then always contain their current solution.
     */
    void TestOdeSolvingUsingStateVariable()
    {
        /* Define an instance of the ODE. See the class definition above.
         * Note that this ODE has a variable called {{{mStateVariables}}}, which has
         * been set to be a vector of size one, containing the value 1.0. */
        MyOdeUsingStateVariables my_ode_using_state_vars;

        /* To solve updating the state variable, just call appropriate method with
         * a chosen solver. Note that no initial condition is required, no
         * {{{OdeSolution}}} is returned, and no sampling timestep is given. */
        EulerIvpOdeSolver euler_solver;
        euler_solver.SolveAndUpdateStateVariable(&my_ode_using_state_vars, 0.0, 1.0, 0.01);

        /* To see what the solution was at the end, let's print out the state variable. */
        std::cout << "Solution at end time = " << my_ode_using_state_vars.rGetStateVariables()[0] << "\n";
    }

    /*
     * EMPTYLINE
     *
     * == Solving n-dimensional ODEs ==
     *
     * Finally, here's a simple test showing how to solve a 2d ODE using the first method.
     * All that is different is the initial condition has be a 2d vector, and returned
     * solution is 2d at every timestep.
     */
    void TestWith2dOde()
    {
        My2dOde my_2d_ode;
        EulerIvpOdeSolver euler_solver;

        /* Define a 2d initial condition. */
        std::vector<double> initial_condition;
        initial_condition.push_back(1.0);
        initial_condition.push_back(0.0);

        /* Solve, and print the solution as [time, y1, y2]. */
        OdeSolution solutions = euler_solver.Solve(&my_2d_ode, initial_condition, 0, 1, 0.01, 0.1);
        for(unsigned i=0; i<solutions.rGetTimes().size(); i++)
        {
            std::cout << solutions.rGetTimes()[i] << " "
                      << solutions.rGetSolutions()[i][0] << " "
                      << solutions.rGetSolutions()[i][1] << "\n";
        }
    }
};
#endif /*TESTSOLVINGODESTUTORIAL_HPP_*/
