/*
 * 
 *  Chaste Pde Tutorial - this page get automatically changed to a wiki page
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
 * In this tutorial we show how Chaste can be solve an ODE .
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
/*
 * EMPTYLINE
 * 
 * Let us solve the ODE dy/dt = y^2^+t^2^, with y(0) = 1. To do so, we have to define
 * our own ODE class, inheriting from {{{AbstractOdeSystem}}}, which implements that
 * {{{EvaluateYDerivatives()}}} method.
 */
class MyOde : public AbstractOdeSystem
{
public:
/* The constructor does nothing, except calling the base constructor, with the number of 
 * state variables in the ODE system (here, 1, i.e. y is a 1d vector).
 */
    MyOde() : AbstractOdeSystem(1)
    {
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

/* Now we can define the test, where the ODEs are solved. */
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
         * method on the {{{OdeSolution}}} class. To do this, we need to provide the
         * ODE system and the units time is in, since a header line is written. '''TODO'''*/
        //solutions.WriteToFile("SolvingOdesTutorial", "ode1.txt", &my_ode, "sec");
        
        /* We can see from the printed out results that y goes above 2.5 somewhere just 
         * before 0.6. To solve only up until y=2.5, we can solve the ODE that has the
         * stopping event defined, using the same solver as before. */
        MyOdeWithStoppingEvent my_ode_stopping;
        /* '''TODO:''' shouldn't have to redefine the initial_conditon. */
        initial_condition[0] = 1.0;
        solutions = euler_solver.Solve(&my_ode_stopping, initial_condition, 0, 1, 0.01, 0.1);
        /* We can check with the solver that it stopped because of the stopping event, rather than because
         * it reached to end time. */
        assert(euler_solver.StoppingEventOccured()==true);
        /* Finally, let's print the time of the stopping event (to the nearest dt or so). */
        std::cout << "Stopping event occured at t="<<solutions.rGetTimes().back()<<"\n";
    }
};
#endif /*TESTSOLVINGODESTUTORIAL_HPP_*/
