/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/ErrorPrinter.h>

int main() {
 return CxxTest::ErrorPrinter().run();
}
#include "odes/tests/TestAbstractIvpOdeSolver.hpp"

static TestAbstractIvpOdeSolver suite_TestAbstractIvpOdeSolver;

static CxxTest::List Tests_TestAbstractIvpOdeSolver = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestAbstractIvpOdeSolver( "odes/tests/TestAbstractIvpOdeSolver.hpp", 19, "TestAbstractIvpOdeSolver", suite_TestAbstractIvpOdeSolver, Tests_TestAbstractIvpOdeSolver );

static class TestDescription_TestAbstractIvpOdeSolver_testAddition : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractIvpOdeSolver_testAddition() : CxxTest::RealTestDescription( Tests_TestAbstractIvpOdeSolver, suiteDescription_TestAbstractIvpOdeSolver, 23, "testAddition" ) {}
 void runTest() { suite_TestAbstractIvpOdeSolver.testAddition(); }
} testDescription_TestAbstractIvpOdeSolver_testAddition;

static class TestDescription_TestAbstractIvpOdeSolver_testEulerSolver : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractIvpOdeSolver_testEulerSolver() : CxxTest::RealTestDescription( Tests_TestAbstractIvpOdeSolver, suiteDescription_TestAbstractIvpOdeSolver, 28, "testEulerSolver" ) {}
 void runTest() { suite_TestAbstractIvpOdeSolver.testEulerSolver(); }
} testDescription_TestAbstractIvpOdeSolver_testEulerSolver;

static class TestDescription_TestAbstractIvpOdeSolver_testAdamsBashforthSolver : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractIvpOdeSolver_testAdamsBashforthSolver() : CxxTest::RealTestDescription( Tests_TestAbstractIvpOdeSolver, suiteDescription_TestAbstractIvpOdeSolver, 49, "testAdamsBashforthSolver" ) {}
 void runTest() { suite_TestAbstractIvpOdeSolver.testAdamsBashforthSolver(); }
} testDescription_TestAbstractIvpOdeSolver_testAdamsBashforthSolver;

static class TestDescription_TestAbstractIvpOdeSolver_testRungeKutta2Solver : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractIvpOdeSolver_testRungeKutta2Solver() : CxxTest::RealTestDescription( Tests_TestAbstractIvpOdeSolver, suiteDescription_TestAbstractIvpOdeSolver, 72, "testRungeKutta2Solver" ) {}
 void runTest() { suite_TestAbstractIvpOdeSolver.testRungeKutta2Solver(); }
} testDescription_TestAbstractIvpOdeSolver_testRungeKutta2Solver;

static class TestDescription_TestAbstractIvpOdeSolver_testRungeKutta4Solver : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractIvpOdeSolver_testRungeKutta4Solver() : CxxTest::RealTestDescription( Tests_TestAbstractIvpOdeSolver, suiteDescription_TestAbstractIvpOdeSolver, 94, "testRungeKutta4Solver" ) {}
 void runTest() { suite_TestAbstractIvpOdeSolver.testRungeKutta4Solver(); }
} testDescription_TestAbstractIvpOdeSolver_testRungeKutta4Solver;

static class TestDescription_TestAbstractIvpOdeSolver_testLastTimeStep : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractIvpOdeSolver_testLastTimeStep() : CxxTest::RealTestDescription( Tests_TestAbstractIvpOdeSolver, suiteDescription_TestAbstractIvpOdeSolver, 116, "testLastTimeStep" ) {}
 void runTest() { suite_TestAbstractIvpOdeSolver.testLastTimeStep(); }
} testDescription_TestAbstractIvpOdeSolver_testLastTimeStep;

static class TestDescription_TestAbstractIvpOdeSolver_testGlobalError : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractIvpOdeSolver_testGlobalError() : CxxTest::RealTestDescription( Tests_TestAbstractIvpOdeSolver, suiteDescription_TestAbstractIvpOdeSolver, 140, "testGlobalError" ) {}
 void runTest() { suite_TestAbstractIvpOdeSolver.testGlobalError(); }
} testDescription_TestAbstractIvpOdeSolver_testGlobalError;

static class TestDescription_TestAbstractIvpOdeSolver_testGlobalErrorSystemOf2 : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractIvpOdeSolver_testGlobalErrorSystemOf2() : CxxTest::RealTestDescription( Tests_TestAbstractIvpOdeSolver, suiteDescription_TestAbstractIvpOdeSolver, 203, "testGlobalErrorSystemOf2" ) {}
 void runTest() { suite_TestAbstractIvpOdeSolver.testGlobalErrorSystemOf2(); }
} testDescription_TestAbstractIvpOdeSolver_testGlobalErrorSystemOf2;

static class TestDescription_TestAbstractIvpOdeSolver_testGlobalErrorSystemOf3 : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractIvpOdeSolver_testGlobalErrorSystemOf3() : CxxTest::RealTestDescription( Tests_TestAbstractIvpOdeSolver, suiteDescription_TestAbstractIvpOdeSolver, 290, "testGlobalErrorSystemOf3" ) {}
 void runTest() { suite_TestAbstractIvpOdeSolver.testGlobalErrorSystemOf3(); }
} testDescription_TestAbstractIvpOdeSolver_testGlobalErrorSystemOf3;

#include "odes/tests/TestAbstractOdeSystem.hpp"

static TestAbstractOdeSystem suite_TestAbstractOdeSystem;

static CxxTest::List Tests_TestAbstractOdeSystem = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestAbstractOdeSystem( "odes/tests/TestAbstractOdeSystem.hpp", 17, "TestAbstractOdeSystem", suite_TestAbstractOdeSystem, Tests_TestAbstractOdeSystem );

static class TestDescription_TestAbstractOdeSystem_testAddition : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractOdeSystem_testAddition() : CxxTest::RealTestDescription( Tests_TestAbstractOdeSystem, suiteDescription_TestAbstractOdeSystem, 21, "testAddition" ) {}
 void runTest() { suite_TestAbstractOdeSystem.testAddition(); }
} testDescription_TestAbstractOdeSystem_testAddition;

static class TestDescription_TestAbstractOdeSystem_TestOdeSystemOne : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractOdeSystem_TestOdeSystemOne() : CxxTest::RealTestDescription( Tests_TestAbstractOdeSystem, suiteDescription_TestAbstractOdeSystem, 26, "TestOdeSystemOne" ) {}
 void runTest() { suite_TestAbstractOdeSystem.TestOdeSystemOne(); }
} testDescription_TestAbstractOdeSystem_TestOdeSystemOne;

static class TestDescription_TestAbstractOdeSystem_TestOdeSystemTwo : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractOdeSystem_TestOdeSystemTwo() : CxxTest::RealTestDescription( Tests_TestAbstractOdeSystem, suiteDescription_TestAbstractOdeSystem, 40, "TestOdeSystemTwo" ) {}
 void runTest() { suite_TestAbstractOdeSystem.TestOdeSystemTwo(); }
} testDescription_TestAbstractOdeSystem_TestOdeSystemTwo;

static class TestDescription_TestAbstractOdeSystem_TestOdeSystemThree : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractOdeSystem_TestOdeSystemThree() : CxxTest::RealTestDescription( Tests_TestAbstractOdeSystem, suiteDescription_TestAbstractOdeSystem, 50, "TestOdeSystemThree" ) {}
 void runTest() { suite_TestAbstractOdeSystem.TestOdeSystemThree(); }
} testDescription_TestAbstractOdeSystem_TestOdeSystemThree;

#include <cxxtest/Root.cpp>
