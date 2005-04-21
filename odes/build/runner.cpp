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
CxxTest::StaticSuiteDescription suiteDescription_TestAbstractIvpOdeSolver( "odes/tests/TestAbstractIvpOdeSolver.hpp", 13, "TestAbstractIvpOdeSolver", suite_TestAbstractIvpOdeSolver, Tests_TestAbstractIvpOdeSolver );

static class TestDescription_TestAbstractIvpOdeSolver_testAddition : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractIvpOdeSolver_testAddition() : CxxTest::RealTestDescription( Tests_TestAbstractIvpOdeSolver, suiteDescription_TestAbstractIvpOdeSolver, 17, "testAddition" ) {}
 void runTest() { suite_TestAbstractIvpOdeSolver.testAddition(); }
} testDescription_TestAbstractIvpOdeSolver_testAddition;

static class TestDescription_TestAbstractIvpOdeSolver_testEulerSolver : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractIvpOdeSolver_testEulerSolver() : CxxTest::RealTestDescription( Tests_TestAbstractIvpOdeSolver, suiteDescription_TestAbstractIvpOdeSolver, 22, "testEulerSolver" ) {}
 void runTest() { suite_TestAbstractIvpOdeSolver.testEulerSolver(); }
} testDescription_TestAbstractIvpOdeSolver_testEulerSolver;

static class TestDescription_TestAbstractIvpOdeSolver_testRK2Solver : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractIvpOdeSolver_testRK2Solver() : CxxTest::RealTestDescription( Tests_TestAbstractIvpOdeSolver, suiteDescription_TestAbstractIvpOdeSolver, 43, "testRK2Solver" ) {}
 void runTest() { suite_TestAbstractIvpOdeSolver.testRK2Solver(); }
} testDescription_TestAbstractIvpOdeSolver_testRK2Solver;

static class TestDescription_TestAbstractIvpOdeSolver_testLastTimeStep : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractIvpOdeSolver_testLastTimeStep() : CxxTest::RealTestDescription( Tests_TestAbstractIvpOdeSolver, suiteDescription_TestAbstractIvpOdeSolver, 67, "testLastTimeStep" ) {}
 void runTest() { suite_TestAbstractIvpOdeSolver.testLastTimeStep(); }
} testDescription_TestAbstractIvpOdeSolver_testLastTimeStep;

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
