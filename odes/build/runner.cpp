/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#define _CXXTEST_HAVE_EH
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

static class TestDescription_TestAbstractIvpOdeSolver_testconstructOdeSolver : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractIvpOdeSolver_testconstructOdeSolver() : CxxTest::RealTestDescription( Tests_TestAbstractIvpOdeSolver, suiteDescription_TestAbstractIvpOdeSolver, 27, "testconstructOdeSolver" ) {}
 void runTest() { suite_TestAbstractIvpOdeSolver.testconstructOdeSolver(); }
} testDescription_TestAbstractIvpOdeSolver_testconstructOdeSolver;

static class TestDescription_TestAbstractIvpOdeSolver_testIvpOdeSolver : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractIvpOdeSolver_testIvpOdeSolver() : CxxTest::RealTestDescription( Tests_TestAbstractIvpOdeSolver, suiteDescription_TestAbstractIvpOdeSolver, 33, "testIvpOdeSolver" ) {}
 void runTest() { suite_TestAbstractIvpOdeSolver.testIvpOdeSolver(); }
} testDescription_TestAbstractIvpOdeSolver_testIvpOdeSolver;

#include "odes/tests/TestAbstractOdeSystem.hpp"

static TestAbstractOdeSystem suite_TestAbstractOdeSystem;

static CxxTest::List Tests_TestAbstractOdeSystem = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestAbstractOdeSystem( "odes/tests/TestAbstractOdeSystem.hpp", 17, "TestAbstractOdeSystem", suite_TestAbstractOdeSystem, Tests_TestAbstractOdeSystem );

static class TestDescription_TestAbstractOdeSystem_testAddition : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractOdeSystem_testAddition() : CxxTest::RealTestDescription( Tests_TestAbstractOdeSystem, suiteDescription_TestAbstractOdeSystem, 21, "testAddition" ) {}
 void runTest() { suite_TestAbstractOdeSystem.testAddition(); }
} testDescription_TestAbstractOdeSystem_testAddition;

static class TestDescription_TestAbstractOdeSystem_TestOdeOne : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractOdeSystem_TestOdeOne() : CxxTest::RealTestDescription( Tests_TestAbstractOdeSystem, suiteDescription_TestAbstractOdeSystem, 29, "TestOdeOne" ) {}
 void runTest() { suite_TestAbstractOdeSystem.TestOdeOne(); }
} testDescription_TestAbstractOdeSystem_TestOdeOne;

static class TestDescription_TestAbstractOdeSystem_TestOdeTwo : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractOdeSystem_TestOdeTwo() : CxxTest::RealTestDescription( Tests_TestAbstractOdeSystem, suiteDescription_TestAbstractOdeSystem, 45, "TestOdeTwo" ) {}
 void runTest() { suite_TestAbstractOdeSystem.TestOdeTwo(); }
} testDescription_TestAbstractOdeSystem_TestOdeTwo;

static class TestDescription_TestAbstractOdeSystem_TestOdeThree : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractOdeSystem_TestOdeThree() : CxxTest::RealTestDescription( Tests_TestAbstractOdeSystem, suiteDescription_TestAbstractOdeSystem, 60, "TestOdeThree" ) {}
 void runTest() { suite_TestAbstractOdeSystem.TestOdeThree(); }
} testDescription_TestAbstractOdeSystem_TestOdeThree;

#include <cxxtest/Root.cpp>
