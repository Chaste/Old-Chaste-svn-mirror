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
#include "pdes/paralleltests/TestParallelLinearSystem.hpp"

static TestParallelLinearSystem suite_TestParallelLinearSystem;

static CxxTest::List Tests_TestParallelLinearSystem = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestParallelLinearSystem( "pdes/paralleltests/TestParallelLinearSystem.hpp", 9, "TestParallelLinearSystem", suite_TestParallelLinearSystem, Tests_TestParallelLinearSystem );

static class TestDescription_TestParallelLinearSystem_testLinearSystem1 : public CxxTest::RealTestDescription {
public:
 TestDescription_TestParallelLinearSystem_testLinearSystem1() : CxxTest::RealTestDescription( Tests_TestParallelLinearSystem, suiteDescription_TestParallelLinearSystem, 12, "testLinearSystem1" ) {}
 void runTest() { suite_TestParallelLinearSystem.testLinearSystem1(); }
} testDescription_TestParallelLinearSystem_testLinearSystem1;

static class TestDescription_TestParallelLinearSystem_testLinearSystem2 : public CxxTest::RealTestDescription {
public:
 TestDescription_TestParallelLinearSystem_testLinearSystem2() : CxxTest::RealTestDescription( Tests_TestParallelLinearSystem, suiteDescription_TestParallelLinearSystem, 48, "testLinearSystem2" ) {}
 void runTest() { suite_TestParallelLinearSystem.testLinearSystem2(); }
} testDescription_TestParallelLinearSystem_testLinearSystem2;

#include "pdes/paralleltests/TestPetSCSetup.hpp"

static TestPetSCSetup suite_TestPetSCSetup;

static CxxTest::List Tests_TestPetSCSetup = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestPetSCSetup( "pdes/paralleltests/TestPetSCSetup.hpp", 38, "TestPetSCSetup", suite_TestPetSCSetup, Tests_TestPetSCSetup );

static class TestDescription_TestPetSCSetup_testPetscIsThere : public CxxTest::RealTestDescription {
public:
 TestDescription_TestPetSCSetup_testPetscIsThere() : CxxTest::RealTestDescription( Tests_TestPetSCSetup, suiteDescription_TestPetSCSetup, 41, "testPetscIsThere" ) {}
 void runTest() { suite_TestPetSCSetup.testPetscIsThere(); }
} testDescription_TestPetSCSetup_testPetscIsThere;

#include <cxxtest/Root.cpp>
