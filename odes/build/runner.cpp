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
#include "odes/tests/TestAbstractOdeSystem.hpp"

static TestAbstractOdeSystem suite_TestAbstractOdeSystem;

static CxxTest::List Tests_TestAbstractOdeSystem = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestAbstractOdeSystem( "odes/tests/TestAbstractOdeSystem.hpp", 19, "TestAbstractOdeSystem", suite_TestAbstractOdeSystem, Tests_TestAbstractOdeSystem );

static class TestDescription_TestAbstractOdeSystem_testAddition : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractOdeSystem_testAddition() : CxxTest::RealTestDescription( Tests_TestAbstractOdeSystem, suiteDescription_TestAbstractOdeSystem, 32, "testAddition" ) {}
 void runTest() { suite_TestAbstractOdeSystem.testAddition(); }
} testDescription_TestAbstractOdeSystem_testAddition;

static class TestDescription_TestAbstractOdeSystem_testPetsc : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractOdeSystem_testPetsc() : CxxTest::RealTestDescription( Tests_TestAbstractOdeSystem, suiteDescription_TestAbstractOdeSystem, 39, "testPetsc" ) {}
 void runTest() { suite_TestAbstractOdeSystem.testPetsc(); }
} testDescription_TestAbstractOdeSystem_testPetsc;

static class TestDescription_TestAbstractOdeSystem_TestMyAbstractOdeSystem : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractOdeSystem_TestMyAbstractOdeSystem() : CxxTest::RealTestDescription( Tests_TestAbstractOdeSystem, suiteDescription_TestAbstractOdeSystem, 69, "TestMyAbstractOdeSystem" ) {}
 void runTest() { suite_TestAbstractOdeSystem.TestMyAbstractOdeSystem(); }
} testDescription_TestAbstractOdeSystem_TestMyAbstractOdeSystem;

static class TestDescription_TestAbstractOdeSystem_TestOdeOne : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractOdeSystem_TestOdeOne() : CxxTest::RealTestDescription( Tests_TestAbstractOdeSystem, suiteDescription_TestAbstractOdeSystem, 94, "TestOdeOne" ) {}
 void runTest() { suite_TestAbstractOdeSystem.TestOdeOne(); }
} testDescription_TestAbstractOdeSystem_TestOdeOne;

static class TestDescription_TestAbstractOdeSystem_TestOdeTwo : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractOdeSystem_TestOdeTwo() : CxxTest::RealTestDescription( Tests_TestAbstractOdeSystem, suiteDescription_TestAbstractOdeSystem, 114, "TestOdeTwo" ) {}
 void runTest() { suite_TestAbstractOdeSystem.TestOdeTwo(); }
} testDescription_TestAbstractOdeSystem_TestOdeTwo;

static class TestDescription_TestAbstractOdeSystem_TestOdeThree : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAbstractOdeSystem_TestOdeThree() : CxxTest::RealTestDescription( Tests_TestAbstractOdeSystem, suiteDescription_TestAbstractOdeSystem, 129, "TestOdeThree" ) {}
 void runTest() { suite_TestAbstractOdeSystem.TestOdeThree(); }
} testDescription_TestAbstractOdeSystem_TestOdeThree;

#include <cxxtest/Root.cpp>
