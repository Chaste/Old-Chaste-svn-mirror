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
#include "cancer/tests/TestNonnecroticTumourNonlinearEllipticEquation.hpp"

static TestNonnecroticTumourNonlinearEllipticEquation suite_TestNonnecroticTumourNonlinearEllipticEquation;

static CxxTest::List Tests_TestNonnecroticTumourNonlinearEllipticEquation = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestNonnecroticTumourNonlinearEllipticEquation( "cancer/tests/TestNonnecroticTumourNonlinearEllipticEquation.hpp", 12, "TestNonnecroticTumourNonlinearEllipticEquation", suite_TestNonnecroticTumourNonlinearEllipticEquation, Tests_TestNonnecroticTumourNonlinearEllipticEquation );

static class TestDescription_TestNonnecroticTumourNonlinearEllipticEquation_testNonnecroticTumourC : public CxxTest::RealTestDescription {
public:
 TestDescription_TestNonnecroticTumourNonlinearEllipticEquation_testNonnecroticTumourC() : CxxTest::RealTestDescription( Tests_TestNonnecroticTumourNonlinearEllipticEquation, suiteDescription_TestNonnecroticTumourNonlinearEllipticEquation, 15, "testNonnecroticTumourC" ) {}
 void runTest() { suite_TestNonnecroticTumourNonlinearEllipticEquation.testNonnecroticTumourC(); }
} testDescription_TestNonnecroticTumourNonlinearEllipticEquation_testNonnecroticTumourC;

#include "cancer/tests/TestPracticalOne.hpp"

static TestPracticalOne suite_TestPracticalOne;

static CxxTest::List Tests_TestPracticalOne = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestPracticalOne( "cancer/tests/TestPracticalOne.hpp", 31, "TestPracticalOne", suite_TestPracticalOne, Tests_TestPracticalOne );

static class TestDescription_TestPracticalOne_testPrac1Question1 : public CxxTest::RealTestDescription {
public:
 TestDescription_TestPracticalOne_testPrac1Question1() : CxxTest::RealTestDescription( Tests_TestPracticalOne, suiteDescription_TestPracticalOne, 35, "testPrac1Question1" ) {}
 void runTest() { suite_TestPracticalOne.testPrac1Question1(); }
} testDescription_TestPracticalOne_testPrac1Question1;

static class TestDescription_TestPracticalOne_testPrac1Question1Nonlinear : public CxxTest::RealTestDescription {
public:
 TestDescription_TestPracticalOne_testPrac1Question1Nonlinear() : CxxTest::RealTestDescription( Tests_TestPracticalOne, suiteDescription_TestPracticalOne, 82, "testPrac1Question1Nonlinear" ) {}
 void runTest() { suite_TestPracticalOne.testPrac1Question1Nonlinear(); }
} testDescription_TestPracticalOne_testPrac1Question1Nonlinear;

static class TestDescription_TestPracticalOne_testPrac1Question2 : public CxxTest::RealTestDescription {
public:
 TestDescription_TestPracticalOne_testPrac1Question2() : CxxTest::RealTestDescription( Tests_TestPracticalOne, suiteDescription_TestPracticalOne, 174, "testPrac1Question2" ) {}
 void runTest() { suite_TestPracticalOne.testPrac1Question2(); }
} testDescription_TestPracticalOne_testPrac1Question2;

static class TestDescription_TestPracticalOne_testPrac1Questions3to5 : public CxxTest::RealTestDescription {
public:
 TestDescription_TestPracticalOne_testPrac1Questions3to5() : CxxTest::RealTestDescription( Tests_TestPracticalOne, suiteDescription_TestPracticalOne, 265, "testPrac1Questions3to5" ) {}
 void runTest() { suite_TestPracticalOne.testPrac1Questions3to5(); }
} testDescription_TestPracticalOne_testPrac1Questions3to5;

static class TestDescription_TestPracticalOne_testPrac1Question6 : public CxxTest::RealTestDescription {
public:
 TestDescription_TestPracticalOne_testPrac1Question6() : CxxTest::RealTestDescription( Tests_TestPracticalOne, suiteDescription_TestPracticalOne, 513, "testPrac1Question6" ) {}
 void runTest() { suite_TestPracticalOne.testPrac1Question6(); }
} testDescription_TestPracticalOne_testPrac1Question6;

#include "cancer/tests/TestPracticalTwo.hpp"

static TestPracticalTwo suite_TestPracticalTwo;

static CxxTest::List Tests_TestPracticalTwo = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestPracticalTwo( "cancer/tests/TestPracticalTwo.hpp", 27, "TestPracticalTwo", suite_TestPracticalTwo, Tests_TestPracticalTwo );

static class TestDescription_TestPracticalTwo_testPrac2Question1 : public CxxTest::RealTestDescription {
public:
 TestDescription_TestPracticalTwo_testPrac2Question1() : CxxTest::RealTestDescription( Tests_TestPracticalTwo, suiteDescription_TestPracticalTwo, 32, "testPrac2Question1" ) {}
 void runTest() { suite_TestPracticalTwo.testPrac2Question1(); }
} testDescription_TestPracticalTwo_testPrac2Question1;

#include <cxxtest/Root.cpp>
