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
CxxTest::StaticSuiteDescription suiteDescription_TestNonnecroticTumourNonlinearEllipticEquation( "cancer/tests/TestNonnecroticTumourNonlinearEllipticEquation.hpp", 11, "TestNonnecroticTumourNonlinearEllipticEquation", suite_TestNonnecroticTumourNonlinearEllipticEquation, Tests_TestNonnecroticTumourNonlinearEllipticEquation );

static class TestDescription_TestNonnecroticTumourNonlinearEllipticEquation_testNonnecroticTumourC : public CxxTest::RealTestDescription {
public:
 TestDescription_TestNonnecroticTumourNonlinearEllipticEquation_testNonnecroticTumourC() : CxxTest::RealTestDescription( Tests_TestNonnecroticTumourNonlinearEllipticEquation, suiteDescription_TestNonnecroticTumourNonlinearEllipticEquation, 14, "testNonnecroticTumourC" ) {}
 void runTest() { suite_TestNonnecroticTumourNonlinearEllipticEquation.testNonnecroticTumourC(); }
} testDescription_TestNonnecroticTumourNonlinearEllipticEquation_testNonnecroticTumourC;

#include "cancer/tests/TestPracticalOne.hpp"

static TestPracticalOne suite_TestPracticalOne;

static CxxTest::List Tests_TestPracticalOne = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestPracticalOne( "cancer/tests/TestPracticalOne.hpp", 30, "TestPracticalOne", suite_TestPracticalOne, Tests_TestPracticalOne );

static class TestDescription_TestPracticalOne_testQuestion1 : public CxxTest::RealTestDescription {
public:
 TestDescription_TestPracticalOne_testQuestion1() : CxxTest::RealTestDescription( Tests_TestPracticalOne, suiteDescription_TestPracticalOne, 34, "testQuestion1" ) {}
 void runTest() { suite_TestPracticalOne.testQuestion1(); }
} testDescription_TestPracticalOne_testQuestion1;

static class TestDescription_TestPracticalOne_testQuestion1Nonlinear : public CxxTest::RealTestDescription {
public:
 TestDescription_TestPracticalOne_testQuestion1Nonlinear() : CxxTest::RealTestDescription( Tests_TestPracticalOne, suiteDescription_TestPracticalOne, 81, "testQuestion1Nonlinear" ) {}
 void runTest() { suite_TestPracticalOne.testQuestion1Nonlinear(); }
} testDescription_TestPracticalOne_testQuestion1Nonlinear;

static class TestDescription_TestPracticalOne_testQuestion2 : public CxxTest::RealTestDescription {
public:
 TestDescription_TestPracticalOne_testQuestion2() : CxxTest::RealTestDescription( Tests_TestPracticalOne, suiteDescription_TestPracticalOne, 183, "testQuestion2" ) {}
 void runTest() { suite_TestPracticalOne.testQuestion2(); }
} testDescription_TestPracticalOne_testQuestion2;

static class TestDescription_TestPracticalOne_testQuestions3to5 : public CxxTest::RealTestDescription {
public:
 TestDescription_TestPracticalOne_testQuestions3to5() : CxxTest::RealTestDescription( Tests_TestPracticalOne, suiteDescription_TestPracticalOne, 274, "testQuestions3to5" ) {}
 void runTest() { suite_TestPracticalOne.testQuestions3to5(); }
} testDescription_TestPracticalOne_testQuestions3to5;

#include <cxxtest/Root.cpp>
