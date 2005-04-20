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
#include "heart/tests/TestLR91.hpp"

static TestLR91 suite_TestLR91;

static CxxTest::List Tests_TestLR91 = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestLR91( "heart/tests/TestLR91.hpp", 19, "TestLR91", suite_TestLR91, Tests_TestLR91 );

static class TestDescription_TestLR91_testCurrentsMagnitude : public CxxTest::RealTestDescription {
public:
 TestDescription_TestLR91_testCurrentsMagnitude() : CxxTest::RealTestDescription( Tests_TestLR91, suiteDescription_TestLR91, 24, "testCurrentsMagnitude" ) {}
 void runTest() { suite_TestLR91.testCurrentsMagnitude(); }
} testDescription_TestLR91_testCurrentsMagnitude;

#include "heart/tests/TestIonicCurrent.hpp"

static TestIonicCurrent suite_TestIonicCurrent;

static CxxTest::List Tests_TestIonicCurrent = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestIonicCurrent( "heart/tests/TestIonicCurrent.hpp", 6, "TestIonicCurrent", suite_TestIonicCurrent, Tests_TestIonicCurrent );

static class TestDescription_TestIonicCurrent_testIonicCurrent : public CxxTest::RealTestDescription {
public:
 TestDescription_TestIonicCurrent_testIonicCurrent() : CxxTest::RealTestDescription( Tests_TestIonicCurrent, suiteDescription_TestIonicCurrent, 11, "testIonicCurrent" ) {}
 void runTest() { suite_TestIonicCurrent.testIonicCurrent(); }
} testDescription_TestIonicCurrent_testIonicCurrent;

#include <cxxtest/Root.cpp>
