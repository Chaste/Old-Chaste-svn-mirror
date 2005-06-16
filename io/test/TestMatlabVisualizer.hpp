// TestMatlabVisualizer.hpp

/**
 * Test suite for the MatlabVisualizer class
 * 
 */

#ifndef _TESTMATLABVISUALIZER_HPP_
#define _TESTMATLABVISUALIZER_HPP_

#include <cxxtest/TestSuite.h>
#include "../AbstractVisualizer.hpp"
#include "../MatlabVisualizer.cpp"

static AbstractVisualizer<1> *spViewer;
static AbstractVisualizer<2> *sp2DViewer;

class TestMatlabVisualizer : public CxxTest::TestSuite
{
	public:
	
	/**
	 * Check that input files are opened correctly.
	 * 
	 */

//	void Test1DVisualizationFileBuild(void)
//	{
//		TS_ASSERT_THROWS_NOTHING(
//		                  spViewer=new MatlabVisualizer<1>(
//		                  "pdes/tests/meshdata/trivial_1d_mesh"));
//		spViewer->CreateNodesFileForVisualization();		                  
//		delete spViewer;
//	}
//	
//	void Test2DVisualizationFileBuild(void)
//	{
//		TS_ASSERT_THROWS_NOTHING(
//		                  sp2DViewer=new MatlabVisualizer<2>(
//		                  "pdes/tests/meshdata/disk_522_elements"));
//		TS_ASSERT_THROWS_NOTHING(
//					sp2DViewer->CreateNodesFileForVisualization());
//		delete sp2DViewer;
//	}
//	
//	/**
//	 * Check whether it can read time file correctly.	 * 
//	 */	
//	void TestVisualizationTimeFileRead(void)
//	{
//		TS_ASSERT_THROWS_NOTHING(
//		                  spViewer=new MatlabVisualizer<1>(
//		                  "pdes/tests/meshdata/trivial_1d_mesh"));
//		TS_ASSERT_THROWS_ANYTHING(
//					spViewer->CreateOutputFileForVisualization());	
//		
//		TS_ASSERT_THROWS_NOTHING(
//		                  sp2DViewer=new MatlabVisualizer<2>(
//		                  "pdes/tests/meshdata/disk_522_elements"));
//		TS_ASSERT_THROWS_ANYTHING(
//						sp2DViewer->CreateOutputFileForVisualization());	   
//		delete spViewer;
//		delete sp2DViewer;
//	}
	
	void NotestVisualizationMonodomainLR91DataRead(void)
	{
		TS_ASSERT_THROWS_NOTHING(
		                  spViewer=new MatlabVisualizer<1>(
		                  "pdes/tests/meshdata/NewMonodomainLR91_1d"));
		TS_ASSERT_THROWS_NOTHING(
					spViewer->CreateOutputFileForVisualization());	
		TS_ASSERT_THROWS_NOTHING(
					spViewer->CreateNodesFileForVisualization()
					);
		delete spViewer;

	}
	
};

#endif //_TESTMATLABVISUALIZER_HPP_
