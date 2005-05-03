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

	void test1DVisualizationFileBuild(void)
	{
		TS_ASSERT_THROWS_NOTHING(
		                  spViewer=new MatlabVisualizer<1>(
		                  "pdes/tests/meshdata/trivial_1d_mesh"));
		spViewer->CreateNodesFileForVisualization();		                  
		
	}
	
	void test2DVisualizationFileBuild(void)
	{
		TS_ASSERT_THROWS_NOTHING(
		                  sp2DViewer=new MatlabVisualizer<2>(
		                  "pdes/tests/meshdata/disk_522_elements"));
		sp2DViewer->CreateNodesFileForVisualization();
	}
	
	/**
	 * Check whether it can read time file correctly.	 * 
	 */	
	void testVisualizationTimeFileRead(void)
	{
		TS_ASSERT_THROWS_NOTHING(
		                  spViewer=new MatlabVisualizer<1>(
		                  "pdes/tests/meshdata/trivial_1d_mesh"));
		spViewer->CreateOutputFileForVisualization();	
		
		TS_ASSERT_THROWS_NOTHING(
		                  sp2DViewer=new MatlabVisualizer<2>(
		                  "pdes/tests/meshdata/disk_522_elements"));
		sp2DViewer->CreateOutputFileForVisualization();			         
	}
	
};

#endif //_TESTMATLABVISUALIZER_HPP_
