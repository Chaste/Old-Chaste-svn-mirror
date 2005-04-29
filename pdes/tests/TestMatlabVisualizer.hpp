// TestMatlabVisualizer.hpp

/**
 * Test suite for the MatlabVisualizer class
 * 
 */

#ifndef _TESTMATLABVISUALIZER_HPP_
#define _TESTMATLABVISUALIZER_HPP_

#include <cxxtest/TestSuite.h>
#include "../AbstractVisualizer.hpp"
#include "../MatlabVisualizer.hpp"

static AbstractVisualizer *spViewer;
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
		                  spViewer=new MatlabVisualizer(
		                  "pdes/tests/meshdata/trivial_1d_mesh",1));
		spViewer->CreateNodesFileForVisualization();		                  
		
	}
	
	void test2DVisualizationFileBuild(void)
	{
		TS_ASSERT_THROWS_NOTHING(
		                  spViewer=new MatlabVisualizer(
		                  "pdes/tests/meshdata/disk_522_elements",2));
		spViewer->CreateNodesFileForVisualization();
	}
	
	/**
	 * Check whether it can read time file correctly.	 * 
	 */	
	void testVisualizationTimeFileRead(void)
	{
		TS_ASSERT_THROWS_NOTHING(
		                  spViewer=new MatlabVisualizer(
		                  "pdes/tests/meshdata/trivial_1d_mesh",1));
		spViewer->CreateOutputFileForVisualization();	
		
		TS_ASSERT_THROWS_NOTHING(
		                  spViewer=new MatlabVisualizer(
		                  "pdes/tests/meshdata/disk_522_elements",2));
		spViewer->CreateOutputFileForVisualization();			         
	}
	
};

#endif //_TESTMATLABVISUALIZER_HPP_
