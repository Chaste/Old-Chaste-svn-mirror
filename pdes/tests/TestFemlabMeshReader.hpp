// TestFemlabMeshReader.hpp

/**
 * Test suite for the FEMLABMESHReader class.
 * 
 */

#ifndef _TESTFEMLABMESHREADER_HPP_
#define _TESTFEMLABMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "../FemlabMeshReader.hpp"

static		AbstractMeshReader *spFemlabMeshReader;
class TestFemlabMeshReaders : public CxxTest::TestSuite
{
	public:
	
	/**
	 * Check that input files are opened correctly.
	 * 
	 */
	
	void testFilesOpen(void)
	{
		TS_ASSERT_THROWS_NOTHING(
		                  spFemlabMeshReader=new FemlabMeshReader(
		                  "pdes/tests/meshdata/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat"));
		
	}
	
	void testNodesDataRead(void)
	{
		spFemlabMeshReader=new FemlabMeshReader(
		                  "pdes/tests/meshdata/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		TS_ASSERT( spFemlabMeshReader->GetNumNodes() == 151); 
		
				
	}
	
	void testDimension(void)
	{
		spFemlabMeshReader=new FemlabMeshReader(
		                  "pdes/tests/meshdata/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		TS_ASSERT( spFemlabMeshReader->GetDimension() == 2); 
	}
	
	void testElementsDataRead(void)
	{
		spFemlabMeshReader=new FemlabMeshReader(
		                  "pdes/tests/meshdata/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		TS_ASSERT( spFemlabMeshReader->GetNumElements() == 260); 	
		
				
	}	

	void testFacesDataRead(void)
	{
		spFemlabMeshReader=new FemlabMeshReader(
		                  "pdes/tests/meshdata/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		TS_ASSERT( spFemlabMeshReader->GetNumFaces() == 54); 	
		
				
	}		
	
	void testGetNextNode(void)
	{
		spFemlabMeshReader=new FemlabMeshReader(
		                  "pdes/tests/meshdata/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		std::vector<double> FirstNode;
		                  
		FirstNode = spFemlabMeshReader->GetNextNode();
		
		TS_ASSERT_DELTA( FirstNode[0] ,  0.0 , 1e-6 );
		TS_ASSERT_DELTA( FirstNode[1] , 1.0 , 1e-6 )
		
		std::vector<double> NextNode;
		                  
		NextNode = spFemlabMeshReader->GetNextNode();
		
		TS_ASSERT_DELTA( NextNode[0] , 0.5 , 1e-6 );
		TS_ASSERT_DELTA( NextNode[1] , 1.0 , 1e-6 )
			    		
		for (int i = 2; i < spFemlabMeshReader->GetNumNodes(); i++)
		{
			TS_ASSERT_THROWS_NOTHING(NextNode = spFemlabMeshReader->GetNextNode());
		}
		
		TS_ASSERT_THROWS_ANYTHING(NextNode = spFemlabMeshReader->GetNextNode());
		
	}
	
	void testGetNextElement(void)
	{
		spFemlabMeshReader=new FemlabMeshReader(
		                  "pdes/tests/meshdata/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		std::vector<int> FirstElement;
		                  
		FirstElement = spFemlabMeshReader->GetNextElement();
		
		TS_ASSERT( FirstElement[0]==15);
		TS_ASSERT( FirstElement[1]==3);
		TS_ASSERT( FirstElement[2]==62);		

		std::vector<int> NextElement;
		                  
		NextElement = spFemlabMeshReader->GetNextElement();
		
		TS_ASSERT( NextElement[0]==8);
		TS_ASSERT( NextElement[1]==0);
		TS_ASSERT( NextElement[2]==53);		
			    		
		for (int i = 2; i < spFemlabMeshReader->GetNumElements(); i++)
		{
			TS_ASSERT_THROWS_NOTHING(NextElement = spFemlabMeshReader->GetNextElement());
		}
		
		TS_ASSERT_THROWS_ANYTHING(NextElement = spFemlabMeshReader->GetNextElement());
		
	}
	
	void testGetNextBoundaryFace(void)
	{
		spFemlabMeshReader=new FemlabMeshReader(
		                  "pdes/tests/meshdata/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		std::vector<int> FirstBoundaryFace;
		                  
		FirstBoundaryFace = spFemlabMeshReader->GetNextBoundaryFace();
		
		TS_ASSERT( FirstBoundaryFace[0]==0);
		TS_ASSERT( FirstBoundaryFace[1]==8);

		std::vector<int> NextBoundaryFace;
		                  
		NextBoundaryFace = spFemlabMeshReader->GetNextBoundaryFace();
		
		TS_ASSERT( NextBoundaryFace[0]==8);
		TS_ASSERT( NextBoundaryFace[1]==9);
			    		
		for (int i = 2; i < spFemlabMeshReader->GetNumBoundaryFaces(); i++)
		{
			TS_ASSERT_THROWS_NOTHING(NextBoundaryFace = spFemlabMeshReader->GetNextBoundaryFace());
		}
		
		TS_ASSERT_THROWS_ANYTHING(NextBoundaryFace = spFemlabMeshReader->GetNextBoundaryFace());
		
	}
		
};

#endif //_TESTFEMLABMESHREADER_HPP_
