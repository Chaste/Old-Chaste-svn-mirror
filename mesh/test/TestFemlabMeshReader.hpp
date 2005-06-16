// TestFemlabMeshReader.hpp

/**
 * Test suite for the FemlabMeshReader class.
 * 
 */

#ifndef _TESTFEMLABMESHREADER_HPP_
#define _TESTFEMLABMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "FemlabMeshReader.hpp"

static		AbstractMeshReader *spFemlabMeshReader;
class TestFemlabMeshReaders : public CxxTest::TestSuite
{
	public:
	
	/**
	 * Check that input files are opened correctly.
	 * 
	 */
	
	void TestFilesOpen(void)
	{
		TS_ASSERT_THROWS_NOTHING(
		                  spFemlabMeshReader=new FemlabMeshReader(
		                  "mesh/test/data/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat"));
	}
	
	/**
	 * Check that the nodes are read correctly. Checks that the output vector
	 * for a given input file is the correct length and that if the input file
	 * is corrupted (missing nodes) then an exception is thrown.
	 * 
	 */	 
	 
	void TestNodesDataRead(void)
	{
		spFemlabMeshReader=new FemlabMeshReader(
		                  "mesh/test/data/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		TS_ASSERT( spFemlabMeshReader->GetNumNodes() == 151); 
		
				
	}
	
	void TestDimension(void)
	{
		spFemlabMeshReader=new FemlabMeshReader(
		                  "mesh/test/data/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		TS_ASSERT( spFemlabMeshReader->GetDimension() == 2); 
	}
	
	/**
	 * Check that the elements are read correctly. Checks that the output vector
	 * for a given input file is the correct length and that if the input file
	 * is corrupted (missing elements) then an exception is thrown.
	 * 
	 */
	 
	void TestElementsDataRead(void)
	{
		spFemlabMeshReader=new FemlabMeshReader(
		                  "mesh/test/data/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		TS_ASSERT( spFemlabMeshReader->GetNumElements() == 260); 	
		
				
	}	

	/**
	 * Check that the faces are read correctly. Checks that the output vector
	 * for a given input file is the correct length and that if the input file
	 * is corrupted (missing faces) then an exception is thrown.
	 * 
	 */
	 
	void TestFacesDataRead(void)
	{
		spFemlabMeshReader=new FemlabMeshReader(
		                  "mesh/test/data/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		TS_ASSERT( spFemlabMeshReader->GetNumFaces() == 54); 	
		
				
	}		
	
	/**
	 * Check that GetNextNode() returns the coordinates of the correct node.
	 * Compares the coordinates of the first two nodes with their known
	 * values, checks that no errors are thrown for the remaining nodes and
	 * that an error is thrown if we try to call the function too many times.
	 * 
	 */
	 
	void TestGetNextNode(void)
	{
		spFemlabMeshReader=new FemlabMeshReader(
		                  "mesh/test/data/",
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
	
	/**
	 * Check that GetNextElement() works. Checks that no errors are thrown for 
	 * all of the elements and that an error is thrown if we try to call the 
	 * function too many times.
	 * 
	 */
	 
	void TestGetNextElement(void)
	{
		spFemlabMeshReader=new FemlabMeshReader(
		                  "mesh/test/data/",
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
	
	/**
	 * Check that GetNextBoundaryFace() works. Checks that no errors are thrown for 
	 * all of the elements and that an error is thrown if we try to call the 
	 * function too many times.
	 * 
	 */
	 
	void TestGetNextBoundaryFace(void)
	{
		spFemlabMeshReader=new FemlabMeshReader(
		                  "mesh/test/data/",
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
