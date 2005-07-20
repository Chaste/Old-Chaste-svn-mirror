// TestFemlabMeshReader.hpp

/**
 * Test suite for the FemlabMeshReader class.
 * 
 */

#ifndef _TESTFEMLABMESHREADER_HPP_
#define _TESTFEMLABMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "FemlabMeshReader.hpp"


class TestFemlabMeshReaders : public CxxTest::TestSuite
{
	public:
     
    AbstractMeshReader *pFemlabMeshReader;	
	/**
	 * Check that input files are opened correctly.
	 * 
	 */
	
	void TestFilesOpen(void)
	{
		TS_ASSERT_THROWS_NOTHING(pFemlabMeshReader=new FemlabMeshReader(
		                  "mesh/test/data/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat"));
		                  
		delete pFemlabMeshReader;
	}
	
	/**
	 * Check that the nodes are read correctly. Checks that the output vector
	 * for a given input file is the correct length and that if the input file
	 * is corrupted (missing nodes) then an exception is thrown.
	 * 
	 */	 
	 
	void TestNodesDataRead(void)
	{
		pFemlabMeshReader=new FemlabMeshReader(
		                  "mesh/test/data/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		TS_ASSERT( pFemlabMeshReader->GetNumNodes() == 151); 
		
		delete pFemlabMeshReader;
	}
	
	void TestDimension(void)
	{
		pFemlabMeshReader=new FemlabMeshReader(
		                  "mesh/test/data/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		TS_ASSERT( pFemlabMeshReader->GetDimension() == 2); 
		
		delete pFemlabMeshReader;
	}
	
	/**
	 * Check that the elements are read correctly. Checks that the output vector
	 * for a given input file is the correct length and that if the input file
	 * is corrupted (missing elements) then an exception is thrown.
	 * 
	 */
	 
	void TestElementsDataRead(void)
	{
		pFemlabMeshReader=new FemlabMeshReader(
		                  "mesh/test/data/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		TS_ASSERT( pFemlabMeshReader->GetNumElements() == 260); 	
		
		delete pFemlabMeshReader;
				
	}	

	/**
	 * Check that the faces are read correctly. Checks that the output vector
	 * for a given input file is the correct length and that if the input file
	 * is corrupted (missing faces) then an exception is thrown.
	 * 
	 */
	 
	void TestFacesDataRead(void)
	{
		pFemlabMeshReader=new FemlabMeshReader(
		                  "mesh/test/data/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		TS_ASSERT( pFemlabMeshReader->GetNumFaces() == 54); 	
		
		delete pFemlabMeshReader;
				
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
		pFemlabMeshReader=new FemlabMeshReader(
		                  "mesh/test/data/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		std::vector<double> FirstNode;
		                  
		FirstNode = pFemlabMeshReader->GetNextNode();
		
		TS_ASSERT_DELTA( FirstNode[0] ,  0.0 , 1e-6 );
		TS_ASSERT_DELTA( FirstNode[1] , 1.0 , 1e-6 )
		
		std::vector<double> NextNode;
		                  
		NextNode = pFemlabMeshReader->GetNextNode();
		
		TS_ASSERT_DELTA( NextNode[0] , 0.5 , 1e-6 );
		TS_ASSERT_DELTA( NextNode[1] , 1.0 , 1e-6 )
			    		
		for (int i = 2; i < pFemlabMeshReader->GetNumNodes(); i++)
		{
			TS_ASSERT_THROWS_NOTHING(NextNode = pFemlabMeshReader->GetNextNode());
		}
		
		TS_ASSERT_THROWS_ANYTHING(NextNode = pFemlabMeshReader->GetNextNode());
		
		delete pFemlabMeshReader;
		
	}
	
	/**
	 * Check that GetNextElement() works. Checks that no errors are thrown for 
	 * all of the elements and that an error is thrown if we try to call the 
	 * function too many times.
	 * 
	 */
	 
	void TestGetNextElement(void)
	{
		pFemlabMeshReader=new FemlabMeshReader(
		                  "mesh/test/data/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		std::vector<int> FirstElement;
		                  
		FirstElement = pFemlabMeshReader->GetNextElement();
		
		TS_ASSERT( FirstElement[0]==15);
		TS_ASSERT( FirstElement[1]==3);
		TS_ASSERT( FirstElement[2]==62);		

		std::vector<int> NextElement;
		                  
		NextElement = pFemlabMeshReader->GetNextElement();
		
		TS_ASSERT( NextElement[0]==8);
		TS_ASSERT( NextElement[1]==0);
		TS_ASSERT( NextElement[2]==53);		
			    		
		for (int i = 2; i < pFemlabMeshReader->GetNumElements(); i++)
		{
			TS_ASSERT_THROWS_NOTHING(NextElement = pFemlabMeshReader->GetNextElement());
		}
		
		TS_ASSERT_THROWS_ANYTHING(NextElement = pFemlabMeshReader->GetNextElement());
		
		delete pFemlabMeshReader;
		
	}
	
	/**
	 * Check that GetNextBoundaryFace() works. Checks that no errors are thrown for 
	 * all of the elements and that an error is thrown if we try to call the 
	 * function too many times.
	 * 
	 */
	 
	void TestGetNextBoundaryFace(void)
	{
		pFemlabMeshReader=new FemlabMeshReader(
		                  "mesh/test/data/",
		                  "femlab_lshape_nodes.dat",
		                  "femlab_lshape_elements.dat",
		                  "femlab_lshape_edges.dat");
		
		std::vector<int> FirstBoundaryFace;
		                  
		FirstBoundaryFace = pFemlabMeshReader->GetNextBoundaryFace();
		
		TS_ASSERT( FirstBoundaryFace[0]==0);
		TS_ASSERT( FirstBoundaryFace[1]==8);

		std::vector<int> NextBoundaryFace;
		                  
		NextBoundaryFace = pFemlabMeshReader->GetNextBoundaryFace();
		
		TS_ASSERT( NextBoundaryFace[0]==8);
		TS_ASSERT( NextBoundaryFace[1]==9);
			    		
		for (int i = 2; i < pFemlabMeshReader->GetNumBoundaryFaces(); i++)
		{
			TS_ASSERT_THROWS_NOTHING(NextBoundaryFace = pFemlabMeshReader->GetNextBoundaryFace());
		}
		
		TS_ASSERT_THROWS_ANYTHING(NextBoundaryFace = pFemlabMeshReader->GetNextBoundaryFace());
		
		delete pFemlabMeshReader;
		
	}
		
};

#endif //_TESTFEMLABMESHREADER_HPP_
