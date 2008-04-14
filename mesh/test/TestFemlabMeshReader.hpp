/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/

// TestFemlabMeshReader.hpp

/**
 * Test suite for the FemlabMeshReader class.
 *
 */

#ifndef _TESTFEMLABMESHREADER_HPP_
#define _TESTFEMLABMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "FemlabMeshReader.hpp"

typedef FemlabMeshReader<2,2> READER_2D;
typedef FemlabMeshReader<1,1> READER_1D;

class TestFemlabMeshReaders : public CxxTest::TestSuite
{
public:

    AbstractMeshReader<2,2> *pFemlabMeshReader;
    /**
     * Check that input files are opened correctly.
     * 
     */
    
    void TestFilesOpen(void)
    {
        TS_ASSERT_THROWS_NOTHING(pFemlabMeshReader=new READER_2D(
                                                       "mesh/test/data/",
                                                       "femlab_lshape_nodes.dat",
                                                       "femlab_lshape_elements.dat",
                                                       "femlab_lshape_edges.dat"));
                                                       
        delete pFemlabMeshReader;
        
        // Coverage test
        
        AbstractMeshReader<1,1> *pWrongFemlabMeshReader;
        
        TS_ASSERT_THROWS_ANYTHING(pWrongFemlabMeshReader=new READER_1D(
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
        pFemlabMeshReader=new FemlabMeshReader<2,2>(
                              "mesh/test/data/",
                              "femlab_lshape_nodes.dat",
                              "femlab_lshape_elements.dat",
                              "femlab_lshape_edges.dat");
                              
        TS_ASSERT( pFemlabMeshReader->GetNumNodes() == 151);
        
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
        pFemlabMeshReader=new FemlabMeshReader<2,2>(
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
        pFemlabMeshReader=new FemlabMeshReader<2,2>(
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
        pFemlabMeshReader=new FemlabMeshReader<2,2>(
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
        
        for (unsigned i = 2; i < pFemlabMeshReader->GetNumNodes(); i++)
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
        pFemlabMeshReader=new FemlabMeshReader<2,2>(
                              "mesh/test/data/",
                              "femlab_lshape_nodes.dat",
                              "femlab_lshape_elements.dat",
                              "femlab_lshape_edges.dat");
                              
        std::vector<unsigned> FirstElement;
        
        FirstElement = pFemlabMeshReader->GetNextElement();
        
        TS_ASSERT( FirstElement[0]==15);
        TS_ASSERT( FirstElement[1]==3);
        TS_ASSERT( FirstElement[2]==62);
        
        std::vector<unsigned> NextElement;
        
        NextElement = pFemlabMeshReader->GetNextElement();
        
        TS_ASSERT( NextElement[0]==8);
        TS_ASSERT( NextElement[1]==0);
        TS_ASSERT( NextElement[2]==53);
        
        for (unsigned i = 2; i < pFemlabMeshReader->GetNumElements(); i++)
        {
            TS_ASSERT_THROWS_NOTHING(NextElement = pFemlabMeshReader->GetNextElement());
        }
        
        TS_ASSERT_THROWS_ANYTHING(NextElement = pFemlabMeshReader->GetNextElement());
        
        delete pFemlabMeshReader;
        
    }
    
    /**
     * Check that GetNextFace() works. Checks that no errors are thrown for 
     * all of the elements and that an error is thrown if we try to call the 
     * function too many times.
     * 
     */
    
    void TestGetNextFace(void)
    {
        pFemlabMeshReader=new FemlabMeshReader<2,2>(
                              "mesh/test/data/",
                              "femlab_lshape_nodes.dat",
                              "femlab_lshape_elements.dat",
                              "femlab_lshape_edges.dat");
                              
        std::vector<unsigned> FirstFace;
        
        FirstFace = pFemlabMeshReader->GetNextFace();
        
        TS_ASSERT( FirstFace[0]==0);
        TS_ASSERT( FirstFace[1]==8);
        
        std::vector<unsigned> NextFace;
        
        NextFace = pFemlabMeshReader->GetNextFace();
        
        TS_ASSERT( NextFace[0]==8);
        TS_ASSERT( NextFace[1]==9);
        
        for (unsigned i = 2; i < pFemlabMeshReader->GetNumFaces(); i++)
        {
            TS_ASSERT_THROWS_NOTHING(NextFace = pFemlabMeshReader->GetNextFace());
        }
        
        TS_ASSERT_THROWS_ANYTHING(NextFace = pFemlabMeshReader->GetNextFace());
        
        delete pFemlabMeshReader;
        
    }
    
};

#endif //_TESTFEMLABMESHREADER_HPP_
