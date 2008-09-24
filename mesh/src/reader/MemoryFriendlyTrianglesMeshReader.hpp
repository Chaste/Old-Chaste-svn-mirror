/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef _MEMORYFRIENDLYTRIANGLESMESHREADER_HPP_
#define _MEMORYFRIENDLYTRIANGLESMESHREADER_HPP_

#include <vector>
#include <string>
#include <fstream>
#include <cassert>
#include "AbstractMeshReader.hpp"

const static char* NODES_FILE_EXTENSION = ".node";
const static char* ELEMENTS_FILE_EXTENSION = ".ele";
const static char* FACES_FILE_EXTENSION = ".face";
const static char* EDGES_FILE_EXTENSION = ".edge";

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MemoryFriendlyTrianglesMeshReader : public AbstractMeshReader<ELEMENT_DIM,SPACE_DIM>
{
private:

    bool mIndexFromZero; /**< True if input data is numbered from zero, false otherwise */

	std::string mFilesBaseName;

	std::ifstream mNodesFile;
	std::ifstream mElementsFile;
	std::ifstream mFacesFile;

	unsigned mNumNodes;
	unsigned mNumElements;
	unsigned mNumFaces;
	
	unsigned mNodesRead;
	unsigned mElementsRead;
	unsigned mFacesRead;

    unsigned mNumNodeAttributes; /**< Is the number of attributes stored at each node */
    unsigned mMaxNodeBdyMarker; /**< Is the maximum node boundary marker */
    unsigned mNumElementNodes; /** Is the number of nodes per element*/
    unsigned mNumElementAttributes; /**< Is the number of attributes stored for each element */
    unsigned mMaxFaceBdyMarker; /**< Is the maximum face (or edge) boundary marker */

	unsigned mOrderOfElements;
    unsigned mNodesPerElement;

public:
    MemoryFriendlyTrianglesMeshReader(std::string pathBaseName, unsigned orderOfElements=1):
        mFilesBaseName(pathBaseName),
    	mNumNodes(0),
    	mNumElements(0),
    	mNumFaces(0),
		mNodesRead(0),
		mElementsRead(0),
		mFacesRead(0),
    	mOrderOfElements(orderOfElements)
    {
        // Only linear and quadratic elements
        assert(orderOfElements==1 || orderOfElements==2);
        if(mOrderOfElements==1)
        {
            mNodesPerElement = ELEMENT_DIM+1;
        }
        else
        {
            assert(SPACE_DIM==ELEMENT_DIM);
            mNodesPerElement = (ELEMENT_DIM+1)*(ELEMENT_DIM+2)/2;
        }

        mIndexFromZero = false; // Initially assume that nodes are not numbered from zero
        

    	OpenFiles();

		ReadHeaders();

    }
    
 	/**< Returns the number of elements in the mesh */
    unsigned GetNumElements() const
    {
        return mNumElements;
    }
    
    /**< Returns the number of nodes in the mesh */
    unsigned GetNumNodes() const
    {
        return mNumNodes;
    } 
    
    /**< Returns the number of faces in the mesh (synonym of GetNumEdges()) */
    unsigned GetNumFaces() const
    {
        return mNumFaces;
    } 
    
    /**< Returns the number of edges in the mesh (synonym of GetNumFaces()) */
    unsigned GetNumEdges() const
    {
        return mNumFaces;
    }


	/**< Returns a vector of the coordinates of each node in turn */
    std::vector<double> GetNextNode()
    {
    	std::vector<double> ret_coords;
        
        std::string buffer;		
		GetNextLineFromStream(mNodesFile, buffer);

		std::stringstream buffer_stream(buffer);

		unsigned index;		
		buffer_stream >> index;
		
        unsigned offset = mIndexFromZero ? 0 : 1;
		if(index != mNodesRead+offset)
		{
			std::stringstream error;
			error << "Data for node " << mNodesRead << " missing";
			EXCEPTION(error.str());
		}

        double coord;
    	
	    for (unsigned i = 0; i < SPACE_DIM; i++)
        {
	        buffer_stream >> coord;
	        ret_coords.push_back(coord);
        }
                
        mNodesRead++;
        
        return ret_coords;
    }
    
    /**< Resets pointers to beginning*/
    void Reset()
    {
        CloseFiles();
        OpenFiles();
        ReadHeaders();
    } 
  
  	/**< Returns a vector of the nodes of each element in turn */
    std::vector<unsigned> GetNextElement()
    {
		std::vector<unsigned> ret_indices;		
        
        std::string buffer;		
		GetNextLineFromStream(mElementsFile, buffer);

		std::stringstream buffer_stream(buffer);

		unsigned element_index;		
		buffer_stream >> element_index;
		
        unsigned offset = mIndexFromZero ? 0 : 1;
		if(element_index != mElementsRead+offset)
		{
			std::stringstream error;
			error << "Data for element " << mElementsRead << " missing";
			EXCEPTION(error.str());
		}

        unsigned node_index;
	    for (unsigned i = 0; i < mNodesPerElement; i++)
        {
	        buffer_stream >> node_index;
	        ret_indices.push_back(node_index - offset);
        }
                
        mElementsRead++;

		return ret_indices;
    } 
        
    /**< Returns a vector of the nodes of each face in turn (synonym of GetNextEdge()) */
    std::vector<unsigned> GetNextFace()
    {
		std::vector<unsigned> ret_indices;		

		if (SPACE_DIM == 1)
		{
			ret_indices.push_back(mFacesRead);
		}
		else
		{
            unsigned offset = mIndexFromZero ? 0 : 1;
            
            unsigned is_boundary;
            unsigned element_dim=ELEMENT_DIM;//In case ELEMENT_DIM is erroneously instatiated to zero
            assert(element_dim != 0); //Covered in earlier exception, but needed in loop guard here.
            do
            {       
                ret_indices.clear();
                     
		        std::string buffer;		
				GetNextLineFromStream(mFacesFile, buffer);
		
				std::stringstream buffer_stream(buffer);
		
				unsigned face_index;		
				buffer_stream >> face_index;
		
		        unsigned node_index;
                
			    for (unsigned i = 0; i<element_dim; i++)
		        {
			        buffer_stream >> node_index;
			        ret_indices.push_back(node_index-offset);
		        }

                buffer_stream >> is_boundary; 
            }
            while((mMaxFaceBdyMarker==1) && (is_boundary!=1));
		}
	                
        mFacesRead++;
		return ret_indices;
    }

    /**< Returns a vector of the nodes of each edge in turn (synonym of GetNextFace()) */
    std::vector<unsigned> GetNextEdge()
    {
		return GetNextFace();
    } 


private:

	void OpenFiles()
	{
        OpenNodeFile();
        OpenElementsFile();
        OpenFacesFile();
    }
    
    void OpenNodeFile()
    {        
		// Nodes definition
    	std::string file_name = mFilesBaseName + NODES_FILE_EXTENSION;
		mNodesFile.open(file_name.c_str());
		if (!mNodesFile.is_open())
    	{
        	EXCEPTION("Could not open data file: "+file_name);
    	}
    }
    
    void OpenElementsFile()
    {
		// Elements definition
		std::string file_name;
        if (ELEMENT_DIM == SPACE_DIM)
		{
    		file_name = mFilesBaseName + ELEMENTS_FILE_EXTENSION;
		}
		else
		{
		    if (SPACE_DIM == 2 && ELEMENT_DIM == 1)
		    {
		        file_name = mFilesBaseName + EDGES_FILE_EXTENSION;
		    }
		    else if (SPACE_DIM == 3 && ELEMENT_DIM == 2)
		    {
		        file_name = mFilesBaseName + FACES_FILE_EXTENSION;
		    }
		    else
		    {
		        EXCEPTION("Can't have a zero-dimensional mesh in a one-dimensional space or a one-dimensional mesh in a three-dimensional space");
		    }
		}
    	
		mElementsFile.open(file_name.c_str());
		if (!mElementsFile.is_open())
    	{
        	EXCEPTION("Could not open data file: "+file_name);
    	}
    }

    void OpenFacesFile()
    {		
		// Faces/edges definition
        std::string file_name;
	    if (SPACE_DIM == 3)
    	{
    		if (SPACE_DIM == ELEMENT_DIM)
    		{
	    		file_name = mFilesBaseName + FACES_FILE_EXTENSION;
    		}
    		else
    		{
    			file_name = mFilesBaseName + EDGES_FILE_EXTENSION;
    		}
	    }
	    else if (SPACE_DIM == 2)
	    { 
	    	file_name = mFilesBaseName + EDGES_FILE_EXTENSION;
	    }
	    else if (SPACE_DIM == 1)
	    {
	    	//There is no file, data will be generated instead of read
	    	return;
	    }
		
		mFacesFile.open(file_name.c_str());
		if (!mFacesFile.is_open())
    	{
        	EXCEPTION("Could not open data file: "+file_name);
    	}		    			
	}
	
	void ReadHeaders()
	{
		std::string buffer;
		
		GetNextLineFromStream(mNodesFile, buffer);
		std::stringstream buffer_stream(buffer);
		unsigned dimension;
	    buffer_stream >> mNumNodes >> dimension >> mNumNodeAttributes >> mMaxNodeBdyMarker;
	    if (SPACE_DIM != dimension)
	    {
	        EXCEPTION("SPACE_DIM  != dimension read from file ");
	    }

        // get the next line to see if it is indexed from zero or not        	    
        GetNextLineFromStream(mNodesFile, buffer);
        std::stringstream buffer_stream_ii(buffer);
        unsigned first_index;
        buffer_stream_ii >> first_index;
        if(first_index==0)
        {
            mIndexFromZero = true;
        }
        else
        {
            if(first_index!=1)
            {
                EXCEPTION("Nodes should indexed from zero or one");
            }
            mIndexFromZero = false;
        }
        //close, reopen, skip header
        mNodesFile.close();
        OpenNodeFile();
        GetNextLineFromStream(mNodesFile, buffer);     
        
        
		/// \todo: rename std::stringstream variables
		GetNextLineFromStream(mElementsFile, buffer);
		std::stringstream buffer_stream2(buffer);
		
		if (ELEMENT_DIM == SPACE_DIM)
		{ 
		    buffer_stream2 >> mNumElements >> mNumElementNodes >> mNumElementAttributes;
	    
	    	if( mNumElementNodes != mNodesPerElement )
	    	{
	        	std::stringstream error;
	        	error << "Number of nodes per elem, " << mNumElementNodes << ", does not match "
		              << "expected number, " << mNodesPerElement << " (which is calculated given "
		              << "the order of elements chosen, " << mOrderOfElements << " (1=linear, 2=quadratics)";
		        EXCEPTION(error.str());
		    }
		}
		else
		{
	    	buffer_stream2 >> mNumElements >> mMaxFaceBdyMarker;

        	mNodesPerElement = ELEMENT_DIM+1;			
		}
	    		

		if (SPACE_DIM == 1)
		{
			mNumFaces = mNumNodes;
		}
		else
		{
			GetNextLineFromStream(mFacesFile, buffer);		
			std::stringstream buffer_stream3(buffer);

		    buffer_stream3 >> mNumFaces >> mMaxFaceBdyMarker;
		    assert(mMaxFaceBdyMarker==0 || mMaxFaceBdyMarker==1);
		    
		    // if mMaxFaceBdyMarker=1 then loop over and set mNumFaces to be
		    // the number of faces which are marked as boundary faces
			if((mMaxFaceBdyMarker==1) && (SPACE_DIM!=1))
			{
                unsigned num_boundary_faces = 0;
                bool end_of_file=false;
                while(!end_of_file)
                {
                    try
                    {
         				GetNextFace();
                        num_boundary_faces++;
                    }
                    catch(Exception& e)
                    {
                        end_of_file = true;
                    }
				}
				mNumFaces = num_boundary_faces;
				
				// close the file, reopen, and skip the header again
				mFacesFile.close();
                OpenFacesFile();
                GetNextLineFromStream(mFacesFile, buffer);      
            }
		}		
	}
	
    void CloseFiles()
    {
        mNodesFile.close();
        mElementsFile.close();
        mFacesFile.close();                
    }
    
	void GetNextLineFromStream(std::ifstream& fileStream, std::string& rawLine)
	{
		bool line_is_blank;
		
		do
		{
	    	getline(fileStream, rawLine);

	    	if (fileStream.eof())
	    	{
	    		/// \todo: improve this error message
	    		EXCEPTION("File contains incomplete data");
	    	}
	
            // Get rid of any comment
            rawLine = rawLine.substr(0,rawLine.find('#',0));
            
            line_is_blank = (rawLine.find_first_not_of(" \t",0) == std::string::npos);	
		}
		while (line_is_blank);		
	}    
};

#endif //_MEMORYFRIENDLYTRIANGLESMESHREADER_HPP_
