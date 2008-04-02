/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CUBOIDMESHCONSTRUCTOR_HPP_
#define CUBOIDMESHCONSTRUCTOR_HPP_

#include "TrianglesMeshWriter.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"

//const double mesh_width = 0.2; // cm

template<unsigned DIM>
class CuboidMeshConstructor
{
private:

    void ConstructHyperCube(ConformingTetrahedralMesh<1,1> &rMesh, unsigned width)
    {
        rMesh.ConstructLinearMesh(width);
    }
    void ConstructHyperCube(ConformingTetrahedralMesh<2,2> &rMesh, unsigned width)
    {
        rMesh.ConstructRectangularMesh(width, width);
    }
    void ConstructHyperCube(ConformingTetrahedralMesh<3,3> &rMesh, unsigned width)
    {
        rMesh.ConstructCuboid(width, width, width);
    }
    
public:
    double mMeshWidth;
    unsigned NumElements;   
    unsigned NumNodes; 
    
    std::string Construct(unsigned meshNum, double meshWidth)
    {
        mMeshWidth=meshWidth;
        assert(meshNum < 30); //Sanity
        const std::string mesh_dir = "ConvergenceMesh";
        OutputFileHandler output_file_handler(mesh_dir);
        
        // create the mesh
        unsigned mesh_size = (unsigned) pow(2, meshNum+2); // number of elements in each dimension
        double scaling = mMeshWidth/(double) mesh_size;
        ConformingTetrahedralMesh<DIM,DIM> mesh;
        ConstructHyperCube(mesh, mesh_size);
        mesh.Scale(scaling, scaling, scaling);
        NumElements = mesh.GetNumElements();
        NumNodes = mesh.GetNumNodes();
        std::stringstream file_name_stream;
        file_name_stream<< "cube_" << DIM << "D_2mm_"<< NumElements <<"_elements";
        std::string mesh_filename = file_name_stream.str();
        
        if (output_file_handler.IsMaster())
        {
            TrianglesMeshWriter<DIM,DIM> mesh_writer(mesh_dir, mesh_filename, false);           
            mesh_writer.WriteFilesUsingMesh(mesh);
        }
        PetscTools::Barrier();
        
        std::string mesh_pathname = output_file_handler.GetOutputDirectoryFullPath()
                                  + mesh_filename;
                                    
        return mesh_pathname;
    }
    
    double GetWidth()
    {
        return mMeshWidth;
    }
    
};

#endif /*CUBOIDMESHCONSTRUCTOR_HPP_*/
