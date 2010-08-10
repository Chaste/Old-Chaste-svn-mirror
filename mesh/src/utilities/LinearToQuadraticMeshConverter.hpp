/*

Copyright (C) University of Oxford, 2005-2010

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


#ifndef LINEARTOQUADRATICMESHCONVERTER_HPP_
#define LINEARTOQUADRATICMESHCONVERTER_HPP_

#include "QuadraticMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "CommandLineArguments.hpp"

// Note: currently completely verbose..

/**
 *  Takes in a linear triangle or tetgen format mesh and converts it to a quadratic mesh.
 *  2D or 3D only.
 * 
 *  Does this in THREE-DIMENSIONS by: 
 *   (i)   Calling tetgen with '-o2nn' to convert node and ele files, and rewrite
 *         face file with containing element info; 
 *   (ii)  Calling python/RemoveAttributeColumn.py to remove the zero-attribute 
 *         column from the new face file (which would cause Chaste to think all 
 *         faces are internal); 
 *   (iii) Using QuadraticMesh to read in the new files, during which it converts 
 *         the linear faces to quadratic (bit slow but only has to be done this one time),
 *         and writing a new face file
 *   (iv)  Move the new files to the output directory.
 *  
 *  In TWO-DIMENSIONS, does the same except uses triangle with '-o2e' (e for write edges), 
 *  as triangle doesn't seem to have an '-nn' option. Therefore (ii) is not needed, and
 *  a slower conversion of linear faces to quad faces takes place in (iii). 
 */
template<unsigned DIM>
class LinearToQuadraticMeshConverter
{
public:

    /** 
     *  Constructor does all the work
     * 
     *  @param inMeshDirectory Directory of input mesh, e.g. mesh/test/data
     *  @param inMeshStem Stem for linear mesh, e.g. "cube_136_elements"
     *  @param outputDirectory Where to put the new files. NOT CLEANED. They will be
     *     named after the inMeshStem with "_quadratic" appended, eg 
     *     "cube_136_elements_quadratic"
     */
    LinearToQuadraticMeshConverter(std::string inMeshDirectory, std::string inMeshStem, std::string outputDirectory)
    {
        assert(DIM==2 || DIM==3);
        
        // run tetgen with -o2 (create quadratic elements) and -nn (add containing element info to face file) arguments
        // on the full initial mesh (.node and .ele)
        std::string full_in_mesh_stem = inMeshDirectory + "/" + inMeshStem;
        
        std::string binary = (DIM==2) ? "triangle -o2e " : "tetgen -o2nn "; 
        std::string faceformat = (DIM==2) ? ".edge" : ".face";
        unsigned column = (DIM==2) ? 4 : 5;
        
        std::string command = binary + full_in_mesh_stem + ".node " + full_in_mesh_stem + ".ele";
        std::cout << "\n[running command]: " << command << "\n" << std::flush;
        system(command.c_str());

        if(DIM==3)
        {
            std::stringstream ss;
            ss << "python projects/pras/test/RemoveAttributeColumn.py " << column 
               << " < " << full_in_mesh_stem << ".1" << faceformat << " > " 
            << "temp.face";
            std::cout << "\n[running command]: " << ss.str() << "\n" << std::flush;
            system(ss.str().c_str());

            command = std::string("mv temp.face ") + full_in_mesh_stem + ".1" + faceformat;              
            std::cout << "\n[running command]: " << command << "\n" << std::flush;
            system(command.c_str());
        }

        // read in the new mesh to a QuadraticMesh. The false says the face file is linear (so it should 
        // convert to quad elements. The true says the face file has containing element info (makes the conversion
        // fast).
        std::cout << "\n[Chaste] Reading files using QuadraticMesh, converting face file, and writing.\n" << std::flush;

        bool face_file_has_containing_elements = (DIM==3);
        std::string new_mesh = full_in_mesh_stem + ".1";
        TrianglesMeshReader<DIM,DIM> reader(new_mesh, 2 /* quad elements */, 1 /* linear faces */, face_file_has_containing_elements);
        
        QuadraticMesh<DIM> mesh;
        mesh.ConstructFromMeshReader(reader);
        
        // Write out the quadratic face file to the output dir
        std::string out_face_file = inMeshStem + "_quadratic" + faceformat;
        mesh.WriteBoundaryElementFile(outputDirectory,out_face_file);
                
        // Move the quadratic (tetgen-created) node and ele files to the output dir 
        OutputFileHandler handler(outputDirectory,false);
        std::string full_output_dir = handler.GetOutputDirectoryFullPath();
        
        command =   std::string("mv ") + full_in_mesh_stem + ".1.node "
                  + full_output_dir + "/" + inMeshStem + "_quadratic.node"; 
        std::cout << "\n[running command]: " << command << "\n" << std::flush;
        system(command.c_str());
        
        command =   std::string("mv ") + full_in_mesh_stem + ".1.ele "
                  + full_output_dir + "/" + inMeshStem + "_quadratic.ele"; 
        std::cout << "\n[running command]: " << command << "\n" << std::flush;
        system(command.c_str());

        // Delete the .1.face and .1.neigh files (latter also created by using -nn in tetgen)
        command = "rm -f " + full_in_mesh_stem + ".1" + faceformat + " " + full_in_mesh_stem + ".1.neigh";
        std::cout << "\n[running command]: " << command << "\n" << std::flush;
        system(command.c_str());
    }
};
    


#endif /*LINEARTOQUADRATICMESHCONVERTER_HPP_*/
