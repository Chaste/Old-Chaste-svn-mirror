/*

Copyright (C) University of Oxford, 2005-2011

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


#ifndef TESTFIBREWRITER_HPP_
#define TESTFIBREWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include <vector>

#include "FibreReader.hpp"
#include "FibreWriter.hpp"
#include "PetscTools.hpp"
#include "UblasIncludes.hpp"


class TestFibreWriter : public CxxTest::TestSuite
{
public:
    void TestAxiWriterAscii()
    {
        FileFinder file_finder("heart/test/data/fibre_tests/SimpleAxisymmetric2.axi", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file_finder, AXISYM);
        std::vector< c_vector<double, 3> > fibre_vector;
        fibre_reader.GetAllAxi(fibre_vector);
        
        //Write ascii file
        FibreWriter<3> fibre_writer("TestFibreWriter", "SimpleAxisymmetric2", true);
        fibre_writer.WriteAllAxi(fibre_vector);
        
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestFibreWriter/";
        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + results_dir + "/SimpleAxisymmetric2.axi heart/test/data/fibre_tests/SimpleAxisymmetric2.axi").c_str()), 0);
        
    }
    void TestAxiWriterBinary()
    {
        FileFinder file_finder("heart/test/data/fibre_tests/SimpleAxisymmetric2.axi", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file_finder, AXISYM);
        std::vector< c_vector<double, 3> > fibre_vector;
        fibre_reader.GetAllAxi(fibre_vector);
        
        //Write binary file
        FibreWriter<3> fibre_writer("TestFibreWriter", "SimpleAxisymmetric2Bin", false);
        fibre_writer.SetWriteFileAsBinary();
        fibre_writer.WriteAllAxi(fibre_vector);
        
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestFibreWriter/";
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/SimpleAxisymmetric2Bin.axi heart/test/data/fibre_tests/SimpleAxisymmetric2Bin.axi").c_str()), 0);
    }
    
    void TestOrthoWriterAscii() throw (Exception)
    {
        FileFinder file_finder("heart/test/data/fibre_tests/Orthotropic3D.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file_finder, ORTHO);
        std::vector< c_vector<double, 3> > fibres;
        std::vector< c_vector<double, 3> > second;
        std::vector< c_vector<double, 3> > third;
        fibre_reader.GetAllOrtho(fibres, second, third);
          
        //Write ascii file
        FibreWriter<3> fibre_writer("TestFibreWriter", "Orthotropic3D", false);
        fibre_writer.WriteAllOrtho(fibres, second, third);
        
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestFibreWriter/";
        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + results_dir + "/Orthotropic3D.ortho heart/test/data/fibre_tests/Orthotropic3D.ortho").c_str()), 0);
    }
    
    void TestOrthoWriterBinary() throw (Exception)
    {
        FileFinder file_finder("heart/test/data/fibre_tests/Orthotropic3D.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file_finder, ORTHO);
        std::vector< c_vector<double, 3> > fibres;
        std::vector< c_vector<double, 3> > second;
        std::vector< c_vector<double, 3> > third;
        fibre_reader.GetAllOrtho(fibres, second, third);
          
        //Write binary file
        FibreWriter<3> fibre_writer("TestFibreWriter", "Orthotropic3DBin", false);
        fibre_writer.SetWriteFileAsBinary();
        fibre_writer.WriteAllOrtho(fibres, second, third);
        
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestFibreWriter/";
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/Orthotropic3DBin.ortho heart/test/data/fibre_tests/Orthotropic3DBin.ortho").c_str()), 0);      
    }
};


#endif /*TESTFIBREREADER_HPP_*/
