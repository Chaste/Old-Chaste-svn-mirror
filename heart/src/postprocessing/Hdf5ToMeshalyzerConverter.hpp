/*

Copyright (C) University of Oxford, 2005-2009

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

#ifndef HDF5TOMESHALYZERCONVERTER_HPP_
#define HDF5TOMESHALYZERCONVERTER_HPP_

#include <string>
#include "Hdf5DataReader.hpp"

/**
 *  This class converts from Hdf5 format to meshalyzer format, ie, for
 *  voltage, one file, which looks like
 *
 *  V_node_0_time_0
 *  ..
 *  V_node_N_time_0
 *  V_node_0_time_1
 *  ..
 *  V_node_N_time_1
 *  V_node_0_time_2
 *  ..
 *  V_node_N_time_M
 *
 *  The files that are written are [base_name]_V.dat or [base_name]_Phi_e.dat,
 *  where [base_name] is the base name of the original .h5 file. The new files
 *  are written in the same directory as the .h5 file. All paths are relative
 *  to the CHASTE_TEST_OUTPUT directory.
 */
class Hdf5ToMeshalyzerConverter
{
private:
    Hdf5DataReader* mpReader; /**< Pointer to reader of the file to be converted*/
    std::string mOutputDirectory; /**< Directory name where the Meshalyzer output is going*/
    std::string mFileBaseName; /**< Base name for the files [basename]_V.dat etc.*/

    /** A helper method which takes in a string, which must be 'V' or 'Phi_e'
     *  and reads the data corresponding to that string, writing it out in
     *  meshalyzer format.
     * @param type - the type of data stored in this file (V/Phi_e)
     */
    void Write(std::string type);


public:
    /** Constructor, which does the conversion.
     *  @param inputDirectory The input directory, relative to CHASTE_TEST_OUTPUT, where the .h5 file has been written
     *  @param outputDirectory  The output directory, relative to CHASTE_TEST_OUTPUT, where the output will be placed
     *  @param fileBaseName The base name of the data file.
     */
    Hdf5ToMeshalyzerConverter(std::string inputDirectory,
                              std::string outputDirectory,
                              std::string fileBaseName);

    ~Hdf5ToMeshalyzerConverter();
};

#endif /*HDF5TOMESHALYZERCONVERTER_HPP_*/
