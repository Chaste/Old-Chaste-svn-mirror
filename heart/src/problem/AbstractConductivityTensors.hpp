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
#ifndef ABSTRACTCONDUCTIVITYTENSORS_HPP_
#define ABSTRACTCONDUCTIVITYTENSORS_HPP_

#include <vector>
#include <string>
#include <fstream>
#include "UblasIncludes.hpp"
#include "Exception.hpp"

template<unsigned SPACE_DIM>
class AbstractConductivityTensors
{
protected:
    unsigned mNumElements;

    bool mUseNonConstantConductivities;
    bool mUseFibreOrientation;

    // Constant Conductivities
    c_vector<double, SPACE_DIM> mConstantConductivities; // mS/cm

    // Non-constant Conductivities
    std::vector<c_vector<double, SPACE_DIM> >* mpNonConstantConductivities; // mS/cm

    // Container for Tensors (single or multiple)
    std::vector< c_matrix<double,SPACE_DIM,SPACE_DIM> > mTensors;

    bool mInitialised;

    std::string mFibreOrientationFilename;
    std::ifstream mDataFile;


    void OpenFibreOrientationFile();

    void CloseFibreOrientationFile();

    unsigned GetTokensAtNextLine(std::vector<double>& tokens);

    unsigned GetNumElementsFromFile();

public:

    AbstractConductivityTensors();

    virtual ~AbstractConductivityTensors();

    /**
     *  Sets a file for reading the fibre orientation from.
     *
     *  @param rFibreOrientationFilename Relative path to the file defining the fibre orientation
     */
    void SetFibreOrientationFile(const std::string &rFibreOrientationFilename);

    /**
     *  Sets constant conductivities for all the elements of the mesh.
     *
     *  @param constLongConduc Longitudinal conductivity (x axis)
     *  @param constTransConduc Transverse conductivity (y axis)
     *  @param constNormalConduc Normal conductivity (z axis)
     *
     *  \todo Explain problem with c_vector and SPACE_DIM
     */
    void SetConstantConductivities(c_vector<double, 1> constantConductivities);

    void SetConstantConductivities(c_vector<double, 2> constantConductivities);

    virtual void SetConstantConductivities(c_vector<double, 3> constantConductivities);


    /**
     *  Sets a different longitudinal and transverse conductivity for every elements of the mesh.
     *
     *  @param rLongitudinalConductivities Vector containing longitudinal conductivities of the elements (x axis)
     *  @param rTransverseConductivities Vector containing transverse conductivities of the elements (y axis)
     *  @param rNormalConductivities Vector containing normal conductivities of the elements (z axis)
     */
    void SetNonConstantConductivities(std::vector<c_vector<double, SPACE_DIM> >* pNonConstantConductivities);

    /**
     *  Computes the tensors based in all the info set
     */
    virtual void Init() throw (Exception) = 0;

    /**
     *  Returns the diffussion tensor of the element number "index"
     *
     *  @param index Index of the element of the mesh
     */
    c_matrix<double,SPACE_DIM,SPACE_DIM>& operator[](const unsigned index);
};



#endif /*ABSTRACTCONDUCTIVITYTENSORS_HPP_*/
