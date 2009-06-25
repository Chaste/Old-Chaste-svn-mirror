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
#ifndef AXISYMMETRICCONDUCTIVITYTENSORS_HPP_
#define AXISYMMETRICCONDUCTIVITYTENSORS_HPP_

#include "AbstractConductivityTensors.hpp"

/**
 * The class is templated over SPACE_DIM to keep compatibility with the abstract class.
 * However axisymmetric conductivity only makes sense in 3D, so we check in the constructor
 * for SPACE_DIM to be 3.
 *
 * \todo is the compatibility argument above really worthwhile?  We could just inherit from
 *       AbstractConductivityTensors<3>.  This would however require changes to the
 *       constructors of AbstractCardiacPde and BidomainPde.
 */
template<unsigned SPACE_DIM>
class AxisymmetricConductivityTensors : public AbstractConductivityTensors<SPACE_DIM>
{
private:

    /**
     * Read from mDataFile with GetTokensAtNextLine
     * @param rOrientVector  vector into which to read the orientation
     */
    void ReadOrientationVectorFromFile (c_vector<double,SPACE_DIM>& rOrientVector);

public:
    /**Constructor*/
    AxisymmetricConductivityTensors();

    /**
     *  Sets constant conductivities for all the elements of the mesh.
     *  @param constantConductivities Longitudinal, Transverse (y axis) and Normal conductivity (z axis)
     */
     void SetConstantConductivities(c_vector<double, 3> constantConductivities);

    /** \todo there are extensive comments within the implementation of this method */
    void Init() throw (Exception);
};


#endif /*AXISYMMETRICCONDUCTIVITYTENSORS_HPP_*/
