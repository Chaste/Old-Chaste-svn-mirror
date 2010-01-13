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
#ifndef ORTHOTROPICCONDUCTIVITYTENSORS_HPP_
#define ORTHOTROPICCONDUCTIVITYTENSORS_HPP_

#include <vector>
#include "UblasIncludes.hpp"
#include "AbstractConductivityTensors.hpp"


/**
 *
 *  This class provides an abstraction for the definition of constant/non-constant difussion tensors
 * associated to the different elements of the mesh.
 *
 *  After instanciating the class any of SetFibreOrientationFile() or SetNonConstantConductivities()
 * (or both) can be called to implement fiber orientation or heterogeneous conductivity into the
 * tensors, respectively. If none of them is called a constant tensor (with constant conductivies)
 * will be generated for all the elements of the mesh.
 *
 *  Init() should be called to actually create the tensors.
 *
 *  Initial values for conductivity from "Laminar Arrangement of Ventricular Myocytes Influences Electrical
 * Behavior of the Heart", Hooks et al. 2007
 *
 */
template<unsigned SPACE_DIM>
class OrthotropicConductivityTensors : public AbstractConductivityTensors<SPACE_DIM>
{
private:
    /**
     * Reads a line from the Conductivity Tensor file (see base class)
     * and interprets it as a SPACE_DIM*SPACE_DIM tensor
     * @param orientMatrix matrix into which to read the tensor
     */
    void ReadOrientationMatrixFromFile (c_matrix<double,SPACE_DIM,SPACE_DIM>& orientMatrix);

public:

    /**
     *  Computes the tensors based in all the info set
     */
    void Init() throw (Exception);

};

#endif /*ORTHOTROPICCONDUCTIVITYTENSORS_HPP_*/
