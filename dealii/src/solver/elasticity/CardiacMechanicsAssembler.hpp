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


#ifndef CARDIACMECHANICSASSEMBLER_HPP_
#define CARDIACMECHANICSASSEMBLER_HPP_

#include "FiniteElasticityAssembler.cpp"
#include "AbstractCardiacMechanicsAssembler.hpp"

template<unsigned DIM>
class CardiacMechanicsAssembler : public FiniteElasticityAssembler<DIM>, public AbstractCardiacMechanicsAssembler<DIM>
{
friend class TestImplicitCardiacMechanicsAssembler2;

protected:
    bool mAllocatedMaterialLawMemory;

    /**
     *  Storage space for dTdE when T and E are in the rotated fibre-sheet frame
     */
    FourthOrderTensor<DIM>  mDTdE_fibre;

    /**
     *  The matrix P using JonW's convention. Orthogonal
     */
    Tensor<2,DIM> mFibreSheetMat;

    /**
     *  The transpose of P, which is also the inverse of P
     */
    Tensor<2,DIM> mTransFibreSheetMat;

    /**
     *  The active tension, quad point-wise. NOTE: the i-th entry of this vector is
     *  assumed to be the i-th quad point obtained by looping over cells in the obvious
     *  way and then looping over quad points
     */
    std::vector<double> mActiveTension;

    /** A scale factor by which (dimensional) material parameters are scaled (which requires
     *  the relevant forcing terms, in this case active tension to be scaled likewise). Defaults
     *  to 1
     */
    double mScaleFactor;


    /** Overloaded method for assembling system, which takes into account the active tensions */
    void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter,
                           Vector<double>&       elementRhs,
                           FullMatrix<double>&   elementMatrix,
                           bool                  assembleResidual,
                           bool                  assembleJacobian);


public:
    /**
     *  Constructor
     *
     *  @param pMesh. A pointer to the mesh. Should have a surface set as the fixed surface
     *  @param outputDirectory. The output directory, relative to TEST_OUTPUT
     *  @param pMaterialLaw. The material law for the tissue. Defaults to NULL, in which case
     *   a default material law is used.
     */
    CardiacMechanicsAssembler(Triangulation<DIM>* pMesh,
                              std::string outputDirectory,
                              AbstractIncompressibleMaterialLaw2<DIM>* pMaterialLaw = NULL);
    virtual ~CardiacMechanicsAssembler();


    /**
     *  Specify a constant fibre-sheet rotation matrix
     *
     *  This is really a temporary method until the fibre-sheet direction can be read in
     */
    virtual void SetFibreSheetMatrix(Tensor<2,DIM> fibreSheetMat);

    virtual void Solve(double currentTime, double nextTime, double timestep);

    /**
     *  Set the current active tensions, by quadrature point. Quad points don't have indices,
     *  so these values should be in the order given by looping over cells and then looping
     *  over quad points
     */
    virtual void SetForcingQuantity(std::vector<double>& activeTension);

    /**
     *  Set a scale factor by which (dimensional) material parameters are scaled. For
     *  this assembler the active tension to be scaled likewise when it is used. A scale
     *  factor may be used/needed to improve GMRES convergence.
     */
    void SetScaling(double scaleFactor);
};

#endif /*CARDIACMECHANICSASSEMBLER_HPP_*/
