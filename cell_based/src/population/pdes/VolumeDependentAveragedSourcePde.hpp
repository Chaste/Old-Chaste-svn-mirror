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

#ifndef VOLUMEDEPENDENTAVERAGEDSOURCEPDE_HPP_
#define VOLUMEDEPENDENTAVERAGEDSOURCEPDE_HPP_

#include "NodeBasedCellPopulation.hpp"
#include "AveragedSourcePde.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractLinearEllipticPde.hpp"

/**
 *  A PDE which calculates the source term by adding the number of cells
 *  in the element containing that point and scaling by the element area.
 */

template<unsigned DIM>
class VolumeDependentAveragedSourcePde : public AveragedSourcePde<DIM>
{
private:

    /** The cell population member. */
    AbstractCellPopulation<DIM>& mrCellPopulation;

    /** Static cast of the NodeBasedCellPopulation **/
    NodeBasedCellPopulation<DIM>* mpStaticCastCellPopulation;

    /** Coefficient of consumption of nutrient by cells per unit volume. */
    double mCoefficient;

    /** Vector of averaged cell densities on elements of the coarse mesh. */
    std::vector<double> mCellDensityOnCoarseElements;

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation reference to the cell population
     * @param coefficient the coefficient of consumption of nutrient by cells
     */
    VolumeDependentAveragedSourcePde(AbstractCellPopulation<DIM>& rCellPopulation, double coefficient);

    /**
     * Set up the source terms.
     *
     * @param rCoarseMesh reference to the coarse mesh
     * @param pCellPdeElementMap pointer to the map from cells to coarse elements
     */
    void SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh, std::map< CellPtr, unsigned >* pCellPdeElementMap=NULL);

    /**
     * Overridden ComputeConstantInUSourceTerm() method.
     *
     * @param rX The point in space
     * @param pElement The element
     *
     * @return the constant in u part of the source term, i.e g(x) in
     *  Div(D Grad u)  +  f(x)u + g(x) = 0.
     */
    double ComputeConstantInUSourceTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement);

    /**
     * Overridden ComputeLinearInUCoeffInSourceTerm() method.
     *
     * @param rX The point in space
     * @param pElement the element
     *
     * @return the coefficient of u in the linear part of the source term, i.e f(x) in
     *  Div(D Grad u)  +  f(x)u + g(x) = 0.
     */
    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement);

    /**
     * Overridden ComputeDiffusionTerm() method.
     *
     * @param rX The point in space at which the diffusion term is computed
     *
     * @return a matrix.
     */
    c_matrix<double,DIM,DIM> ComputeDiffusionTerm(const ChastePoint<DIM>& rX);

    /**
     * Returns the uptake rate.
     *
     * @param elementIndex the element we wish to return the uptake rate for
     */
    double GetUptakeRateForElement(unsigned elementIndex);
};

#endif /*VOLUMEDEPENDENTAVERAGEDSOURCEPDE_HPP_*/
