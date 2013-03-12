/*

Copyright (c) 2005-2013, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef NAGAIHONDAFORCE_HPP_
#define NAGAIHONDAFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "VertexBasedCellPopulation.hpp"

#include <iostream>

/**
 * A force class for use in vertex-based simulations, based on a mechanical
 * model proposed by T. Nagai and H. Honda ("A dynamic cell model for the formation
 * of epithelial tissues", Philosophical Magazine Part B 81:699-719).
 *
 * Each of the model parameter member variables are rescaled such that mDampingConstantNormal
 * takes the default value 1, whereas Nagai and Honda (who denote the parameter by
 * nu) take the value 0.01.
 */
template<unsigned DIM>
class NagaiHondaForce  : public AbstractForce<DIM>
{
friend class TestForces;

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mNagaiHondaDeformationEnergyParameter;
        archive & mNagaiHondaMembraneSurfaceEnergyParameter;
        archive & mNagaiHondaCellCellAdhesionEnergyParameter;
        archive & mNagaiHondaCellBoundaryAdhesionEnergyParameter;
        archive & mMatureCellTargetArea;
    }

protected:

    /**
     * Cell deformation energy parameter. Has units of kg s^-2 (cell size at equilibrium rest length)^-1.
     */
    double mNagaiHondaDeformationEnergyParameter;

    /**
     * Cell membrane energy parameter. Has units of kg (cell size at equilibrium rest length) s^-2.
     */
    double mNagaiHondaMembraneSurfaceEnergyParameter;

    /**
     * Cell-cell adhesion energy parameter. Has has units of kg (cell size at equilibrium rest length)^2 s^-2.
     */
    double mNagaiHondaCellCellAdhesionEnergyParameter;

    /**
     * Cell-boundary adhesion energy parameter. Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     */
    double mNagaiHondaCellBoundaryAdhesionEnergyParameter;

    /**
     * Non-dimensional target area of a mature (fully-grown) Cell.
     */
    double mMatureCellTargetArea;

public:

    /**
     * Constructor.
     */
    NagaiHondaForce();

    /**
     * Destructor.
     */
    virtual ~NagaiHondaForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the force on each node in the vertex-based cell population based on the
     * Nagai Honda model.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Get the adhesion parameter for the edge between two given nodes.
     *
     * @param pNodeA one node
     * @param pNodeB the other node
     * @param rVertexCellPopulation reference to the cell population
     *
     * @return the adhesion parameter for this edge.
     */
    virtual double GetAdhesionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation);

    /**
     * @return mNagaiHondaDeformationEnergyParameter
     */
    double GetNagaiHondaDeformationEnergyParameter();

    /**
     * @return mNagaiHondaMembraneSurfaceEnergyParameter
     */
    double GetNagaiHondaMembraneSurfaceEnergyParameter();

    /**
     * @return mCellCellAdhesionEnergyParameter
     */
    double GetNagaiHondaCellCellAdhesionEnergyParameter();

    /**
     * @return mNagaiHondaCellBoundaryAdhesionEnergyParameter
     */
    double GetNagaiHondaCellBoundaryAdhesionEnergyParameter();

    /**
     * Set mNagaiHondaDeformationEnergyParameter.
     *
     * @param nagaiHondaDeformationEnergyParameter the new value of mNagaiHondaDeformationEnergyParameter
     */
    void SetNagaiHondaDeformationEnergyParameter(double nagaiHondaDeformationEnergyParameter);

    /**
     * Set mNagaiHondaMembraneSurfaceEnergyParameter.
     *
     * @param nagaiHondaMembraneSurfaceEnergyParameter the new value of mNagaiHondaMembraneSurfaceEnergyParameter
     */
    void SetNagaiHondaMembraneSurfaceEnergyParameter(double nagaiHondaMembraneSurfaceEnergyParameter);

    /**
     * Set mNagaiHondaCellCellAdhesionEnergyParameter.
     *
     * @param nagaiHondaCellCellAdhesionEnergyEnergyParameter the new value of mNagaiHondaCellCellAdhesionEnergyParameter
     */
    void SetNagaiHondaCellCellAdhesionEnergyParameter(double nagaiHondaCellCellAdhesionEnergyEnergyParameter);

    /**
     * Set mNagaiHondaCellBoundaryAdhesionEnergyParameter.
     *
     * @param nagaiHondaCellBoundaryAdhesionEnergyParameter the new value of mNagaiHondaCellBoundaryAdhesionEnergyParameter
     */
    void SetNagaiHondaCellBoundaryAdhesionEnergyParameter(double nagaiHondaCellBoundaryAdhesionEnergyParameter);

    /**
     * Get the target area of a given cell. This grows linearly from
     * 0.5*A to A during the G1 phase of the cell cycle, then remains
     * at A for the rest of the cell cycle, where A denotes the
     * member variable mMatureCellTargetArea.
     *
     * @param pCell the cell
     * @return the cell's target area
     */
    virtual double GetTargetAreaOfCell(const CellPtr pCell) const;

    /**
     * @return mMatureCellTargetArea
     */
    double GetMatureCellTargetArea() const;

    /**
     * Set mMatureCellTargetArea.
     *
     * @param matureCellTargetArea the new value of mMatureCellTargetArea
     */
    void SetMatureCellTargetArea(double matureCellTargetArea);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NagaiHondaForce)

#endif /*NAGAIHONDAFORCE_HPP_*/
