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

#ifndef ELECTROMECHANICSPROBLEMDEFINITION_HPP_
#define ELECTROMECHANICSPROBLEMDEFINITION_HPP_

#include "SolidMechanicsProblemDefinition.hpp"
#include "ContractionModelName.hpp"
#include "NashHunterPoleZeroLaw.hpp"
#include "CompressibleExponentialLaw.hpp"



/**
 *  Subclass of SolidMechanicsProblemDefinition with some cardiac-electro-mechanics-specific
 *  methods.
 */
template<unsigned DIM>
class ElectroMechanicsProblemDefinition : public SolidMechanicsProblemDefinition<DIM>
{
private:
    /**
     *  The contraction model used (ContractionModelName is an enumeration containing all contraction
     *  models implemented.
     */
    ContractionModelName mContractionModel;

    /** Timestep to use when solving contraction models */
    double mContractionModelOdeTimeStep;

    /** How often a mechanics solve should be done */
    double mMechanicsSolveTimestep;

    /**
     *  Whether the deformation should affect the electrical physiological conductivity
     *  (or whether this effect is neglected)
     */
    bool mDeformationAffectsConductivity;

    /**
     *  Whether the deformation should affect the cardiac cell models, for example if there
     *  are stretch-activated channels in the cell model.
     */
    bool mDeformationAffectsCellModels;

    /**
     *  This member variable is used if SetDefaultCardiacMateriawLaw() is called.
     */
    AbstractMaterialLaw<DIM>* mpDefaultMaterialLaw;

public:
    /**
     * Constructor
     * @param rMesh the mesh
     */
    ElectroMechanicsProblemDefinition(QuadraticMesh<DIM>& rMesh);

    /** Destructor */
    ~ElectroMechanicsProblemDefinition();

    /**
     * Set the contraction model to be used (throughout the tissue).
     * @param contractionModel contraction model (from the enumeration ContractionModelName)
     * @param timestep timestep to be used in solving (ODE-based) contraction models.
     */
    void SetContractionModel(ContractionModelName contractionModel, double timestep);

    /**
     * Use the default material law (NashHunter in the incompressible case, exponential in the
     *  compressible case), throughout the tissue.
     * @param compressibilityType Either INCOMPRESSIBLE or COMPRESSIBLE
     */
    void SetUseDefaultCardiacMaterialLaw(CompressibilityType compressibilityType);

    /**
     * Set how the deformation should affect the electro-physiology
     * @param deformationAffectsConductivity Whether the deformation should affect the electrical
     *   physiological conductivity (or whether this effect is neglected)
     * @param deformationAffectsCellModels Whether the deformation should affect the cardiac cell
     *   models, for example if there are stretch-activated channels in the cell model.
     */
    void SetDeformationAffectsElectrophysiology(bool deformationAffectsConductivity, bool deformationAffectsCellModels);

    /**
     *  Set how often the mechanics is solved for.
     *  @param timestep timestep
     */
    void SetMechanicsSolveTimestep(double timestep);

    /**
     *  Get the contraction model
     */
    ContractionModelName GetContractionModel()
    {
        if(mContractionModelOdeTimeStep < 0.0)
        {
            EXCEPTION("Contraction model hasn't been set yet");
        }
        return mContractionModel;
    }

    /**
     *  Get the contraction model timestep
     */
    double GetContractionModelOdeTimestep()
    {
        if(mContractionModelOdeTimeStep < 0.0)
        {
            EXCEPTION("Contraction model hasn't been set yet");
        }
        return mContractionModelOdeTimeStep;
    }

    /**
     *  Get how often the mechanics is solved
     */
    double GetMechanicsSolveTimestep()
    {
        if(mMechanicsSolveTimestep < 0.0)
        {
            EXCEPTION("Timestep for mechanics solve hasn't been set yet");
        }
        return mMechanicsSolveTimestep;
    }

    /**
     *  Get whether the deformation affects the electrical physiological conductivity
     *  (or whether this effect is neglected).
     */
    bool GetDeformationAffectsConductivity()
    {
        return mDeformationAffectsConductivity;
    }

    /**
     *  Get whether the deformation affects the cardiac cell models, for example if
     *  there are stretch-activated channels in the cell model.
     */
    bool GetDeformationAffectsCellModels()
    {
        return mDeformationAffectsCellModels;
    }
};

#endif // ELECTROMECHANICSPROBLEMDEFINITION_HPP_
