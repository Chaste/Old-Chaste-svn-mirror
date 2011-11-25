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


#ifndef SOLIDMECHANICSPROBLEMDEFINITION_HPP_
#define SOLIDMECHANICSPROBLEMDEFINITION_HPP_

#include "QuadraticMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "AbstractIncompressibleMaterialLaw.hpp"
#include "AbstractCompressibleMaterialLaw.hpp"
#include "CompressibilityType.hpp"




/**
 *  Simple enumeration for denoting the type of body force, constant or defined via a function
 */
typedef enum BodyForceType_
{
    CONSTANT_BODY_FORCE,
    FUNCTIONAL_BODY_FORCE
} BodyForceType;

/**
 *  Simple enumeration for denoting the type of traction (Neumann) boundary condition:
 */
typedef enum TractionBoundaryConditionType_
{
    NO_TRACTIONS,
    ELEMENTWISE_TRACTION,
    FUNCTIONAL_TRACTION,
    PRESSURE_ON_DEFORMED
} TractionBoundaryConditionType;


/**
 *  A class for specifying various parts of a solid mechanics problem, in particular the material
 *  laws for the deforming body, fixed nodes information, the body force (per unit mass)
 *  (usually acceleration due to gravity or zero), the traction boundary conditions, and the density.
 */
template<unsigned DIM>
class SolidMechanicsProblemDefinition
{
public:
    /** Special displacement that indicates a node is free to move (in a particular direction).
     *  Used in SetFixedNodes(). */
	static const double FREE; // set to DBL_MAX

private:
    /** The mesh being solved on */
    QuadraticMesh<DIM>& mrMesh;

    /** Density of the body (constant throughout body) */
    double mDensity;

    //////////////////////////////
    // body force
    //////////////////////////////

    /** The body force type */
    BodyForceType mBodyForceType;

    /** The constant body force, only used if mBodyForceType is set appropriately */
    c_vector<double,DIM> mConstantBodyForce;

    /** The body force as a function of space and time, only used if mBodyForceType is set appropriately */
    c_vector<double,DIM> (*mpBodyForceFunction)(c_vector<double,DIM>& rX, double t);

    //////////////////////////////
    // tractions
    //////////////////////////////

    /** The traction (Neumann) boundary condition type */
    TractionBoundaryConditionType mTractionBoundaryConditionType;

    /** The surface elements on which tractions are applied */
    std::vector<BoundaryElement<DIM-1,DIM>*> mTractionBoundaryElements;

    /** The tractions on each surface element (only used if mTractionBoundaryConditionType is set appropriately) */
    std::vector<c_vector<double,DIM> > mElementwiseTractions;

    /** If the tractions are specified to correspond to a pressure acting on the surface: the pressure for the given
     *  boundary elements (only used if mTractionBoundaryConditionType is set appropriately) */
    double mNormalPressure;

    /** The tractions as a function of space and time (only used if mTractionBoundaryConditionType is set appropriately) */
    c_vector<double,DIM> (*mpTractionBoundaryConditionFunction)(c_vector<double,DIM>& rX, double t);

    //////////////////////////////
    // fixed nodes
    //////////////////////////////

    /** All nodes (including non-vertices) which are fixed. */
    std::vector<unsigned> mFixedNodes;

    /** The displacements of those nodes with displacement boundary conditions. */
    std::vector<c_vector<double,DIM> > mFixedNodeDisplacements;

    /////////////////////////////
    // material law
    /////////////////////////////
    /**
     *  The material law, in the case of incompressible material laws. This vector is either of size 1, representing a
     *  homogeneous material, or of size num_elements, representing a heterogeneous material, with a material law per element.
     *  If he material is compressible, this vector will be of size zero.
     */
    std::vector<AbstractIncompressibleMaterialLaw<DIM>*> mIncompressibleMaterialLaws;
    /**
     *  The material law, in the case of compressible material laws. This vector is either of size 1, representing a
     *  homogeneous material, or of size num_elements, representing a heterogeneous material, with a material law per element.
     *  If the material is incompressible, this vector will be of size zero.
     */
    std::vector<AbstractCompressibleMaterialLaw<DIM>*> mCompressibleMaterialLaws;

    /** Whether the material is homogeneous (same material law everywhere) or heterogeneous */
    bool mIsHomogeneousMaterial;

    /** Whether the material is incompressible or compressible. (CompressibilityType is an enumeration).  */
    CompressibilityType mCompressibilityType;

    /**
     *  Helper function for checking whether a dynamic_cast succeeded or not, and throwing an exception
     *  if it failed.
     *
     *  @param compressibilityType compressibility type
     *  @param pMaterialLaw material law
     */
    void CheckCastSuccess(CompressibilityType compressibilityType,AbstractMaterialLaw<DIM>* pMaterialLaw);

public:
    /**
     * Constructor initialises the body force to zero and density to 1.0
     * @param rMesh  is the mesh being solved on
     */
    SolidMechanicsProblemDefinition(QuadraticMesh<DIM>& rMesh);

    /**
     *  Set the density
     *  @param density
     */
    void SetDensity(double density);

    /**
     * Get the density
     */
    double GetDensity();

    /**
     * Set a list of nodes (indices) to be fixed in space with zero displacement
     * @param rFixedNodes the fixed nodes
     */
    void SetZeroDisplacementNodes(std::vector<unsigned>& rFixedNodes);

    /**
     * Set a list of nodes to be fixed, with their corresponding new LOCATIONS (not displacements)
     * @param rFixedNodes the fixed node indices
     * @param rFixedNodeLocation corresponding locations
     */
    void SetFixedNodes(std::vector<unsigned>& rFixedNodes, std::vector<c_vector<double,DIM> >& rFixedNodeLocation);

    /**
     *  Get the fixed nodes
     */
    std::vector<unsigned>& rGetFixedNodes();

    /**
     * Get the fixed node displacements
     */
    std::vector<c_vector<double,DIM> >& rGetFixedNodeDisplacements();

    /**
     * Set the body force to be used - this version sets a constant body force
     * @param bodyForce the constant body force
     */
    void SetBodyForce(c_vector<double,DIM> bodyForce);

    /**
     * Set the body force to be used - this version sets a functional body force
     * @param pFunction a function of space and time returning a vector
     */
    void SetBodyForce(c_vector<double,DIM> (*pFunction)(c_vector<double,DIM>& rX, double t));

    /**
     * Get the body force type
     */
    BodyForceType GetBodyForceType();

    /**
     * Get the constant body force (error if GetBodyForceType()!=CONSTANT_BODY_FORCE)
     */
    c_vector<double,DIM> GetConstantBodyForce();

    /**
     * Evaluate the body force function (error if GetBodyForceType()!=FUNCTIONAL_BODY_FORCE)
     *
     * @param rX spatial location
     * @param t current time
     */
    c_vector<double,DIM> EvaluateBodyForceFunction(c_vector<double,DIM>& rX, double t);

    /**
     * Get the traction (Neumann) boundary condition type
     */
    TractionBoundaryConditionType GetTractionBoundaryConditionType();

    /**
     * Set traction (Neumann) boundary conditions. This version takes in a vector of
     *  boundary elements, and corresponding tractions for each one.
     * @param rTractionBoundaryElements the boundary elements
     * @param rElementwiseTractions corresponding tractions
     */
    void SetTractionBoundaryConditions(std::vector<BoundaryElement<DIM-1,DIM>*>& rTractionBoundaryElements,
                                       std::vector<c_vector<double,DIM> >& rElementwiseTractions);

    /**
     * Set traction (Neumann) boundary conditions. This version takes in a vector of
     *  boundary elements, and a function to be evaluated at points in these boundary
     *  elements
     * @param rTractionBoundaryElements the boundary elements
     * @param pFunction the traction function (a function of space and time, returning a vector)
     */
    void SetTractionBoundaryConditions(std::vector<BoundaryElement<DIM-1,DIM>*> rTractionBoundaryElements,
                                       c_vector<double,DIM> (*pFunction)(c_vector<double,DIM>& rX, double t));

    /**
     * Set traction (Neumann) boundary conditions. This version says that pressures should be applied
     * on surfaces in the DEFORMED body in the outward normal direction.
     *
     * @param rTractionBoundaryElements The boundary elements
     * @param normalPressure the corresponding pressure
     */
    void SetApplyNormalPressureOnDeformedSurface(std::vector<BoundaryElement<DIM-1,DIM>*> rTractionBoundaryElements,
                                                 double normalPressure);

    /**
     * Get the vector of traction boundary elements
     */
    std::vector<BoundaryElement<DIM-1,DIM>*>& rGetTractionBoundaryElements();

    /**
     *  Get the element-wise tractions vector (corresponding to
     *  vector returned by rGetTractionBoundaryElements())
     *  (error if GetTractionBoundaryConditionType()!=ELEMENTWISE_TRACTION)
     */
    std::vector<c_vector<double,DIM> >& rGetElementwiseTractions();

    /**
     *  Get the pressure for the boundary elements (corresponding to
     *  vector returned by rGetTractionBoundaryElements())
     *  (error if GetTractionBoundaryConditionType()!=PRESSURE_ON_DEFORMED)
     */
    double GetNormalPressure();


    /**
     * Evaluate the traction boundary condition function (error if GetTractionBoundaryConditionType()!=FUNCTIONAL_TRACTION)
     *
     * @param rX spatial location
     * @param t current time
     */
    c_vector<double,DIM> EvaluateTractionFunction(c_vector<double,DIM>& rX, double t);

    /**
     * Set a material law for the entire body (ie the homogeneous case). If compressibilityType==INCOMPRESSIBLE,
     * the material law pointer will be checked at run-time that it is of type `AbstractIncompressibleMaterialLaw`,
     * and similarly for the compressible case. Any previous material information will be deleted.
     *
     * @param compressibilityType either 'INCOMPRESSIBLE' or 'COMPRESSIBLE'
     * @param pMaterialLaw The material law for the entire body
     */
    void SetMaterialLaw(CompressibilityType compressibilityType, AbstractMaterialLaw<DIM>* pMaterialLaw);

    /**
     * Set a vector of material laws for the body, one for each element in the mesh (the heterogeneous case). If
     * compressibilityType==INCOMPRESSIBLE, the material law pointer will be checked at run-time that it is
     * of type `AbstractIncompressibleMaterialLaw`, and similarly for the compressible case.
     * Any previous material information will be deleted.
     *
     * @param compressibilityType either 'INCOMPRESSIBLE' or 'COMPRESSIBLE'
     * @param rMaterialLaws Vector of pointers to material laws
     */
    void SetMaterialLaw(CompressibilityType compressibilityType, std::vector<AbstractMaterialLaw<DIM>*>& rMaterialLaws);

    /**
     * Get whether the material is homogeneous or heterogeneous.
     * SetMaterialLaw() must be called before calling this.
     */
    bool IsHomogeneousMaterial();

    /**
     * Get whether the material is incompressible or compressible.
     * SetMaterialLaw() must be called before calling this.
     */
    CompressibilityType GetCompressibilityType();

    /**
     * Get the material law for a given element, when the body is incompressible. An assertion will
     * fail if GetCompressibilityType()!=INCOMPRESSIBLE. If the material is homogeneous, it doesn't matter what
     * the element index is.
     * @param elementIndex index of element
     */
    AbstractIncompressibleMaterialLaw<DIM>* GetIncompressibleMaterialLaw(unsigned elementIndex);

    /**
     * Get the material law for a given element, when the body is compressible. An assertion will
     * fail if GetCompressibilityType()!=COMPRESSIBLE. If the material is homogeneous, it doesn't matter what
     * the element index is.
     * @param elementIndex index of element
     */
    AbstractCompressibleMaterialLaw<DIM>* GetCompressibleMaterialLaw(unsigned elementIndex);
};




#endif /* SOLIDMECHANICSPROBLEMDEFINITION_HPP_ */
