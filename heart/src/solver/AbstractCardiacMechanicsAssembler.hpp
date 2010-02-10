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

#ifndef ABSTRACTCARDIACMECHANICSASSEMBLER_HPP_
#define ABSTRACTCARDIACMECHANICSASSEMBLER_HPP_

#include "NonlinearElasticityAssembler.hpp"
#include "NashHunterPoleZeroLaw.hpp"
#include "QuadraticBasisFunction.hpp" // not included in NonlinearElasticityAssembler.hpp, just the cpp
#include "LinearBasisFunction.hpp"
#include "AbstractContractionModel.hpp"

/**
 *  AbstractCardiacMechanicsAssembler
 * 
 *  Base class to implicit and explicit cardiac mechanics assemblers. Inherits from NonlinearElasticityAssembler
 *  Main method is the overloaded AssembleOnElement which does the extra work needed for cardiac problems. The
 *  child classes hold the contraction models and need to implement a method for getting the active tension from 
 *  the model. 
 */
template<unsigned DIM>
class AbstractCardiacMechanicsAssembler : public NonlinearElasticityAssembler<DIM>
{
protected:
    static const unsigned STENCIL_SIZE = NonlinearElasticityAssembler<DIM>::STENCIL_SIZE;
    static const unsigned NUM_NODES_PER_ELEMENT = NonlinearElasticityAssembler<DIM>::NUM_NODES_PER_ELEMENT;
    static const unsigned NUM_VERTICES_PER_ELEMENT = NonlinearElasticityAssembler<DIM>::NUM_VERTICES_PER_ELEMENT;

    /**
     *  Vector of contraction model (pointers). One for each quadrature point.
     *  Note the indexing: the i-th entry corresponds to the i-th global quad 
     *  point, when looping over elements and then quad points 
     */
    std::vector<AbstractContractionModel*> mContractionModelSystems;

    /**
     *  Stored stretches (in fibre direction, at each quadrature point). Should be stored
     *  when GetActiveTensionAndTensionDerivs() is called, and can be used either in that
     *  timestep (implicit solver), or the next timestep (explicit solver)
     */
    std::vector<double> mStretches;    

    /** Total number of quad points in the (mechanics) mesh */
    unsigned mTotalQuadPoints;

    /** Whether the material law was passed in or the default used */
    bool mAllocatedMaterialLawMemory;
    
    /** Current time */
    double mCurrentTime;
    /** Time to which the solver has been asked to solve to */
    double mNextTime;
    /** Time used to integrate the contraction model */
    double mOdeTimestep;
    
    /** The fibre-sheet-normal directions (in a matrix), if constant (defaults to the identity, ie fibres in the X-direction, sheet in the XY plane) */
    c_matrix<double,DIM,DIM> mConstantFibreSheetDirections;

    /**
	 *	The fibre-sheet-normal directions (matrices), one for each element. Only non-NULL if SetVariableFibreSheetDirections()
     *  is called, if not mConstantFibreSheetDirections is used instead
     */
    std::vector<c_matrix<double,DIM,DIM> >* mpVariableFibreSheetDirections;
    
    
    /**
     *  Whether the solver is implicit or not (ie whether the contraction model depends on lambda (and depends on
     *  lambda at the current time)). For whether dTa_dLam dependent terms need to be added to the Jacbobian
     */
    virtual bool IsImplicitSolver()=0;

    /**
     * Overloaded AssembleOnElement. Apart from a tiny bit of initial set up and
     * the lack of the body force term in the residual, the bits where this is
     * different to the base class AssembleOnElement are restricted to two bits
     * (see code): getting Ta and using it to compute the stress, and (in when Ta
     * is a function of the stretch) the addition of extra term to the Jacobian.
     * 
     * @param rElement The element to assemble on.
     * @param rAElem The element's contribution to the LHS matrix is returned in this
     *     n by n matrix, where n is the no. of nodes in this element. There is no
     *     need to zero this matrix before calling.
     * @param rAElemPrecond The element's contribution to the matrix passed to PetSC
     *     in creating a preconditioner
     * @param rBElem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     * @param assembleResidual A bool stating whether to assemble the residual vector.
     * @param assembleJacobian A bool stating whether to assemble the Jacobian matrix.
     */
    void AssembleOnElement(Element<DIM, DIM>& rElement,
                           c_matrix<double,STENCIL_SIZE,STENCIL_SIZE>& rAElem,
                           c_matrix<double,STENCIL_SIZE,STENCIL_SIZE>& rAElemPrecond,
                           c_vector<double,STENCIL_SIZE>& rBElem,
                           bool assembleResidual,
                           bool assembleJacobian);


    /**
     *  Pure method called in AbstractCardiacMechanicsAssembler::AssembleOnElement(), which needs to provide
     *  the active tension (and other info if implicit (if the contraction model depends on stretch
     *  or stretch rate)) at a particular quadrature point. Takes in the current fibre stretch.
     * 
     *  @param currentFibreStretch The stretch in the fibre direction
     *  @param currentQuadPointGlobalIndex quadrature point integrand currently being evaluated at in AssembleOnElement
     *  @param assembleJacobian  A bool stating whether to assemble the Jacobian matrix.
     *  @param rActiveTension The returned active tension
     *  @param rDerivActiveTensionWrtLambda The returned dT_dLam, derivative of active tension wrt stretch. Only should be set in implicit assemblers
     *  @param rDerivActiveTensionWrtDLambdaDt The returned dT_dLamDot, derivative of active tension wrt stretch rate.  Only should be set in implicit assemblers
     */
    virtual void GetActiveTensionAndTensionDerivs(double currentFibreStretch, 
                                                  unsigned currentQuadPointGlobalIndex,
                                                  bool assembleJacobian,
                                                  double& rActiveTension,
                                                  double& rDerivActiveTensionWrtLambda,
                                                  double& rDerivActiveTensionWrtDLambdaDt)=0;
public:
    /**
     * Constructor
     *
     * @param pQuadMesh A pointer to the mesh.
     * @param outputDirectory The output directory, relative to TEST_OUTPUT
     * @param rFixedNodes The fixed nodes
     * @param pMaterialLaw The material law for the tissue. If NULL the default
     *   (NashHunterPoleZero) law is used.
     */
    AbstractCardiacMechanicsAssembler(QuadraticMesh<DIM>* pQuadMesh,
                                      std::string outputDirectory,
                                      std::vector<unsigned>& rFixedNodes,
                                      AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw)
        : NonlinearElasticityAssembler<DIM>(pQuadMesh,
                                        pMaterialLaw!=NULL ? pMaterialLaw : new NashHunterPoleZeroLaw<DIM>,
                                        zero_vector<double>(DIM),
                                        DOUBLE_UNSET,
                                        outputDirectory,
                                        rFixedNodes),
        mCurrentTime(DBL_MAX),
        mNextTime(DBL_MAX),
        mOdeTimestep(DBL_MAX)
    {
        // compute total num quad points
        mTotalQuadPoints = pQuadMesh->GetNumElements()*this->mpQuadratureRule->GetNumQuadPoints();

        // note that if pMaterialLaw is NULL a new NashHunter law was sent to the
        // NonlinElas constuctor (see above)
        mAllocatedMaterialLawMemory = (pMaterialLaw==NULL);
        
        // initialise the store of fibre stretches
        mStretches.resize(mTotalQuadPoints, 1.0);
        
		// initialise fibre/sheet direction matrix to be the identity, fibres in X-direction, and sheet in XY-plane
        mConstantFibreSheetDirections = zero_matrix<double>(DIM,DIM);
        for(unsigned i=0; i<DIM; i++)
        {
            mConstantFibreSheetDirections(i,i) = 1.0;
        }
        
        mpVariableFibreSheetDirections = NULL;
    }

    /**
     *  Destructor just deletes memory if it was allocated
     */
    ~AbstractCardiacMechanicsAssembler()
    {
        if(mAllocatedMaterialLawMemory)
        {
            assert(this->mMaterialLaws.size()==1); // haven't implemented heterogeniety yet
            delete this->mMaterialLaws[0];
        }

        if(mpVariableFibreSheetDirections)
        {
            delete mpVariableFibreSheetDirections;
        }
    }


    /** Get the total number of quad points in the mesh. Pure, implemented in concrete assembler */
    unsigned GetTotalNumQuadPoints()
    {
        return mTotalQuadPoints;
    }

    /** Get the quadrature rule used in the elements. */
    virtual GaussianQuadratureRule<DIM>* GetQuadratureRule()
    {
        return this->mpQuadratureRule;
    }
    
    
    /**  
	 *	Set a constant fibre-sheet-normal direction (a matrix) to something other than the default (fibres in X-direction, 
	 *  sheet in the XY plane) 
     *  @param rFibreSheetMatrix The fibre-sheet-normal matrix (fibre dir the first column, normal-to-fibre-in sheet in second
     *  column, sheet-normal in third column).
     */
    void SetConstantFibreSheetDirections(const c_matrix<double,DIM,DIM>& rFibreSheetMatrix)
    {
        mConstantFibreSheetDirections = rFibreSheetMatrix; 
// EMTODO:
// CheckOrthogonality(mConstantFibreSheetDirections);
    }
    
    /** 
	 *	Set a variable fibre-sheet-normal direction (matrices), one for each element, from a file.
     *  The file should be a .ortho file (ie each line has the fibre dir, sheet dir, normal dir for that element).
     *  The number of elements must match the number in the MECHANICS mesh!
     *  @param orthoFile the file containing the fibre/sheet directions
     */
    void SetVariableFibreSheetDirections(std::string orthoFile);
        


    /**
     *  Set the intracellular Calcium concentrations and voltages at each quad point. Pure.
     * 
     *  Implicit solvers (for contraction models which are functions of stretch (and maybe 
     *  stretch rate) would integrate the contraction model with this Ca/V/t using the current
     *  stretch (ie inside AssembleOnElement, ie inside GetActiveTensionAndTensionDerivs).
     *  Explicit solvers (for contraction models which are NOT functions of stretch can immediately
     *  integrate the contraction models to get the active tension.
     * 
     *  @param rCalciumConcentrations Reference to a vector of intracellular calcium concentrations at each quadrature point
     *  @param rVoltages Reference to a vector of voltages at each quadrature point
     */
    void SetCalciumAndVoltage(std::vector<double>& rCalciumConcentrations, 
                              std::vector<double>& rVoltages);

    /**
     *  Solve for the deformation, integrating the contraction model ODEs.
     * 
     *  @param time the current time
     *  @param nextTime the next time
     *  @param odeTimestep the ODE timestep
     */
    virtual void Solve(double time, double nextTime, double odeTimestep)=0;
};



template<unsigned DIM>
void AbstractCardiacMechanicsAssembler<DIM>::SetCalciumAndVoltage(std::vector<double>& rCalciumConcentrations, 
                                                                  std::vector<double>& rVoltages)
                                        
{
    assert(rCalciumConcentrations.size() == this->mTotalQuadPoints);
    assert(rVoltages.size() == this->mTotalQuadPoints);

    ContractionModelInputParameters input_parameters;
    
    for(unsigned i=0; i<rCalciumConcentrations.size(); i++)
    {
        input_parameters.intracellularCalciumConcentration = rCalciumConcentrations[i];
        input_parameters.voltage = rVoltages[i];
        
        mContractionModelSystems[i]->SetInputParameters(input_parameters);
    }
}


template<unsigned DIM>
void AbstractCardiacMechanicsAssembler<DIM>::AssembleOnElement(Element<DIM, DIM>& rElement,
                                                               c_matrix<double,STENCIL_SIZE,STENCIL_SIZE>& rAElem,
                                                               c_matrix<double,STENCIL_SIZE,STENCIL_SIZE>& rAElemPrecond,
                                                               c_vector<double,STENCIL_SIZE>& rBElem,
                                                               bool assembleResidual,
                                                               bool assembleJacobian)
{
    // check these have been set
    assert(mCurrentTime != DBL_MAX);
    assert(mNextTime != DBL_MAX);
    assert(mOdeTimestep != DBL_MAX);
    double mech_dt = mNextTime - mCurrentTime;

    // set up the fibre direction for this element
    c_matrix<double,DIM,DIM> r_fibre_sheet_matrix = mpVariableFibreSheetDirections ? (*mpVariableFibreSheetDirections)[rElement.GetIndex()] : mConstantFibreSheetDirections;
    c_vector<double,DIM> fibre_dir;
    for(unsigned i=0; i<DIM; i++)
    {
        fibre_dir(i) = r_fibre_sheet_matrix(i,0);
    }
//r_fibre_sheet_matrix(0,0) = 1.0; 
//r_fibre_sheet_matrix(0,1) = 0;
//r_fibre_sheet_matrix(1,0) = 0;
//r_fibre_sheet_matrix(1,1) = 1.0;


    c_matrix<double, DIM, DIM> jacobian, inverse_jacobian;
    double jacobian_determinant;
    this->mpQuadMesh->GetInverseJacobianForElement(rElement.GetIndex(), jacobian, jacobian_determinant, inverse_jacobian);

    if (assembleJacobian)
    {
        rAElem.clear();
        rAElemPrecond.clear();
    }

    if (assembleResidual)
    {
        rBElem.clear();
    }

    ///////////////////////////////////////////////
    // Get the current displacement at the nodes
    ///////////////////////////////////////////////
    static c_matrix<double,DIM,NUM_NODES_PER_ELEMENT> element_current_displacements;
    static c_vector<double,NUM_VERTICES_PER_ELEMENT> element_current_pressures;
    for (unsigned II=0; II<NUM_NODES_PER_ELEMENT; II++)
    {
        for (unsigned JJ=0; JJ<DIM; JJ++)
        {
            element_current_displacements(JJ,II) = this->mCurrentSolution[DIM*rElement.GetNodeGlobalIndex(II) + JJ];
        }
    }

    ///////////////////////////////////////////////
    // Get the current pressure at the vertices
    ///////////////////////////////////////////////
    for (unsigned II=0; II<NUM_VERTICES_PER_ELEMENT; II++)
    {
        element_current_pressures(II) = this->mCurrentSolution[DIM*this->mpQuadMesh->GetNumNodes() + rElement.GetNodeGlobalIndex(II)];
    }

    // allocate memory for the basis functions values and derivative values
    c_vector<double, NUM_VERTICES_PER_ELEMENT> linear_phi;
    c_vector<double, NUM_NODES_PER_ELEMENT> quad_phi;
    c_matrix<double, DIM, NUM_NODES_PER_ELEMENT> grad_quad_phi;

    // get the material law
    AbstractIncompressibleMaterialLaw<DIM>* p_material_law;
    if (this->mMaterialLaws.size()==1)
    {
        // homogeneous
        p_material_law = this->mMaterialLaws[0];
    }
    else
    {
        // heterogeneous
        #define COVERAGE_IGNORE // not going to have tests in cts for everything
        p_material_law = this->mMaterialLaws[rElement.GetIndex()];
        #undef COVERAGE_IGNORE
    }

    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    //// loop over Gauss points
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    for (unsigned quadrature_index=0; quadrature_index < this->mpQuadratureRule->GetNumQuadPoints(); quadrature_index++)
    {
        unsigned current_quad_point_global_index =   rElement.GetIndex()*this->mpQuadratureRule->GetNumQuadPoints()
                                                   + quadrature_index;

        double wJ = jacobian_determinant * this->mpQuadratureRule->GetWeight(quadrature_index);

        const ChastePoint<DIM>& quadrature_point = this->mpQuadratureRule->rGetQuadPoint(quadrature_index);

        //////////////////////////////////////
        // set up basis function info
        //////////////////////////////////////
        LinearBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, linear_phi);
        QuadraticBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, quad_phi);
        QuadraticBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(quadrature_point, inverse_jacobian, grad_quad_phi);

        ////////////////////////////////////////////////////
        // (dont get the body force)
        ////////////////////////////////////////////////////

        //////////////////////////////////////
        // interpolate grad_u and p
        //////////////////////////////////////
        static c_matrix<double,DIM,DIM> grad_u; // grad_u = (du_i/dX_M)
        grad_u = zero_matrix<double>(DIM,DIM);  // must be on new line!!

        for (unsigned node_index=0; node_index<NUM_NODES_PER_ELEMENT; node_index++)
        {
            for (unsigned i=0; i<DIM; i++)
            {
                for (unsigned M=0; M<DIM; M++)
                {
                    grad_u(i,M) += grad_quad_phi(M,node_index)*element_current_displacements(i,node_index);
                }
            }
        }

        double pressure = 0;
        for (unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
        {
            pressure += linear_phi(vertex_index)*element_current_pressures(vertex_index);
        }


        ///////////////////////////////////////////////
        // calculate C, inv(C) and T
        ///////////////////////////////////////////////
        static c_matrix<double,DIM,DIM> F;
        static c_matrix<double,DIM,DIM> C;
        static c_matrix<double,DIM,DIM> inv_C;
        static c_matrix<double,DIM,DIM> inv_F;
        static c_matrix<double,DIM,DIM> T;

        for (unsigned i=0; i<DIM; i++)
        {
            for (unsigned M=0; M<DIM; M++)
            {
                F(i,M) = (i==M?1:0) + grad_u(i,M);
            }
        }

        C = prod(trans(F),F);
        inv_C = Inverse(C);
        inv_F = Inverse(F);

        double detF = Determinant(F);


        /*****************************
         * The cardiac-specific code 
         *****************************/
        
        // 1. Compute T and dTdE for the PASSIVE part of the strain energy.
        // This is essentially a one-liner (p_material_law->ComputeStressAndStressDerivative()), but the material laws
        // assume the fibre direction is (1,0,0) and sheet direction is (0,1,0), so we have to transform C,inv(C),and T.
        // Let P be the change-of-basis matrix P = (\mathbf{m}_f, \mathbf{m}_s, \mathbf{m}_n). The transformed C for the
        // fibre/sheet basis is C* = P^T C P. We then compute T* = T*(C*), and compute T = P T* P^T.

 
        static c_matrix<double,DIM,DIM> C_transformed;
        static c_matrix<double,DIM,DIM> invC_transformed;
        static c_matrix<double,DIM,DIM> T_transformed;

        C_transformed = prod(trans(r_fibre_sheet_matrix),(c_matrix<double,DIM,DIM>)prod(C,r_fibre_sheet_matrix));          // C*    = P^T C    P
        invC_transformed = prod(trans(r_fibre_sheet_matrix),(c_matrix<double,DIM,DIM>)prod(inv_C,r_fibre_sheet_matrix));   // invC* = P^T invC P

        p_material_law->ComputeStressAndStressDerivative(C_transformed,invC_transformed,pressure,T_transformed,this->dTdE /*transformed*/,assembleJacobian);
        
        T = prod(r_fibre_sheet_matrix, (c_matrix<double,DIM,DIM>)prod(T_transformed, trans(r_fibre_sheet_matrix)));  // T = P T* P^T 

//std::cout << r_fibre_sheet_matrix << C << C_transformed << T_transformed << T <<"\n";


        // compute un-transformed dTdE: dTdE_{MNPQ}  =   P_{Mm}P_{Nn}P_{Pp}P_{Qq} dT*dE*_{mnpq} 
        if(assembleJacobian)
        {
            static FourthOrderTensor<DIM> temp;
            temp.SetAsProduct(this->dTdE, r_fibre_sheet_matrix, 0);
            this->dTdE.SetAsProduct(temp, r_fibre_sheet_matrix, 1);
            temp.SetAsProduct(this->dTdE, r_fibre_sheet_matrix, 2);
            this->dTdE.SetAsProduct(temp, r_fibre_sheet_matrix, 3);
        }

        // 2. Compute the active tension and add to the stress and stress-derivative
        // This is done in the original coordinate system (ie using lambda = m^T C m, where m is the fibre direction, 
        // rather than lam = C(0,0)), so must be done after the passive part has been untransformed.

        double I4 = inner_prod(fibre_dir, prod(C, fibre_dir));
        double lambda = sqrt(I4);

        double active_tension = 0;
        double d_act_tension_dlam = 0.0;     // Set and used if assembleJacobian==true
        double d_act_tension_d_dlamdt = 0.0; // Set and used if assembleJacobian==true

        GetActiveTensionAndTensionDerivs(lambda, current_quad_point_global_index, assembleJacobian,
                                         active_tension, d_act_tension_dlam, d_act_tension_d_dlamdt);
	
        // amend the stress and dTdE using the active tension
        double dTdE_coeff = -2*active_tension/(I4*I4); // note: I4*I4 = lam^4
        if(IsImplicitSolver())
        {     
            dTdE_coeff += (d_act_tension_dlam + d_act_tension_d_dlamdt/mech_dt)/(lambda*I4); // note: I4*lam = lam^3
        }

        T += (active_tension/I4)*outer_prod(fibre_dir,fibre_dir);

        for (unsigned M=0; M<DIM; M++)
        {
            for (unsigned N=0; N<DIM; N++)
            {
                for (unsigned P=0; P<DIM; P++)
                {
                    for (unsigned Q=0; Q<DIM; Q++)
                    {
                        this->dTdE(M,N,P,Q) +=  dTdE_coeff*fibre_dir(M)*fibre_dir(N)*fibre_dir(P)*fibre_dir(Q);
                    }
                }
            }
        }

        /*******************************************************************
         * End of cardiac specific code 
         * 
         * The following, excluding the lack of body force term, should be 
         * identical to NonlinearElasticityAssembler::AssembleOnElement
         *******************************************************************/
         
         
        static FourthOrderTensor<DIM> dTdE_F;
        static FourthOrderTensor<DIM> dTdE_FF1;
        static FourthOrderTensor<DIM> dTdE_FF2;
  
        dTdE_F.SetAsProduct(this->dTdE, F, 1);  // B^{MdPQ}  = F^d_N * dTdE^{MdPQ}
        dTdE_FF1.SetAsProduct(dTdE_F, F, 3);    // B1^{MdPe} = F^d_N * F^e_Q * dTdE^{MNPQ} 
        dTdE_FF2.SetAsProduct(dTdE_F, F, 2);    // B2^{MdeQ} = F^d_N * F^e_P * dTdE^{MNPQ}


        /////////////////////////////////////////
        // residual vector
        /////////////////////////////////////////
        if (assembleResidual)
        {
            for (unsigned index=0; index<NUM_NODES_PER_ELEMENT*DIM; index++)
            {
                unsigned spatial_dim = index%DIM;
                unsigned node_index = (index-spatial_dim)/DIM;

                assert(node_index < NUM_NODES_PER_ELEMENT);

                // no body force bit as body force = 0

                for (unsigned M=0; M<DIM; M++)
                {
                    for (unsigned N=0; N<DIM; N++)
                    {
                        rBElem(index) +=   T(M,N)
                                         * F(spatial_dim,M)
                                         * grad_quad_phi(N,node_index)
                                         * wJ;
                    }
                }
            }

            for (unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
            {
                rBElem( NUM_NODES_PER_ELEMENT*DIM + vertex_index ) +=   linear_phi(vertex_index)
                                                                      * (detF - 1)
                                                                      * wJ;
            }
        }

        /////////////////////////////////////////
        // Jacobian matrix
        /////////////////////////////////////////
        if (assembleJacobian)
        {
            for (unsigned index1=0; index1<NUM_NODES_PER_ELEMENT*DIM; index1++)
            {
                unsigned spatial_dim1 = index1%DIM;
                unsigned node_index1 = (index1-spatial_dim1)/DIM;


                for (unsigned index2=0; index2<NUM_NODES_PER_ELEMENT*DIM; index2++)
                {
                    unsigned spatial_dim2 = index2%DIM;
                    unsigned node_index2 = (index2-spatial_dim2)/DIM;

                    for (unsigned M=0; M<DIM; M++)
                    {
                        for (unsigned N=0; N<DIM; N++)
                        {
                            rAElem(index1,index2) +=   T(M,N)
                                                     * grad_quad_phi(M,node_index1)
                                                     * grad_quad_phi(N,node_index2)
                                                     * (spatial_dim1==spatial_dim2?1:0)
                                                     * wJ;
                        }
                    }
                    
                    for (unsigned M=0; M<DIM; M++)
                    {
                        for (unsigned P=0; P<DIM; P++)
                        {
                            rAElem(index1,index2)  +=   0.5
                                                      * dTdE_FF1(M,spatial_dim1,P,spatial_dim2)
                                                      * grad_quad_phi(P,node_index2)
                                                      * grad_quad_phi(M,node_index1)
                                                      * wJ;
                        }
                        
                        for (unsigned Q=0; Q<DIM; Q++)
                        {
                           rAElem(index1,index2)  +=   0.5
                                                     * dTdE_FF2(M,spatial_dim1,spatial_dim2,Q)
                                                     * grad_quad_phi(Q,node_index2)
                                                     * grad_quad_phi(M,node_index1)
                                                     * wJ;
                        }
                    }
                }

                for (unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
                {
                    unsigned index2 = NUM_NODES_PER_ELEMENT*DIM + vertex_index;

                    for (unsigned M=0; M<DIM; M++)
                    {
                         rAElem(index1,index2)  +=  - inv_F(M,spatial_dim1)
                                                    * grad_quad_phi(M,node_index1)
                                                    * linear_phi(vertex_index)
                                                    * wJ;
                    }
                }
            }

            for (unsigned vertex_index=0; vertex_index<NUM_VERTICES_PER_ELEMENT; vertex_index++)
            {
                unsigned index1 = NUM_NODES_PER_ELEMENT*DIM + vertex_index;

                for (unsigned index2=0; index2<NUM_NODES_PER_ELEMENT*DIM; index2++)
                {
                    unsigned spatial_dim2 = index2%DIM;
                    unsigned node_index2 = (index2-spatial_dim2)/DIM;

                    for (unsigned M=0; M<DIM; M++)
                    {
                        // same as (negative of) the opposite block (ie a few lines up), except for detF
                        rAElem(index1,index2) +=   detF
                                                 * inv_F(M,spatial_dim2)
                                                 * grad_quad_phi(M,node_index2)
                                                 * linear_phi(vertex_index)
                                                 * wJ;
                    }
                }

                /////////////////////////////////////////////////////
                // Preconditioner matrix
                // Fill the mass matrix (ie \intgl phi_i phi_j) in the
                // pressure-pressure block. Note, the rest of the
                // entries are filled in at the end
                /////////////////////////////////////////////////////
                for (unsigned vertex_index2=0; vertex_index2<NUM_VERTICES_PER_ELEMENT; vertex_index2++)
                {
                    unsigned index2 = NUM_NODES_PER_ELEMENT*DIM + vertex_index2;
                    rAElemPrecond(index1,index2) +=   linear_phi(vertex_index)
                                                    * linear_phi(vertex_index2)
                                                    * wJ;
                }
            }
        }
    }


    if (assembleJacobian)
    {
        // Fill in the other blocks of the preconditioner matrix. (This doesn't
        // effect the pressure-pressure block of the rAElemPrecond as the
        // pressure-pressure block of rAElem is zero
        rAElemPrecond = rAElemPrecond + rAElem;
    }
}


template<unsigned DIM>
void AbstractCardiacMechanicsAssembler<DIM>::SetVariableFibreSheetDirections(std::string orthoFile)
{
    if((orthoFile.length()<7) || orthoFile.substr(orthoFile.length()-6,orthoFile.length()) != ".ortho")
    {
        EXCEPTION("Fibre file must be a .ortho file");
    }

    std::ifstream ifs(orthoFile.c_str());
    if (!ifs.is_open())
    {
        EXCEPTION("Could not open file: " + orthoFile);
    }
    
    unsigned num_elem_read_from_file;
    ifs >> num_elem_read_from_file;
    assert(num_elem_read_from_file == this->mpQuadMesh->GetNumElements());
    
    mpVariableFibreSheetDirections = new std::vector<c_matrix<double,DIM,DIM> >(this->mpQuadMesh->GetNumElements(), zero_matrix<double>(DIM,DIM));
    for(unsigned elem_index=0; elem_index<this->mpQuadMesh->GetNumElements(); elem_index++)
    {
        for(unsigned j=0; j<DIM*DIM; j++) 
        {
            double data;
            ifs >> data; 
            if(ifs.fail())
            {
                std::stringstream error_message;
                error_message << "Error occurred when reading file " << orthoFile
                              << ". Expected " << this->mpQuadMesh->GetNumElements() << " rows and "
                              << "three (not DIM!) columns, ie three components for each fibre, whichever "
                              << "dimension you are in";
                delete mpVariableFibreSheetDirections;           
                EXCEPTION(error_message.str());
            }
            
            (*mpVariableFibreSheetDirections)[elem_index](j/DIM,j%DIM) = data;
        }

//EMTODO:
//CheckOrthogonality((*mpVariableFibreSheetDirections)[elem_index]);
    }
   
    ifs.close();
}


#endif /*ABSTRACTCARDIACMECHANICSASSEMBLER_HPP_*/
