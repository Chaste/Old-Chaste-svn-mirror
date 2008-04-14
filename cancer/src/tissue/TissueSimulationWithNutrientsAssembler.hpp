#ifndef _TISSUESIMULATIONWITHNUTRIENTSASSEMBLER_HPP_
#define _TISSUESIMULATIONWITHNUTRIENTSASSEMBLER_HPP_

#include <vector>
#include <petscvec.h>

#include "ConformingTetrahedralMesh.hpp"
#include "SimpleLinearEllipticAssembler.hpp"
#include "GaussianQuadratureRule.hpp"


/**
 *  TissueSimulationWithNutrientsAssembler
 * 
 *  This is a purpose made elliptic assembler that interpolates the source terms 
 *  from node onto gauss points, as for a nutrients simulation the source will only
 *  be known at the cells (nodes), not the gauss points.
 */
template<unsigned DIM>
class TissueSimulationWithNutrientsAssembler
    : public SimpleLinearEllipticAssembler<DIM, DIM, TissueSimulationWithNutrientsAssembler<DIM> >
{
    typedef SimpleLinearEllipticAssembler<DIM, DIM, TissueSimulationWithNutrientsAssembler<DIM> > BaseClassType;
    friend class AbstractStaticAssembler<DIM, DIM, 1u, true, BaseClassType>;
    
private:

    /** The constant in u part of the source term, interpolated onto
     *  the current point 
     */
    double mConstantInUSourceTerm;
    
    /** The linear in u part of the source term, interpolated onto
     *  the current point 
     */
    double mLinearInUCoeffInSourceTerm;
    
protected:

    /**
     *  The SimpleLinearEllipticAssembler version of this method is 
     *  overloaded using the interpolated source term
     */
    virtual c_vector<double,1*(DIM+1)> ComputeVectorTerm(
        c_vector<double, DIM+1> &rPhi,
        c_matrix<double, DIM, DIM+1> &rGradPhi,
        ChastePoint<DIM> &rX,
        c_vector<double,1>& rU,
        c_matrix<double, 1, DIM> &rGradU /* not used */,
        Element<DIM,DIM>* pElement);

    /**
     *  The SimpleLinearEllipticAssembler version of this method is 
     *  overloaded using the interpolated source term
     */
    virtual c_matrix<double,1*(DIM+1),1*(DIM+1)> ComputeMatrixTerm(
        c_vector<double, DIM+1> &rPhi,
        c_matrix<double, DIM, DIM+1> &rGradPhi,
        ChastePoint<DIM> &rX,
        c_vector<double,1> &u,
        c_matrix<double,1,DIM> &rGradU,
        Element<DIM,DIM>* pElement);
    
    void ResetInterpolatedQuantities();    
    
    void IncrementInterpolatedQuantities(double phiI, const Node<DIM> *pNode);    
    
public:

    /**
     * Constructor stores the mesh and pde and boundary conditions.
     */
    TissueSimulationWithNutrientsAssembler(ConformingTetrahedralMesh<DIM,DIM>* pMesh,
                                  AbstractLinearEllipticPde<DIM>* pPde,
                                  BoundaryConditionsContainer<DIM,DIM,1>* pBoundaryConditions,
                                  unsigned numQuadPoints = 2);
    
    /**
     *  Destructor
     */
    ~TissueSimulationWithNutrientsAssembler();
    
};

template<unsigned DIM>
TissueSimulationWithNutrientsAssembler<DIM>::TissueSimulationWithNutrientsAssembler(ConformingTetrahedralMesh<DIM,DIM>* pMesh,
                              AbstractLinearEllipticPde<DIM>* pPde,
                              BoundaryConditionsContainer<DIM,DIM,1>* pBoundaryConditions,
                              unsigned numQuadPoints) :
        BaseClassType(pMesh, pPde, pBoundaryConditions, numQuadPoints)
{
}

template<unsigned DIM>
TissueSimulationWithNutrientsAssembler<DIM>::~TissueSimulationWithNutrientsAssembler()
{
}

template<unsigned DIM>
c_vector<double,1*(DIM+1)> TissueSimulationWithNutrientsAssembler<DIM>::ComputeVectorTerm(
        c_vector<double, DIM+1> &rPhi,
        c_matrix<double, DIM, DIM+1> &rGradPhi,
        ChastePoint<DIM> &rX,
        c_vector<double,1>& rU,
        c_matrix<double, 1, DIM> &rGradU /* not used */,
        Element<DIM,DIM>* pElement)
{
    return mConstantInUSourceTerm * rPhi;
}

template<unsigned DIM>
c_matrix<double,1*(DIM+1),1*(DIM+1)> TissueSimulationWithNutrientsAssembler<DIM>::ComputeMatrixTerm(
        c_vector<double, DIM+1> &rPhi,
        c_matrix<double, DIM, DIM+1> &rGradPhi,
        ChastePoint<DIM> &rX,
        c_vector<double,1> &u,
        c_matrix<double,1,DIM> &rGradU,
        Element<DIM,DIM>* pElement)
{
    c_matrix<double, DIM, DIM> pde_diffusion_term = this->mpEllipticPde->ComputeDiffusionTerm(rX);
    
    // if statement just saves computing phi*phi^T if it is to be multiplied by zero
    if(mLinearInUCoeffInSourceTerm!=0)
    {
        return   prod( trans(rGradPhi), c_matrix<double, DIM, DIM+1>(prod(pde_diffusion_term, rGradPhi)) )
               - mLinearInUCoeffInSourceTerm * outer_prod(rPhi,rPhi);
    }
    else
    {
        return   prod( trans(rGradPhi), c_matrix<double, DIM, DIM+1>(prod(pde_diffusion_term, rGradPhi)) );
    }
}

template<unsigned DIM>
void TissueSimulationWithNutrientsAssembler<DIM>::ResetInterpolatedQuantities()
{
    mConstantInUSourceTerm = 0;
    mLinearInUCoeffInSourceTerm = 0;
}

template<unsigned DIM>
void TissueSimulationWithNutrientsAssembler<DIM>::IncrementInterpolatedQuantities(double phiI, const Node<DIM> *pNode)
{
    mConstantInUSourceTerm += phiI * this->mpEllipticPde->ComputeConstantInUSourceTermAtNode(*pNode);
    mLinearInUCoeffInSourceTerm += phiI * this->mpEllipticPde->ComputeLinearInUCoeffInSourceTermAtNode(*pNode);
}

template<unsigned DIM>
struct AssemblerTraits<TissueSimulationWithNutrientsAssembler<DIM> >
{
    typedef TissueSimulationWithNutrientsAssembler<DIM> CVT_CLS;
    typedef TissueSimulationWithNutrientsAssembler<DIM> CMT_CLS;
    typedef TissueSimulationWithNutrientsAssembler<DIM> INTERPOLATE_CLS;
};

#endif //_TISSUESIMULATIONWITHNUTRIENTSASSEMBLER_HPP_
