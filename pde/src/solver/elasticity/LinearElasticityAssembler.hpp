#ifndef LINEARELASTICITYASSEMBLER_HPP_
#define LINEARELASTICITYASSEMBLER_HPP_

#include "AbstractLinearStaticProblemAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"

template <int DIM>
class LinearElasticityAssembler : public AbstractLinearStaticProblemAssembler<DIM,DIM,DIM>
{
private :
    double mYoungsModulus;            // units??!!!!!!!!!
    double mPoissonsRatio;            // dimensionless
    c_vector<double, DIM>  mGravity;  // units??
    
    virtual c_matrix<double,DIM*(DIM+1),DIM*(DIM+1)> ComputeLhsTerm(
        const c_vector<double, DIM+1> &rPhi,
        const c_matrix<double, DIM, DIM+1> &rGradPhi,
        const Point<DIM> &rX,
        const c_vector<double,DIM> &u)
    {
        return zero_matrix<double>(DIM*(DIM+1)); // todo: fill in
    }
    
    
    virtual c_vector<double,DIM*(DIM+1)> ComputeRhsTerm(
        const c_vector<double, DIM+1> &rPhi,
        const Point<DIM> &rX,
        const c_vector<double,DIM> &u)
    {
        return zero_vector<double>(DIM*(DIM+1)); // todo: fill in
    }  

    virtual c_vector<double, DIM*DIM> ComputeSurfaceRhsTerm(
        const BoundaryElement<DIM-1,DIM> &rSurfaceElement,
        const c_vector<double, DIM> &phi,
        const Point<DIM> &x )
    {
        assert(0);
        return zero_vector<double>(DIM*(DIM+1)); 
    }
    
public :
    LinearElasticityAssembler(ConformingTetrahedralMesh<DIM,DIM>* pMesh,
                              BoundaryConditionsContainer<DIM,DIM,DIM>* pBoundaryConditions)
    {
        assert(pMesh!=NULL);
        assert(pBoundaryConditions!=NULL);
        
        this->mpMesh = pMesh;
        this->mpBoundaryConditions = pBoundaryConditions;
    }
    
    void SetParameters(double youngsModulus, double poissonsRatio, c_vector<double, DIM> gravity)
    {
        if(youngsModulus < 0.0)
        {
            EXCEPTION("Young's Modulus, E, should be positive");
        }
        if(poissonsRatio > 0.5)
        {
            EXCEPTION("Poisson's ratio, nu, must be less than 0.5");
        }
        
        mYoungsModulus = youngsModulus;
        mPoissonsRatio = poissonsRatio;
        
        mGravity = gravity;
   }
};

#endif /*LINEARELASTICITYASSEMBLER_HPP_*/
