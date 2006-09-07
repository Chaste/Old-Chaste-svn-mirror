#ifndef LINEARELASTICITYASSEMBLER_HPP_
#define LINEARELASTICITYASSEMBLER_HPP_

#include "AbstractLinearStaticProblemAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "ElasticityBoundaryConditionsContainer.hpp"

template <int DIM>
class LinearElasticityAssembler : public AbstractLinearStaticProblemAssembler<DIM,DIM,DIM>
{
private :
    double mLambda;                  // units??!!!!!!!!!
    double mMu;                      // units??
    double mDensity;                 // units??
    c_vector<double, DIM> mGravity;  // units??
    
    
    
    virtual c_matrix<double,DIM*(DIM+1),DIM*(DIM+1)> ComputeLhsTerm(
        const c_vector<double, DIM+1> &rPhi,
        const c_matrix<double, DIM, DIM+1> &rGradPhi,
        const Point<DIM> &rX,
        const c_vector<double,DIM> &u)
    {
        c_matrix<double,DIM*(DIM+1),DIM*(DIM+1)> ret;
        
        for(unsigned I=0; I<DIM+1; I++) // I=node index
        {
            for(unsigned s=0; s<DIM; s++) // s=spatial dimension index
            {
                for(unsigned J=0; J<DIM+1; J++) //J=node index
                {
                    for(unsigned t=0; t<DIM; t++) //t=spatial dimension index
                    {
                        ret(DIM*I+s, DIM*J+t) = mLambda*rGradPhi(s,I)*rGradPhi(t,J) + mMu*rGradPhi(t,I)*rGradPhi(s,J);
                        
                        for(unsigned j=0; j<DIM; j++)
                        {
                            ret(DIM*I+s, DIM*J+t) += mMu*(s==t)*rGradPhi(j,I)*rGradPhi(j,J);
                        }
                    }
                }
            }
        }
        return ret;
    }
    
    
    virtual c_vector<double,DIM*(DIM+1)> ComputeRhsTerm(
        const c_vector<double, DIM+1> &rPhi,
        const Point<DIM> &rX,
        const c_vector<double,DIM> &u)
    {
        c_vector<double,DIM*(DIM+1)> ret;
        for(unsigned I=0; I<DIM+1; I++) // I = node_index
        {
            for(unsigned s=0; s<DIM; s++) // s = spatial dimension index
            {
                ret(DIM * I + s) = mDensity * mGravity(s) * rPhi(I);
            }
        }
        return ret;
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
                              ElasticityBoundaryConditionsContainer<DIM>* pBoundaryConditions)
    {
        assert(pMesh!=NULL);
        assert(pBoundaryConditions!=NULL);
        
        this->mpMesh = pMesh;
        this->mpBoundaryConditions = pBoundaryConditions;
        
        mGravity = zero_vector<double>(DIM);
    }
    
    
    ///\todo: finish this 
    /*
    void SetYoungsModulusAndPoissonsRatio(double youngsModulus, double poissonsRatio)
    {
        if(youngsModulus < 0.0)
        {
            EXCEPTION("Young's Modulus, E, should be positive");
        }
        if(poissonsRatio > 0.5)
        {
            EXCEPTION("Poisson's ratio, nu, must be less than 0.5");
        }
        
        mLambda = .....;
        mMu = .....;
    }
    */
    
    void SetLameCoefficients(double lambda, double mu)
    {
        mLambda = lambda;
        mMu = mu;
    }
    
    
    void SetDensityAndGravity(double density, c_vector<double, DIM> gravity)
    {
        if(density < 0.0)
        {
            EXCEPTION("Density must be positive");
        }
        mDensity = density;
        mGravity = gravity;
    }
};

#endif /*LINEARELASTICITYASSEMBLER_HPP_*/
