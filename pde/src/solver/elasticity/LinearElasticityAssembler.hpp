#ifndef LINEARELASTICITYASSEMBLER_HPP_
#define LINEARELASTICITYASSEMBLER_HPP_

#include "AbstractLinearStaticProblemAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "ElasticityBoundaryConditionsContainer.hpp"

#define LAME_COEFF_UNSET -1e10



/**
 *  LinearElasticityAssembler
 *
 *  Solve the isotropic linear elasticity equations
 *
 *      sigma_{ij,j} + rho g_i = 0,
 *
 *  where
 *
 *      x_i  = old position,
 *      rho  = density,
 *      g_i  = body force (per unit volume) (eg gravity)
 *
 *  and the stress sigma_ij is given by
 *
 *      sigma_ij = mu (u_{i,j} + u_{j,i}) +  lambda u_{k,k} delta_{ij}
 *               =  2 mu e_{i,j} + lamda u_{k,k} delta_{ij}
 *
 *  where u_i is the displacement, and lambda and mu are the Lame coefficients (and
 *  e_{i,j} = 0.5*(u_{i,j}+u_{j,i}) is the infinitessimal strain).
 *
 *  NOTE: currently only solves heterogeneous problems.
 */
template <unsigned DIM>
class LinearElasticityAssembler : public AbstractLinearStaticProblemAssembler<DIM,DIM,DIM>
{
private :
    double mLambda;                  // units: anything, but {lambda, mu, rho, g} must be consistent
    double mMu;                      // units: anything, but {lambda, mu, rho, g} must be consistent
    double mDensity;                 // units: anything, but {lambda, mu, rho, g} must be consistent
    c_vector<double, DIM> mGravity;  // units: anything, but {lambda, mu, rho, g} must be consistent
    
    bool mLameCoefficientsSet;  // bool saying whether the lame coeffs have been set yet
    
    virtual c_matrix<double,DIM*(DIM+1),DIM*(DIM+1)> ComputeMatrixTerm(
        c_vector<double, DIM+1> &rPhi,
        c_matrix<double, DIM, DIM+1> &rGradPhi,
        Point<DIM> &rX,
        c_vector<double,DIM> &u,
        c_matrix<double,DIM,DIM> &rGradU)
    {
        c_matrix<double,DIM*(DIM+1),DIM*(DIM+1)> ret;
        
        for (unsigned I=0; I<DIM+1; I++) // I=node index
        {
            for (unsigned s=0; s<DIM; s++) // s=spatial dimension index
            {
                for (unsigned J=0; J<DIM+1; J++) //J=node index
                {
                    for (unsigned t=0; t<DIM; t++) //t=spatial dimension index
                    {
                        ret(DIM*I+s, DIM*J+t) = mLambda*rGradPhi(s,I)*rGradPhi(t,J) + mMu*rGradPhi(t,I)*rGradPhi(s,J);
                        
                        for (unsigned j=0; j<DIM; j++)
                        {
                            ret(DIM*I+s, DIM*J+t) += mMu*(s==t)*rGradPhi(j,I)*rGradPhi(j,J);
                        }
                    }
                }
            }
        }
        return ret;
    }
    
    
    virtual c_vector<double,DIM*(DIM+1)> ComputeVectorTerm(
        c_vector<double, DIM+1> &rPhi,
        c_matrix<double, DIM, DIM+1> &rGradPhi,
        Point<DIM> &rX,
        c_vector<double,DIM> &u,
        c_matrix<double,DIM,DIM> &rGradU)
    {
        c_vector<double,DIM*(DIM+1)> ret;
        for (unsigned I=0; I<DIM+1; I++) // I = node_index
        {
            for (unsigned s=0; s<DIM; s++) // s = spatial dimension index
            {
                ret(DIM*I+s) = mDensity * mGravity(s) * rPhi(I);
            }
        }
        return ret;
    }
    
    virtual c_vector<double, DIM*DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<DIM-1,DIM> &rSurfaceElement,
        c_vector<double, DIM> &rPhi,
        Point<DIM> &rX )
    {
        c_vector<double,DIM*DIM> ret;
        for (unsigned I=0; I<DIM; I++) // DIM = number of nodes, in this context (as element is a surface element)
        {
            for (unsigned s=0; s<DIM; s++) // s = spatial dimension index
            {
                // GetNeumannBCValue(&rSurfaceElement,x,s) = s-component of traction
                ret(DIM*I+s) = rPhi(I) * this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement,rX,s);
            }
        }
        return ret;
    }
    
    
    /**
     *  This method is called at the beginning of Solve() and used here to check
     *  we are ready to go
     */
    void PrepareForSolve()
    {
        AbstractLinearStaticProblemAssembler<DIM,DIM,DIM>::PrepareForSolve();
        
        if (mLameCoefficientsSet == false)
        {
            EXCEPTION("Lame coefficients have not been set yet");
        }
    }
    
public :
    LinearElasticityAssembler(ConformingTetrahedralMesh<DIM,DIM>* pMesh,
                              ElasticityBoundaryConditionsContainer<DIM>* pBoundaryConditions)
            :  AbstractLinearStaticProblemAssembler<DIM,DIM,DIM>()
    {
        assert(pMesh != NULL);
        assert(pBoundaryConditions != NULL);
        
        this->mpMesh = pMesh;
        this->mpBoundaryConditions = pBoundaryConditions;
        
        mGravity = zero_vector<double>(DIM);
        mDensity = 0;
        
        mLameCoefficientsSet = false;
    }
    
    /**
     *  SetLameCoefficients
     * 
     *  One of three methods allowing the Lame coeffients (lambda and mu) to be set. 
     *  Alternatively, SetYoungsModulusAndPoissonsRatio() (E and nu) or 
     *  SetBulkModulusAndShearModulus() (K and G) can be called.
     */
    void SetLameCoefficients(double lambda, double mu)
    {
        mLambda = lambda;
        mMu = mu;
        
        mLameCoefficientsSet = true;
    }
    
    /**
     *  SetYoungsModulusAndPoissonsRatio
     * 
     *  Set the Young's modulus (usually denoted E) and Poisson's ratio (usually 
     *  denoted nu). The equivalent Lame coefficients will automatically be 
     *  computed.
     * 
     *  One of three methods allowing the Lame coeffients (lambda and mu) to be set. 
     *  Alternatively, SetLameCoefficients() or SetBulkModulusAndShearModulus() 
     *  (K and G) can be called.
     */
    void SetYoungsModulusAndPoissonsRatio(double youngsModulus, double poissonsRatio)
    {
        if (youngsModulus < 0.0)
        {
            EXCEPTION("Young's Modulus, E, should be positive");
        }
        if (poissonsRatio > 0.5)
        {
            EXCEPTION("Poisson's ratio, nu, must be less than 0.5");
        }
        
        mLambda = youngsModulus * poissonsRatio / ( (1+poissonsRatio)*(1-2*poissonsRatio));
        mMu = youngsModulus / (2*(1+poissonsRatio));
        
        mLameCoefficientsSet = true;
    }
    
    
    /**
     *  SetBulkModulusAndShearModulus
     * 
     *  Set the Bulk modulus (usually denoted K) and Shear Modulus (usually 
     *  denoted G). The equivalent Lame coefficients will automatically be 
     *  computed.
     * 
     *  One of three methods allowing the Lame coeffients (lambda and mu) to be set. 
     *  Alternatively, SetLameCoefficients() or SetYoungsModulusAndPoissonsRatio() 
     *  (E and nu) can be called.
     */
    void SetBulkModulusAndShearModulus(double bulkModulus, double shearModulus)
    {
        mLambda = bulkModulus - (2.0*shearModulus/3);
        mMu = shearModulus;
        
        mLameCoefficientsSet = true;
    }
    
    
    /**
     *  SetDensityAndGravity
     * 
     *  Set the total body force by setting the body density and 
     */
    void SetDensityAndGravity(double density, c_vector<double, DIM> gravity)
    {
        if (density < 0.0)
        {
            EXCEPTION("Density must be positive");
        }
        mDensity = density;
        mGravity = gravity;
    }
    
    
    /**
     *  Returns the Lame coefficient lambda. Use after calling 
     *  SetYoungsModulusAndPoissonsRatio() say 
     */
    double GetLambda()
    {
        return mLambda;
    }
    
    /**
     *  Returns the Lame coefficient mu. Use after calling 
     *  SetYoungsModulusAndPoissonsRatio() say
     */
    double GetMu()
    {
        return mMu;
    }
    
};

#endif /*LINEARELASTICITYASSEMBLER_HPP_*/
