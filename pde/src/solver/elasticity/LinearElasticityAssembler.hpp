#ifndef LINEARELASTICITYASSEMBLER_HPP_
#define LINEARELASTICITYASSEMBLER_HPP_

#include "AbstractLinearStaticProblemAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "ElasticityBoundaryConditionsContainer.hpp"

#define LAME_COEFF_UNSET -1e10

template <int DIM>
class LinearElasticityAssembler : public AbstractLinearStaticProblemAssembler<DIM,DIM,DIM>
{
private :
    double mLambda;                  // units: anything, but must be {lambda, mu, rho, g} must be consistent
    double mMu;                      // units: anything, but must be {lambda, mu, rho, g} must be consistent
    double mDensity;                 // units: anything, but must be {lambda, mu, rho, g} must be consistent
    c_vector<double, DIM> mGravity;  // units: anything, but must be {lambda, mu, rho, g} must be consistent 
   
    bool mLameCoefficientsSet;  // bool saying whether the lame coeffs have been set yet
    
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
                ret(DIM*I+s) = mDensity * mGravity(s) * rPhi(I);
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


    /**
     *  This method is called at the beginning of Solve() and used here to check
     *  we are ready to go
     */    
    void PrepareForSolve()
    {
        if(mLameCoefficientsSet == false)
        {
            EXCEPTION("Lame coefficients have not been set yet");
        }
    }
    
public :
    LinearElasticityAssembler(ConformingTetrahedralMesh<DIM,DIM>* pMesh,
                              ElasticityBoundaryConditionsContainer<DIM>* pBoundaryConditions)
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
        if(youngsModulus < 0.0)
        {
            EXCEPTION("Young's Modulus, E, should be positive");
        }
        if(poissonsRatio > 0.5)
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
        if(density < 0.0)
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
