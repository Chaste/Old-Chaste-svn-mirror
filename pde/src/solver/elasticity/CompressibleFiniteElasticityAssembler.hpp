#ifndef COMPRESSIBLEFINITEELASTICITYASSEMBLER_
#define COMPRESSIBLEFINITEELASTICITYASSEMBLER_

#include "AbstractNonlinearStaticAssembler.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "ElasticityBoundaryConditionsContainer.hpp"
#include "FourthOrderTensor.hpp"


template<unsigned DIM>
class CompressibleFiniteElasticityAssembler : public AbstractNonlinearStaticAssembler<DIM,DIM,DIM>
{
private:
    double mDensity;
    c_vector<double, DIM> mGravity; 
    
    virtual c_matrix<double,DIM*(DIM+1),DIM*(DIM+1)> ComputeMatrixTerm(
        c_vector<double, DIM+1> &rPhi,
        c_matrix<double, DIM, DIM+1> &rGradPhi,
        Point<DIM> &X,                           // old position
        c_vector<double,DIM> &u,                 // displacement
        c_matrix<double,DIM,DIM> &rGradU)        // du/dX (displacement gradient) 
    {
        c_matrix<double,DIM*(DIM+1),DIM*(DIM+1)> ret;
        ret.clear();
 
        double c1 = 2;
        double c2 = 1;
        double c3 = -4;
        
        c_matrix<double,DIM,DIM> F = rGradU+identity_matrix<double>(DIM,DIM);
        c_matrix<double,DIM,DIM> C = prod(trans(F),F);
        c_matrix<double,DIM,DIM> invC = Inverse(C);

        double I_1 = Trace(C);
//        double I_2 = SecondInvariant(C);
        double I_3 = Determinant(C);
        
        c_matrix<double,DIM,DIM> T =   2*c1*identity_matrix<double>(DIM,DIM)
                                     + 2*c2*(I_1*identity_matrix<double>(DIM,DIM) - C)
                                     + 2*c3*I_3*Inverse(C);
                                     
        FourthOrderTensor<DIM> dTdE;   // BAD!!!!!!!!!!!!!!!!!!!!!!!! change 4thOrderTensor to not use std::vectors, or use some better written code
        c_matrix<double,DIM,DIM> eye = identity_matrix<double>(DIM);
        
        for(unsigned M=0;M<DIM;M++)
        {
            for(unsigned N=0;N<DIM;N++)
            {
                for(unsigned P=0;P<DIM;P++)
                {
                    for(unsigned Q=0;Q<DIM;Q++)
                    {
                        dTdE.rGetVal()[M][N][P][Q] =   4*c2*eye(M,N)*eye(P,Q)
                                                     - 4*c2*eye(M,P)*eye(N,Q)
                                                     + 4*c3*I_3*invC(M,N)*invC(P,Q)
                                                     - 4*c3*I_3*invC(M,P)*invC(Q,N);
                    }
                }
            }
        }
 
 
 
        for (unsigned I=0; I<DIM+1; I++) // I=node index
        {
            for (unsigned s=0; s<DIM; s++) // s=spatial dimension index
            {
                for (unsigned J=0; J<DIM+1; J++) //J=node index
                {
                    for (unsigned t=0; t<DIM; t++) //t=spatial dimension index
                    {
                        for(unsigned M=0; M<DIM; M++)
                        {
                            for(unsigned N=0; N<DIM; N++)
                            {
                                ret(DIM*I+s, DIM*J+t) += T(M,N)*rGradPhi(M,J)*rGradPhi(N,I)*eye(s,t);
                                
                                for(unsigned P=0; P<DIM; P++)
                                {
                                    for(unsigned Q=0; Q<DIM; Q++)
                                    {
                                        ret(DIM*I+s, DIM*J+t) +=  dTdE.rGetVal()[M][N][P][Q]  
                                                                 * 0.5*(rGradPhi(Q,J)*F(t,P) + rGradPhi(P,J)*F(t,Q))
                                                                 * F(s,M) * rGradPhi(N,I);
                                    }
                                }                                        
                            }
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
        Point<DIM> &X,
        c_vector<double,DIM> &u,
        c_matrix<double,DIM,DIM> &rGradU)
    {
        c_vector<double,DIM*(DIM+1)> ret;
        ret.clear();
        
        double c1 = 2;
        double c2 = 1;
        double c3 = -4;
        
        c_matrix<double,DIM,DIM> F = rGradU+identity_matrix<double>(DIM,DIM);
        c_matrix<double,DIM,DIM> C = prod(trans(F),F);


        double I_1 = Trace(C);
//        double I_2 = SecondInvariant(C);
        double I_3 = Determinant(C);
        
        c_matrix<double,DIM,DIM> T =   2*c1*identity_matrix<double>(DIM,DIM)
                                     + 2*c2*(I_1*identity_matrix<double>(DIM,DIM) - C)
                                     + 2*c3*I_3*Inverse(C);
        
        
        c_matrix<double,DIM,DIM> temp = prod(F,T);  // F_iM * T_MN = temp_iN
        
        for (unsigned I=0; I<DIM+1; I++) // I = node_index
        {
            for (unsigned s=0; s<DIM; s++) // s = spatial dimension index
            {
                ret(DIM*I+s) = mDensity * mGravity(s) * rPhi(I);
                for(unsigned M=0; M<DIM; M++)
                {
                    for(unsigned N=0; N<DIM; N++)
                    {
                        ret(DIM*I+s) += T(M,N)*F(s,M)*rGradPhi(N,I);   //temp(s,N)*rGradPhi(N,I);
                    }
                }
            }
        }
        return ret;
    }
    
    
    virtual c_vector<double, DIM*DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<DIM-1,DIM> &rSurfaceElement,
        c_vector<double, DIM> &rPhi,
        Point<DIM> &X )
    {
        assert(0);
        c_vector<double,DIM*DIM> ret;
        for (unsigned I=0; I<DIM; I++) // DIM = number of nodes, in this context (as element is a surface element)
        {
            for (unsigned s=0; s<DIM; s++) // s = spatial dimension index
            {
                // GetNeumannBCValue(&rSurfaceElement,x,s) = s-component of traction
                ret(DIM*I+s) = rPhi(I) * this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement,X,s);
            }
        }
        return ret;
    }
    
   
    
public :
    CompressibleFiniteElasticityAssembler(ConformingTetrahedralMesh<DIM,DIM>* pMesh,
                                          ElasticityBoundaryConditionsContainer<DIM>* pBoundaryConditions)
             :  AbstractNonlinearStaticAssembler<DIM,DIM,DIM>()
    {
        assert(pMesh != NULL);
        assert(pBoundaryConditions != NULL);
        
        this->mpMesh = pMesh;
        this->mpBoundaryConditions = pBoundaryConditions;
        
        mGravity = zero_vector<double>(DIM);
        mDensity = 0;
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
    
};

#endif /*COMPRESSIBLEFINITEELASTICITYASSEMBLER_*/
