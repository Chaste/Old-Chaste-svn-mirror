#ifndef _ABSTRACTMATERIAL_HPP_
#define _ABSTRACTMATERIAL_HPP_

#include "MatrixDouble.hpp"
#include <vector>

//typedef std::vector< std::vector < std::vector< std::vector<double> > > > FourthOrderTensor;
#include "FourthOrderTensor.hpp"

/** AbstractMaterial
 * 
 *  The user should write a derived class which specfies an isotropic strain energy by specifying
 *  the derivatives and second derivatives of W.
 *  
 *  NOTE: make sure the boolean mIsIsotropicLaw is set appropiately in the constructor of any 
 *  derived class
 * 
 * 	This abstract class can compute the stress T and stress derivative dTdE, and Get/Set densities
 * 
 *  \todo Implement the non-isotropic versions of this
 *  \todo Move ComputeStress and Compute_dTdE etc into a new 'helper' class
 *  \todo MakeFourthOrderTensor a class with overloaded multiplication operator
 *  \todo Test the results of Compute_dTdE properly!
 */
 
template <int SPACE_DIM>
class AbstractMaterial 
{

// functions which need to be implemented by derived classes	
private:

	// if isotropic overload these
    virtual double GetdW_by_dI1    (double I1,double I2=0,double I3=0)=0;
    virtual double GetdW_by_dI2    (double I1,double I2=0,double I3=0)=0;
    virtual double GetdW_by_dI3    (double I1,double I2=0,double I3=0)=0;

    virtual double Getd2W_by_dI1   (double I1,double I2=0,double I3=0)=0;
    virtual double Getd2W_by_dI2   (double I1,double I2=0,double I3=0)=0;
    virtual double Getd2W_by_dI3   (double I1,double I2=0,double I3=0)=0;
    virtual double Getd2W_by_dI1dI2(double I1,double I2=0,double I3=0)=0;
    virtual double Getd2W_by_dI2dI3(double I1,double I2=0,double I3=0)=0;
    virtual double Getd2W_by_dI1dI3(double I1,double I2=0,double I3=0)=0;

	// if anisotropic overload this 
    //virtual double GetdW_by_dE(MatrixDouble E, int index1, int index2)=0;
    
    
    
    
public:
    /**
     * Virtual destructor
     */
    virtual ~AbstractMaterial() {};
    
    
	/** Get the density of this material
	 */
    double GetDensity()
    {
		if(!mDensitySet)
		{ 
			std::cout << std::endl << "Error: density not set" << std::endl;
			assert(0);
		}
		return mDensity;
    }
    
	/** Set the density of this material
	 */
    void SetDensity(double density)
	{
		mDensity = density;
		mDensitySet = true;
	}
    
    
    /** Compute this stress in this material at a given value of deformation gradient
     * 	NOTE: The input parameter is deformation gradient F, not strain E
     */
    MatrixDouble ComputeStress(MatrixDouble F)
    {
    	if(mIsIsotropicLaw)
		{
			return ComputeIsotropicStress(F);
		}
		else
		{
			assert(0);
			return F; // Avoid compiler warning
		}
	}
	

    /** Compute this stress derivative in this material at a given value of deformation gradient
     *  The stress derivative is the fourth order tensor dT^{MN}/dE_{PQ} = dW/(dE^{MN}dE^{PQ})
     * 	NOTE: The input parameter is deformation gradient F, not strain E
     */
  	FourthOrderTensor<SPACE_DIM> Compute_dTdE (MatrixDouble F)
    {
    	if(mIsIsotropicLaw)
		{
			return ComputeIsotropic_dTdE(F);
		}
		else
		{
			assert(0);
		}
    }



protected:
    bool   mDensitySet;
    double mDensity;
    bool   mIsIsotropicLaw;    

private:

	MatrixDouble ComputeIsotropicStress(MatrixDouble F)
	{
		assert(F.Rows()==SPACE_DIM);
		
		MatrixDouble C = (F.Transpose())*F;
		switch(SPACE_DIM)
		{
			case 1:
			{	
				double I1 = C.GetFirstInvariant();
				double dW_by_dI1 = GetdW_by_dI1(I1);		
				
				MatrixDouble identity = MatrixDouble::Identity(1);		
				MatrixDouble T(1,1);		
					
				// using dI1/dE = 2*identity 
				T = 2*dW_by_dI1*identity;
				return T;	
			}	
			break;	
			
			case 2:
			{
				double I1 = C.GetFirstInvariant();
				double I2 = C.GetSecondInvariant();
								
				double dW_by_dI1 = GetdW_by_dI1(I1,I2);		
				double dW_by_dI2 = GetdW_by_dI2(I1,I2);		
				
				MatrixDouble identity = MatrixDouble::Identity(2);		
				MatrixDouble T(2,2);		
					
				// using dI1/dE = 2*identity and dI2/dE = 2*I2*inv(C)
				T = (2 * dW_by_dI1 * identity)   +   (2 * dW_by_dI2 * I2 * C.Inverse());
				return T;
			}
			break;

			case 3:
			{		
				double I1 = C.GetFirstInvariant();
				double I2 = C.GetSecondInvariant();
				double I3 = C.GetThirdInvariant();
								
				double dW_by_dI1 = GetdW_by_dI1(I1,I2,I3);		
				double dW_by_dI2 = GetdW_by_dI2(I1,I2,I3);		
				double dW_by_dI3 = GetdW_by_dI3(I1,I2,I3);		
				
				MatrixDouble identity = MatrixDouble::Identity(3);		
				MatrixDouble T(3,3);		
					
				// using dI1/dE = 2*identity and dI2/dE = 2*(I1*Idenity - C) and dI3/dE = 2*I3*inv(C)
				T = (2*dW_by_dI1*identity)  +  (2*dW_by_dI2*(I1*identity - C))  +  (2*dW_by_dI3 * I3 * C.Inverse());
				return T;
			}	
			break;
		
			default:
			{
				assert(0);
				//return F; // Avoid compiler warning!
			}	
		}
        return F;
	}
	
	FourthOrderTensor<SPACE_DIM> ComputeIsotropic_dTdE(MatrixDouble F)
	{
		assert(F.Rows()==SPACE_DIM);

		FourthOrderTensor<SPACE_DIM> dTdE;
		

		FourthOrderTensor<SPACE_DIM> IdentityCrossIdentity;
		IdentityCrossIdentity.SetAsTensorProductOfIdentities();		
		
		
		MatrixDouble    C = (F.Transpose())*F;
		MatrixDouble invC = C.Inverse();

		switch(SPACE_DIM)
		{
			case 1:
			{	
				double I1 = C.GetFirstInvariant();
				//double  dW_by_dI1  = GetdW_by_dI1(I1);	
				double  d2W_by_dI1 = Getd2W_by_dI1(I1);
				
				for(int M=0; M<SPACE_DIM; M++)
				{
					for(int N=0; N<SPACE_DIM; N++)
					{	
						for(int P=0; P<SPACE_DIM; P++)
						{
							for(int Q=0; Q<SPACE_DIM; Q++)
							{
								dTdE.rGetVal()[M][N][P][Q] = 4*d2W_by_dI1*IdentityCrossIdentity.rGetVal()[M][N][P][Q];
							}
						}
					}
				}				
			}	
			break;	
			
			case 2:
			{	

				double I1 = C.GetFirstInvariant();
				double I2 = C.GetSecondInvariant();

				//double   dW_by_dI1   =    GetdW_by_dI1(I1,I2);
				double   dW_by_dI2   =    GetdW_by_dI2(I1,I2);
				double  d2W_by_dI1   =   Getd2W_by_dI1(I1,I2);
				double  d2W_by_dI2   =   Getd2W_by_dI2(I1,I2);
				double  d2W_by_dI1dI2 = Getd2W_by_dI1dI2(I1,I2);
				
				MatrixDouble dI1_dE = 2*MatrixDouble::Identity(2);
				MatrixDouble dI2_dE = 2 * I2 * C.Inverse();
				
				// note: could have FourthOrderTensor dI1_dEdE, this is equal to zero
				FourthOrderTensor<SPACE_DIM> dI2_dEdE;
				
				for(int M=0; M<SPACE_DIM; M++)
				{
					for(int N=0; N<SPACE_DIM; N++)
					{	
						for(int P=0; P<SPACE_DIM; P++)
						{
							for(int Q=0; Q<SPACE_DIM; Q++)
							{
								dI2_dEdE.rGetVal()[M][N][P][Q] = -4*I2*invC(M,P)*invC(Q,N) + 2*invC(M,N)*dI2_dE(P,Q);
							}
						}
					}
				}
				
				
				for(int M=0; M<SPACE_DIM; M++)
				{
					for(int N=0; N<SPACE_DIM; N++)
					{	
						for(int P=0; P<SPACE_DIM; P++)
						{
							for(int Q=0; Q<SPACE_DIM; Q++)
							{
								dTdE.rGetVal()[M][N][P][Q] =   d2W_by_dI1    * dI1_dE(M,N) * dI1_dE(P,Q)
								                        + d2W_by_dI1dI2 * dI1_dE(M,N) * dI2_dE(P,Q)
								                        + 0 // dW_by_dI1 * dI1_dEdE
								                        
								                        + d2W_by_dI2    * dI2_dE(M,N) * dI2_dE(P,Q)
								                        + d2W_by_dI1dI2 * dI2_dE(M,N) * dI1_dE(P,Q)
								                        + dW_by_dI2     * dI2_dEdE.rGetVal()[M][N][P][Q];
							}
						}
					}
				}	
		
			}
			break;

			case 3:
			{		
				double I1 = C.GetFirstInvariant();
				double I2 = C.GetSecondInvariant();
				double I3 = C.GetThirdInvariant();

				//double dW_by_dI1 = GetdW_by_dI1(I1,I2,I3);		
				double dW_by_dI2 = GetdW_by_dI2(I1,I2,I3);		
				//double dW_by_dI3 = GetdW_by_dI2(I1,I2,I3);		

				double  d2W_by_dI1   =   Getd2W_by_dI1(I1,I2,I3);
				double  d2W_by_dI2   =   Getd2W_by_dI2(I1,I2,I3);
				//double  d2W_by_dI3   =   Getd2W_by_dI3(I1,I2,I3);

				double  d2W_by_dI1dI2 = Getd2W_by_dI1dI2(I1,I2,I3);
				double  d2W_by_dI1dI3 = Getd2W_by_dI2dI3(I1,I2,I3);
				double  d2W_by_dI2dI3 = Getd2W_by_dI2dI3(I1,I2,I3);

				MatrixDouble dI1_dE = 2 * (MatrixDouble::Identity(3));	
				MatrixDouble dI2_dE = 2 * (I1 * MatrixDouble::Identity(3) - C);
				MatrixDouble dI3_dE = 2 * I3 * invC;

				// note: could have FourthOrderTensor dI1_dEdE, this is equal to zero
				FourthOrderTensor<SPACE_DIM> dI2_dEdE;
				FourthOrderTensor<SPACE_DIM> dI3_dEdE;

				for(int M=0; M<SPACE_DIM; M++)
				{
					for(int N=0; N<SPACE_DIM; N++)
					{	
						for(int P=0; P<SPACE_DIM; P++)
						{
							for(int Q=0; Q<SPACE_DIM; Q++)
							{
								dI2_dEdE.rGetVal()[M][N][P][Q] =  4*(M==N)*(P==Q) - 4*(M==P)*(N==Q);
								dI3_dEdE.rGetVal()[M][N][P][Q] = -4*I3*invC(M,P)*invC(Q,N) + 2*invC(M,N)*dI3_dE(P,Q);
							}
						}
					}
				}


				for(int M=0; M<SPACE_DIM; M++)
				{
					for(int N=0; N<SPACE_DIM; N++)
					{	
						for(int P=0; P<SPACE_DIM; P++)
						{
							for(int Q=0; Q<SPACE_DIM; Q++)
							{
								dTdE.rGetVal()[M][N][P][Q] =   d2W_by_dI1    * dI1_dE(M,N) * dI1_dE(P,Q)
								                        + d2W_by_dI1dI2 * dI1_dE(M,N) * dI2_dE(P,Q)
								                        + d2W_by_dI1dI3 * dI1_dE(M,N) * dI3_dE(P,Q)
								                        + 0 // dW_by_dI1 * dI1_dEdE
								                        
								                        + d2W_by_dI2    * dI2_dE(M,N) * dI2_dE(P,Q)
								                        + d2W_by_dI1dI2 * dI2_dE(M,N) * dI1_dE(P,Q)
								                        + d2W_by_dI2dI3 * dI2_dE(M,N) * dI3_dE(P,Q)
								                        + dW_by_dI2     * dI2_dEdE.rGetVal()[M][N][P][Q];

								                        //+ d2W_by_dI3    * dI3_dE(M,N) * dI3_dE(P,Q)
								                        //+ d2W_by_dI1dI3 * dI3_dE(M,N) * dI1_dE(P,Q)
								                        //+ d2W_by_dI2dI3 * dI3_dE(M,N) * dI2_dE(P,Q)
								                        //+ dW_by_dI3     * dI3_dEdE.rGetVal()[M][N][P][Q];
							}
						}
					}
				}
				
			}	
			break;
		
			default:
			{
				assert(0);
			}	
		}
		
		return dTdE;
	}	
};


#endif //_ABSTRACTMATERIAL_HPP_
