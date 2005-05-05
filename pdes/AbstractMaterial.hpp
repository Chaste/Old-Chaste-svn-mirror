#ifndef _ABSTRACTMATERIAL_HPP_
#define _ABSTRACTMATERIAL_HPP_

#include "MatrixDouble.hpp"
#include <vector>

typedef std::vector< std::vector < std::vector< std::vector<double> > > > FourthOrderTensor;

template <int SPACE_DIM>
class AbstractMaterial 
{

// functions which need to be implemented by derived classes	
private:
    virtual double GetdW_by_dI1   (double I1,double I2=0,double I3=0)=0;
    virtual double GetdW_by_dI2   (double I1,double I2=0,double I3=0)=0;
    virtual double GetdW_by_dI3   (double I1,double I2=0,double I3=0)=0;

    virtual double Getd2W_by_dI1  (double I1,double I2=0,double I3=0)=0;
    virtual double Getd2W_by_dI2  (double I1,double I2=0,double I3=0)=0;
    virtual double Getd2W_by_dI3  (double I1,double I2=0,double I3=0)=0;
    virtual double Getd2W_by_dI1I2(double I1,double I2=0,double I3=0)=0;
    virtual double Getd2W_by_dI2I3(double I1,double I2=0,double I3=0)=0;
    virtual double Getd2W_by_dI1I3(double I1,double I2=0,double I3=0)=0;

    virtual double GetdW_by_dE(MatrixDouble E, int index1, int index2)=0;
    
    
public:
    double GetDensity()
    {
		if(!mDensitySet)
		{ 
			std::cout << std::endl << "Error: density not set" << std::endl;
			assert(0);
		}
		return mDensity;
    }
    
    void SetDensity(double density)
	{
		mDensity = density;
		mDensitySet = true;
	}
    
    
    MatrixDouble ComputeStress(MatrixDouble F)
    {
    	if(mIsIsotropicLaw)
		{
			return ComputeIsotropicStress(F);
		}
		else
		{
			assert(0);
		}
	}
	
	
	FourthOrderTensor Compute_dTdE (MatrixDouble F)
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
					
				T = (2*dW_by_dI1*identity)  +  (2*dW_by_dI2*(I1*identity - C))  +  (2*dW_by_dI3 * I3 * C.Inverse());
				return T;
			}	
			break;
		
			default:
			{
				assert(0);
			}	
		}
	}
	
	FourthOrderTensor ComputeIsotropic_dTdE(MatrixDouble F)
	{
		assert(F.Rows()==SPACE_DIM);

		FourthOrderTensor dTdE;
		
		dTdE.resize(SPACE_DIM);		
		for(M=0; M<SPACE_DIM; M++)
		{
			dTdE[M].resize(SPACE_DIM);
			for(N=0; N<SPACE_DIM; N++)
			{	
				dTdE[M][N].resize(SPACE_DIM);
				for(P=0; P<SPACE_DIM; P++)
				{
					dTdE[M][N][P].resize(SPACE_DIM);
				}
			}
		}				


		for(M=0; M<SPACE_DIM; M++)
		{
			for(N=0; N<SPACE_DIM; N++)
			{	
				for(P=0; P<SPACE_DIM; P++)
				{
					for(Q=0; Q<SPACE_DIM; Q++)
					{
						dTdE[M][N][P][Q] = 0;
					}
				}
			}
		}				

		MatrixDouble C = (F.Transpose())*F;

		switch(SPACE_DIM)
		{
			case 1:
			{	
				double I1 = C.GetFirstInvariant();
				double  dW_by_dI1  = GetdW_by_dI1(I1);	
				double  d2W_by_dI1 = Getd2W_by_dI1(I1);
				
				for(M=0; M<SPACE_DIM; M++)
				{
					for(N=0; N<SPACE_DIM; N++)
					{	
						for(P=0; P<SPACE_DIM; P++)
						{
							for(Q=0; Q<SPACE_DIM; Q++)
							{
								dTdE[M][N][P][Q] = 4*d2W_by_dI1*(M==N)*(P==Q);
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

				double   dW_by_dI1   =    GetdW_by_dI1(I1,I2);
				double   dW_by_dI2   =    GetdW_by_dI2(I1,I2);
				double  d2W_by_dI1   =   Getd2W_by_dI1(I1,I2);
				double  d2W_by_dI2   =   Getd2W_by_dI2(I1,I2);
				double  d2W_by_dI1I2 = Getd2W_by_dI1I2(I1,I2);
				
				// fill in
				assert(0);				
			}
			break;

			case 3:
			{		
				double I1 = C.GetFirstInvariant();
				double I2 = C.GetSecondInvariant();
				double I3 = C.GetThirdInvariant();

				double dW_by_dI1 = GetdW_by_dI1(I1,I2,I3);		
				double dW_by_dI2 = GetdW_by_dI2(I1,I2,I3);		
				double dW_by_dI3 = GetdW_by_dI2(I1,I2,I3);		

				double  d2W_by_dI1   =   Getd2W_by_dI1(I1,I2,I3);
				double  d2W_by_dI2   =   Getd2W_by_dI2(I1,I2,I3);
				double  d2W_by_dI3   =   Getd2W_by_dI3(I1,I2,I3);

				double  d2W_by_dI1I2 = Getd2W_by_dI1I2(I1,I2,I3);
				double  d2W_by_dI1I3 = Getd2W_by_dI2I3(I1,I2,I3);
				double  d2W_by_dI2I3 = Getd2W_by_dI2I3(I1,I2,I3);

				// fill in
				assert(0);						
			}	
			break;
		
			default:
			{
				assert(0);
			}	
		}
	}	
};


#endif //_ABSTRACTMATERIAL_HPP_
