#ifndef _COMPRESSIBLEISOTROPICMOONEYRIVLINMATERIAL_HPP_
#define _COMPRESSIBLEISOTROPICMOONEYRIVLINMATERIAL_HPP_

/** CompressibleIsotropicMooneyRivlinMaterial
 *  
 *  Implementation of an compressible Mooney-Rivlin type material
 *  
 *  In 2D:
 *  This is strain energy W(I1,I2) = c1(I1-2) + c2(I2-1)
 *  The stress at zero strain must be zero, so the user can only specify the first 
 *  parameter c1 (in the constructor), and the second parameter is computed 
 *  accordingly
 *  
 *  This is W(I1,I2) = c1(I1-3) + c2(I2-3) + c3*(I3-1)
 *  The stress at zero strain must be zero, so the user can only specify the first 
 *  and second parameters c1 and c2 (in the constructor), and the third parameter is 
 *  computed accordingly. 
 * 
 *  The 1D equivalent of this cannot satisfy zero stress at zero strain, so this
 *  isn't a valid strain energy and the user is not allowed to instantiated a 1D 
 *  version of this class
 */
template <int SPACE_DIM>
class CompressibleIsotropicMooneyRivlinMaterial : public AbstractMaterial<SPACE_DIM>
{
public:
	CompressibleIsotropicMooneyRivlinMaterial(double c1, double c2=0)
	{
		AbstractMaterial<SPACE_DIM>::mDensitySet = false;

		/** Set as isotropic law
		 */
		AbstractMaterial<SPACE_DIM>::mIsIsotropicLaw = true;

		mC1 = c1;	
		
		/** This isn't a valid strain energy in 1D
		 */
		if(SPACE_DIM==1)
		{
			assert(0);
		}
		/** set up last parameter so that the stress is zero when strain is zero
		 */
		if(SPACE_DIM==2)
		{
			mC2 = -mC1;
		}		
		if(SPACE_DIM==3)
		{
			mC2 = c2;	
			mC3 = - mC1 - (SPACE_DIM-1)*mC2;
		}
	}	


private:
	double mC1;
	double mC2;
	double mC3;
	
	double GetdW_by_dI1   (double I1,double I2=0,double I3=0)
	{
		return mC1;
	}
    double GetdW_by_dI2   (double I1,double I2=0,double I3=0)
    {
    	assert(SPACE_DIM>1);
		return mC2;
	}
    double GetdW_by_dI3   (double I1,double I2=0,double I3=0)
    {
    	assert(SPACE_DIM>2);
		return mC3;
	}
    

	double Getd2W_by_dI1  (double I1,double I2=0,double I3=0)
	{
		return 0;
	}
	
    double Getd2W_by_dI2  (double I1,double I2=0,double I3=0)
	{
    	assert(SPACE_DIM>1);
		return 0;
	}

    double Getd2W_by_dI3  (double I1,double I2=0,double I3=0)
   	{
    	assert(SPACE_DIM>2);
		return 0;
	}

    double Getd2W_by_dI1I2(double I1,double I2=0,double I3=0)
	{
    	assert(SPACE_DIM>1);
		return 0;
	}

    double Getd2W_by_dI2I3(double I1,double I2=0,double I3=0)
	{
    	assert(SPACE_DIM>2);
		return 0;
	}

    double Getd2W_by_dI1I3(double I1,double I2=0,double I3=0)
	{
    	assert(SPACE_DIM>2);
		return 0;
	}	
	
};

#endif //_COMPRESSIBLEISOTROPICMOONEYRIVLINMATERIAL_HPP_
