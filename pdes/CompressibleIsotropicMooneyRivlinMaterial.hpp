#ifndef _COMPRESSIBLEISOTROPICMOONEYRIVLINMATERIAL_HPP_
#define _COMPRESSIBLEISOTROPICMOONEYRIVLINMATERIAL_HPP_

template <int SPACE_DIM>
class CompressibleIsotropicMooneyRivlinMaterial : public AbstractMaterial<SPACE_DIM>
{
public:
	CompressibleIsotropicMooneyRivlinMaterial(double c1, double c2=0)
	{
		AbstractMaterial<SPACE_DIM>::mDensitySet = false;
		AbstractMaterial<SPACE_DIM>::mIsIsotropicLaw = true;

		mC1 = c1;	
		
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
	
	double GetdW_by_dE(MatrixDouble E, int index1, int index2)
	{
		// this is an isotropic law so this function is not implemented
		assert(0);
	}
	
};

#endif //_COMPRESSIBLEISOTROPICMOONEYRIVLINMATERIAL_HPP_
