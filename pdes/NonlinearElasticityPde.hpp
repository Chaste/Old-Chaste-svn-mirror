#ifndef _NONLINEARELASTICITYPDE_HPP_
#define _NONLINEARELASTICITYPDE_HPP_

#include "AbstractMaterial.hpp"
#include "Element.hpp"
#include "MatrixDouble.hpp"
#include "VectorDouble.hpp"
#include "FourthOrderTensor.hpp"

template <int ELEM_DIM, int SPACE_DIM>
class NonlinearElasticityPde
{
private:
	VectorDouble* mpGravity;

public:
	NonlinearElasticityPde(VectorDouble gravity)
	{
		assert(gravity.Size()==SPACE_DIM);
		mpGravity = new VectorDouble(SPACE_DIM);
		
		for(int i=0;i<SPACE_DIM;i++)
		{
			(*mpGravity)(i) = gravity(i);
		}
	}


	MatrixDouble ComputeStress(const Element<ELEM_DIM, SPACE_DIM> & rElement, const MatrixDouble F)
	{
		return rElement.GetMaterial()->ComputeStress(F);
	}

	FourthOrderTensor<SPACE_DIM> Compute_dTdE(const Element<ELEM_DIM, SPACE_DIM> & rElement, const MatrixDouble F)
	{
		return rElement.GetMaterial().Compute_dTdE(F);
	}	
	
	VectorDouble ComputeGravityForceTerm(const Element<ELEM_DIM, SPACE_DIM> & rElement)
	{
		return rElement.GetMaterial()->GetDensity() * (*mpGravity);
	}
};

#endif //_NONLINEARELASTICITYPDE_HPP_
