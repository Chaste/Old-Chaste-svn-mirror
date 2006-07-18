#ifndef _NONLINEARELASTICITYPDE_HPP_
#define _NONLINEARELASTICITYPDE_HPP_

#include "AbstractMaterial.hpp"
#include "Element.cpp"
#include "FourthOrderTensor.hpp"

template <int ELEM_DIM, int SPACE_DIM>
class NonlinearElasticityPde
{
private:
	c_vector<double, SPACE_DIM> mpGravity;

public:
	NonlinearElasticityPde(c_vector<double, SPACE_DIM> gravity)
	{
		assert(gravity.size()==SPACE_DIM);
		mpGravity = gravity;
	}
	
 //\todo:  Ticket 96: we need to use ublas here (and compile-time mumble)    
//	MatrixDouble ComputeStress(const Element<ELEM_DIM, SPACE_DIM> & rElement, const MatrixDouble F)
//	{
//		return rElement.GetMaterial()->ComputeStress(F);
//	}

//	FourthOrderTensor<SPACE_DIM> Compute_dTdE(const Element<ELEM_DIM, SPACE_DIM> & rElement, const MatrixDouble F)
//	{
//		return rElement.GetMaterial().Compute_dTdE(F);
//	}	
	
	c_vector<double, SPACE_DIM> ComputeGravityForceTerm(const Element<ELEM_DIM, SPACE_DIM> & rElement)
	{
		return rElement.GetMaterial()->GetDensity() * mpGravity;
	}
};

#endif //_NONLINEARELASTICITYPDE_HPP_
