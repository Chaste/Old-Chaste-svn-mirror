#ifndef _VARYINGDIFFUSIONANDSOURCETERMPDE_
#define _VARYINGDIFFUSIONANDSOURCETERMPDE_

#include "AbstractLinearEllipticPde.hpp"
#include "VectorDouble.hpp"
#include "Point.hpp"
#include <cmath>

template <int SPACE_DIM>
class VaryingDiffusionAndSourceTermPde:public AbstractLinearEllipticPde<SPACE_DIM>
{
	private:
	double DistanceFromOrigin(Point<SPACE_DIM> x)
	{
		double sum=0;
		for (int i=0; i<SPACE_DIM; i++)
		{
			sum += x[i]*x[i];
		}
	return sqrt(sum);
	}
	
	public:
	double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
	{
		return 0.0;
	}
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x, double u)
    {
    	return -pow(DistanceFromOrigin(x),3);
    }

    MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x)
    {
    	return pow(DistanceFromOrigin(x),2)*MatrixDouble::Identity(SPACE_DIM);
    }
};


#endif //_VARYINGDIFFUSIONANDSOURCETERMPDE_
