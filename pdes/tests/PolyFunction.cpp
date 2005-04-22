#include "PolyFunction.hpp"
#include <cmath>


	PolyFunction::PolyFunction(int degree)
	{
		mDegree=degree;
	}
	
	double PolyFunction::Evaluate(Point<1> evaluation_point)
	{
		return pow(evaluation_point[0], mDegree);
	}

