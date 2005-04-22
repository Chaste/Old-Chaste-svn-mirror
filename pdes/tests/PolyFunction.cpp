#include "PolyFunction.hpp"
#include <cmath>


	poly_function::poly_function(int degree)
	{
		mDegree=degree;
	}
	
	double poly_function::Evaluate(Point<1> evaluation_point)
	{
		return pow(evaluation_point[0], mDegree);
	}

