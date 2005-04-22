#ifndef _POLYFUNCTION_HPP_
#define _POLYFUNCTION_HPP_

#include "AbstractIntegrand.hpp"

class poly_function : public AbstractIntegrand<1> {
	private:
	int mDegree;
	
	public:
	poly_function(int degree);
	
	double Evaluate(Point<1> evaluation_point);
};

#endif //_POLYFUNCTION_HPP_
