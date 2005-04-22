#ifndef _POLYFUNCTION_HPP_
#define _POLYFUNCTION_HPP_

#include "AbstractIntegrand.hpp"

class PolyFunction : public AbstractIntegrand<1> {
	private:
	int mDegree;
	
	public:
	PolyFunction(int degree);
	
	double Evaluate(Point<1> evaluation_point);
};

#endif //_POLYFUNCTION_HPP_
