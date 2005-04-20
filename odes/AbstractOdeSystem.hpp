#ifndef _ABSTRACTODESYSTEM_HPP_
#define _ABSTRACTODESYSTEM_HPP_

#include "petscvec.h"

// AbstractOdeSystem.hpp

class AbstractOdeSystem
{
	public:
	int mNumberOfEquations;
	double mTInit;
	double * mYInit;
	
	AbstractOdeSystem(const int& rNumberOfEquations,const double& rTInit, double * rYInit);
	~AbstractOdeSystem();
	
	virtual void EvaluateYPrime (double rTime, double * rY, double * rYPrime) = 0;
	
	
};

#endif //_ABSTRACTODESYSTEM_HPP_
