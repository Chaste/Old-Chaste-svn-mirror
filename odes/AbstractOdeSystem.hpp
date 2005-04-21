#ifndef _ABSTRACTODESYSTEM_HPP_
#define _ABSTRACTODESYSTEM_HPP_

#include "petscvec.h"

// AbstractOdeSystem.hpp

class AbstractOdeSystem
{
	public:
	int mNumberOfEquations;
	// Gary's dodgy update
	//double mTInit;
	//double * mYInit;
	
	AbstractOdeSystem(const int& rNumberOfEquations);
	~AbstractOdeSystem();
	
	virtual void EvaluateYDerivatives (double rTime, double * rY, double * rYDerivatives) = 0;
	
	
};

#endif //_ABSTRACTODESYSTEM_HPP_
