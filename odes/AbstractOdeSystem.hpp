#ifndef _ABSTRACTODESYSTEM_HPP_
#define _ABSTRACTODESYSTEM_HPP_

//#include "petscvec.h"
#include <vector>
// AbstractOdeSystem.hpp

class AbstractOdeSystem
{
	public:
	int mNumberOfEquations;
	
	AbstractOdeSystem(const int& rNumberOfEquations);
	~AbstractOdeSystem();
	
	virtual std::vector<double> EvaluateYDerivatives (double rTime, std::vector<double> &rY) = 0;
	
	
};

#endif //_ABSTRACTODESYSTEM_HPP_
