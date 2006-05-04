/**
 * Concrete Ode2 class
 */ 
#ifndef _ODE2_HPP_
#define _ODE2_HPP_
#include "AbstractOdeSystem.hpp"


class Ode2 : public AbstractOdeSystem
{
	public :
	Ode2();
	
	std::vector<double> EvaluateYDerivatives (double time, const std::vector<double> &rY);

};

#endif //_ODE2_HPP_
