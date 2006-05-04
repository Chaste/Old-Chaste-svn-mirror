/**
 * Concrete Ode3 class
 */ 
#ifndef _ODE3_HPP_
#define _ODE3_HPP_
#include "AbstractOdeSystem.hpp"


class Ode3 : public AbstractOdeSystem
{
	public :
	Ode3();
	
	std::vector<double> EvaluateYDerivatives (double time, const std::vector<double> &rY);

};

#endif //_ODE3_HPP_
