#ifndef SIMULATIONTIME_HPP_
#define SIMULATIONTIME_HPP_

class SimulationTime 
{
public:
	static SimulationTime* Instance();
	double GetTime();
	void IncrementTime(double);
protected:
 	SimulationTime();
private:
 	static SimulationTime* mInstance;
 	double mTime;
};

#endif /*SIMULATIONTIME_HPP_*/
