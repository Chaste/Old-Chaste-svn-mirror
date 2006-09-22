#ifndef SIMULATIONTIME_HPP_
#define SIMULATIONTIME_HPP_

class SimulationTime 
{
public:
	static SimulationTime* Instance(double, int);
	static SimulationTime* Instance();
	double GetTimeStep();
	void IncrementTimeOneStep();
	int GetTimeStepsElapsed();
	double GetDimensionalisedTime();
	static void Destroy();
protected:
 	SimulationTime(double, int);
private:
 	static SimulationTime* mInstance;
 	double mDurationOfSimulation;
 	int mTotalTimeStepsInSimulation;
 	int mTimeStepsElapsed;
};

#endif /*SIMULATIONTIME_HPP_*/
