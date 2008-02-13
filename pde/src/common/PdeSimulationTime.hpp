#ifndef PDESIMULATIONTIME_HPP_
#define PDESIMULATIONTIME_HPP_
class PdeSimulationTime
{
public:
    static void SetTime(double time);
    static double GetTime();
private:
    static double mTime;
};        

#endif /*PDESIMULATIONTIME_HPP_*/
