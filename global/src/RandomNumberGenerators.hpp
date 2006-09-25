#ifndef RANDOMNUMBERGENERATORS_HPP_
#define RANDOMNUMBERGENERATORS_HPP_

class RandomNumberGenerators
{ 
public:
	double StandardNormalRandomDeviate(void);
	double NormalRandomDeviate(double mean, double sd);
	double ranf(void);
	int randMod(int base);
	RandomNumberGenerators()
	{
		mSeeded = false;
	}
private:
	bool mSeeded;
};
#endif /*RANDOMNUMBERGENERATORS_HPP_*/
