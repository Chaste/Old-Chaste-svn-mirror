#ifndef RANDOMNUMBERGENERATORS_HPP_
#define RANDOMNUMBERGENERATORS_HPP_
#include <cmath>
#include <time.h>
#include <stdlib.h>

class RandomNumberGenerator
{
public:
    double StandardNormalRandomDeviate(void);
    double NormalRandomDeviate(double mean, double sd);
    double ranf(void);
    unsigned randMod(unsigned base);
    
 
    
    static RandomNumberGenerator* Instance();
    static void Destroy();
	void Reseed(int seed)
	{
		srandom(seed);
	}
protected:
   /**
     * @param seed Is the new seed which defaults to zero.
     */
    RandomNumberGenerator()
    {
        srandom(0);
    }
    
private:


    static RandomNumberGenerator* mpInstance;
};
#endif /*RANDOMNUMBERGENERATORS_HPP_*/
