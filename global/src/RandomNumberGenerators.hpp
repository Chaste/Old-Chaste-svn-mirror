#ifndef RANDOMNUMBERGENERATORS_HPP_
#define RANDOMNUMBERGENERATORS_HPP_
#include <cmath>
#include <time.h>
#include <stdlib.h>

class RandomNumberGenerators
{
public:
    double StandardNormalRandomDeviate(void);
    double NormalRandomDeviate(double mean, double sd);
    double ranf(void);
    int randMod(int base);
    
    /**
     * @param seed Is the new seed which defaults to zero.
     */
    RandomNumberGenerators(int seed=0)
    {
        srandom(seed);
    }
    
};
#endif /*RANDOMNUMBERGENERATORS_HPP_*/
