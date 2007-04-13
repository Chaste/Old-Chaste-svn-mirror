#ifndef PROBLEMPORTION_HPP_
#define PROBLEMPORTION_HPP_

#include "ParallelIterator.hpp"

class ParallelProblem
{
private:
    unsigned mLo;
    unsigned mHi;
public:
    void SetProblemSize(unsigned problemSize);
    unsigned Size();
    ParallelIterator Begin();
    ParallelIterator End();
    unsigned Global(ParallelIterator iterator);
    unsigned Local(ParallelIterator iterator);
};

#endif /*PROBLEMPORTION_HPP_*/
