#ifndef PROBLEMPORTION_HPP_
#define PROBLEMPORTION_HPP_

#include "ParallelIterator.hpp"

class ParallelProblem
{
private:
    unsigned mLo;
    unsigned mHi;
    unsigned mProblemSize;
public:
    void SetProblemSize(unsigned problemSize);
    unsigned Size();
    unsigned GetProblemSize();
    ParallelIterator Begin();
    ParallelIterator End();
    unsigned Global(ParallelIterator iterator);
    unsigned Local(ParallelIterator iterator);
};

#endif /*PROBLEMPORTION_HPP_*/
