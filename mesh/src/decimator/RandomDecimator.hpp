#ifndef _RANDOM_DECIMATOR_HPP_
#define _RANDOM_DECIMATOR_HPP_


#include "Decimator.hpp"
#include "RandomNumberGenerator.hpp"
///Note ELEMENT_DIM matches SPACE_DIM
template <int SPACE_DIM>
class RandomDecimator : public Decimator<SPACE_DIM>
{
    RandomNumberGenerator * mpRanGen;
protected:

    void CalculateLocalMeasure(NodeInfo<SPACE_DIM> *pNodeInfo, bool before)
    {
        if (before)
        {
            this->mMeasureBefore=mpRanGen->ranf();
        } else {
            this->mMeasureAfter=this->mMeasureBefore;
        }

    }
    
    //double CalculateScore returns mMeasureAfter as in the base class
public:
	void Initialise(ConformingTetrahedralMesh<SPACE_DIM, SPACE_DIM> *pMesh, RandomNumberGenerator *pRanGen=NULL)
    {
        if (pRanGen)
        {
        	mpRanGen=pRanGen;
        }
        else
        {
         	mpRanGen=new RandomNumberGenerator;
        }
        Decimator<SPACE_DIM>::Initialise(pMesh);
    }
};


#endif //_RANDOM_DECIMATOR_HPP_
