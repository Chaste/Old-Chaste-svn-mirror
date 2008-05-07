/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef ABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_
#define ABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_

#include "FiniteElasticityAssembler.cpp"

// todo: lots of doxygen, refactor to take in mesh at the beginning?


template<unsigned DIM>
struct EvaluationPointInfo
{
    Point<DIM> OldPosition;
    Point<DIM> NewPosition;
    double SourceValue;
};


template<unsigned DIM>
class AbstractGrowingTumourSourceModel
{
private:
    friend class TestAbstractGrowingTumourSourceModel;
    
protected :
    std::map< unsigned, EvaluationPointInfo<DIM> > mEvaluationPoints;
    std::vector<int> mVertexToEvalPointMap;
        
public :
    AbstractGrowingTumourSourceModel()
    {
    }
    
    virtual ~AbstractGrowingTumourSourceModel()
    {}
        
    virtual void Run(double tStart, double tEnd, FiniteElasticityAssembler<DIM>* pFiniteElasticityAssembler)=0;

    
    virtual double GetSourceValue(unsigned vertexIndex)
    {
        // could just do:
        //  return mEvaluationPoints[index].SourceValue;
        // but this doesn't check for existence

        typename std::map< unsigned, EvaluationPointInfo<DIM> >::iterator iter
        = mEvaluationPoints.find(vertexIndex);
        
        if (iter == mEvaluationPoints.end())
        {
            std::stringstream ss;
            ss << "No evaluation point corresponding to index=" << vertexIndex;
            EXCEPTION(ss.str());
        }
        
        return iter->second.SourceValue;
    }
    
    

    void AddEvaluationPoint(unsigned vertexIndex, Point<DIM> initialPosition)
    {
        // check an evaluation point with this index doesn't already 
        // exist

        typename std::map< unsigned, EvaluationPointInfo<DIM> >::iterator iter
            = mEvaluationPoints.find(vertexIndex);
        if (iter != mEvaluationPoints.end())
        {
            std::stringstream ss;
            ss << "An evaluation point corresponding to index=" << vertexIndex << " already exists";
            EXCEPTION(ss.str());
        }
        
        EvaluationPointInfo<DIM> evaluation_point_info; 
        mEvaluationPoints[vertexIndex] = evaluation_point_info;
        
        mEvaluationPoints[vertexIndex].OldPosition = initialPosition;
        mEvaluationPoints[vertexIndex].NewPosition = initialPosition;
        mEvaluationPoints[vertexIndex].SourceValue = 0.0;
    }
    
    

    void UpdateEvaluationPointsNewPosition(std::vector<Vector<double> >& deformedPositions)
    {
        typename std::map<unsigned,EvaluationPointInfo<DIM> >::iterator iter 
           = mEvaluationPoints.begin();
        
        while(iter!=mEvaluationPoints.end())
        {
            unsigned vertex_index = iter->first;
            for(unsigned i=0; i<DIM; i++)
            {
                iter->second.NewPosition[i] = deformedPositions[i](vertex_index);
            }
            iter++;
        }
    }
        
    
    unsigned GetNumEvaluationPoints()
    {
        return mEvaluationPoints.size();
    }
};

#endif /*ABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_*/
