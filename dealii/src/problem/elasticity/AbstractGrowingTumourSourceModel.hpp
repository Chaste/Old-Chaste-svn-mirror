#ifndef ABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_
#define ABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_

#include "FiniteElasticityAssembler.cpp"


template<int DIM>
struct EvaluationPointInfo
{
    Point<DIM> OldPosition;
    Point<DIM> NewPosition;
    unsigned MeshGlobalIndex;
    double SourceValue;
}; 


template<int DIM>
class AbstractGrowingTumourSourceModel
{
protected :
    std::map< unsigned, EvaluationPointInfo<DIM> > mEvaluationPoints; 

//    std::vector<bool> mVertexToEvalPointMap;
    
    
public :
    virtual void Run(double tStart, double tEnd)=0;

    
    virtual double GetSourceValue(unsigned index)
    {
        // could just do:
        //  return mEvaluationPoints[index].SourceValue;
        // but this doesn't check for existence
                
        typename std::map< unsigned, EvaluationPointInfo<DIM> >::iterator iter 
           = mEvaluationPoints.find(index);
           
        if(iter == mEvaluationPoints.end())
        {
            std::stringstream ss;
            ss << "No evaluation point corresponding to index=" << index;
            EXCEPTION(ss.str());
        }
        
        return iter->second.SourceValue;
    }
    

    void AddEvaluationPoint(unsigned index, Point<DIM> initialPosition, unsigned meshGlobalIndex)
    {
        // check an evaluation point with this index doesn't already 
        // exist
        typename std::map< unsigned, EvaluationPointInfo<DIM> >::iterator iter 
           = mEvaluationPoints.find(index);
        if(iter != mEvaluationPoints.end())
        {
            std::stringstream ss;
            ss << "An evaluation point corresponding to index=" << index << " already exists";
            EXCEPTION(ss.str());
        }             
        
        EvaluationPointInfo<DIM> evaluation_point_info;
        mEvaluationPoints[index] = evaluation_point_info;
 
        mEvaluationPoints[index].OldPosition = initialPosition;
        mEvaluationPoints[index].NewPosition = initialPosition;
        mEvaluationPoints[index].MeshGlobalIndex = meshGlobalIndex;
        mEvaluationPoints[index].SourceValue = 0.0;
    }


    void UpdateEvaluationPointsNewPosition(FiniteElasticityAssembler<DIM>* pFiniteElasticityAssembler)
    {
//        // loop over tumour points and update new position
//        Triangulation<DIM>* p_mesh = pFiniteElasticityAssembler->GetMesh();
//
//        // this isn't efficient is it?!
//        mVertexToEvalPointMap.resize(p_mesh->n_vertices());
//        for(unsigned vertex_index=0; vertex_index<mVertexToEvalPointMap.size(); vertex_index++)
//        {
//            mVertexToEvalPointMap[vertex_index] = -1;
//        }
//
//        typename std::map<unsigned,EvaluationPointInfo<DIM> >::iterator iter 
//           = mEvaluationPoints.begin();
//        while(iter!=mEvaluationPoints.end())
//        {
//            unsigned eval_pt_index = iter->first;
//            unsigned vertex_index = iter->second.MeshGlobalIndex;
//            mVertexToEvalPointMap[vertex_index] = eval_pt_index;
// 
//            iter++;
//        }
//
//        Vector<double>& solution = pFiniteElasticityAssembler->GetSolutionVector();
//        DoFHandler<DIM>& dof_handler = pFiniteElasticityAssembler->GetDofHandler();
//
//        DofVertexIterator<DIM> vertex_iter(p_mesh, &dof_handler);
//        
//        while(!vertex_iter.ReachedEnd())
//        {
//            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
//            if(mVertexToEvalPointMap[vertex_index]!=-1)
//            {
//                unsigned eval_pt_index = mVertexToEvalPointMap[vertex_index];
//                
//                Point<DIM> old_posn = vertex_iter.GetVertex();
//            
//                for(unsigned i=0; i<DIM; i++)
//                {
//                    mEvaluationPoints[eval_pt_index].NewPosition[i]
//                        = old_posn(i)+solution(vertex_iter.GetDof(i));
//                }
//            }
//            vertex_iter.Next();
//        }
    }


    
    unsigned GetNumEvaluationPoints()
    {
        return mEvaluationPoints.size();
    }
};

#endif /*ABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_*/
