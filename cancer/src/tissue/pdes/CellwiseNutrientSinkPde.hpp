#ifndef CELLWISENUTRIENTSINKPDE_HPP_
#define CELLWISENUTRIENTSINKPDE_HPP_

#include "MeshBasedTissue.hpp"
#include "AbstractLinearEllipticPde.hpp"

/**
 *  A nutrient PDE which has a sink at each non-necrotic cell.
 */
template<unsigned DIM>
class CellwiseNutrientSinkPde : public AbstractLinearEllipticPde<DIM>
{
private:

    MeshBasedTissue<DIM>& mrTissue;
    
    double mCoefficient;
    
public:

    CellwiseNutrientSinkPde(MeshBasedTissue<DIM>& rTissue, double coefficient);

    double ComputeConstantInUSourceTerm(const ChastePoint<DIM>& x);
    
    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<DIM>& x, Element<DIM,DIM>*);
   
    double ComputeLinearInUCoeffInSourceTermAtNode(const Node<DIM>& rNode);
    
    c_matrix<double,DIM,DIM> ComputeDiffusionTerm(const ChastePoint<DIM>& );

};

template<unsigned DIM>
CellwiseNutrientSinkPde<DIM>::CellwiseNutrientSinkPde(MeshBasedTissue<DIM>& rTissue, double coefficient)
        : mrTissue(rTissue),
          mCoefficient(coefficient)
{
}

template<unsigned DIM>
double CellwiseNutrientSinkPde<DIM>::ComputeConstantInUSourceTerm(const ChastePoint<DIM>& x)
{
    return 0.0;
}

template<unsigned DIM>
double CellwiseNutrientSinkPde<DIM>::ComputeLinearInUCoeffInSourceTerm(const ChastePoint<DIM>& x, Element<DIM,DIM>*)
{
    NEVER_REACHED;
    return 0.0;
}

template<unsigned DIM>
double CellwiseNutrientSinkPde<DIM>::ComputeLinearInUCoeffInSourceTermAtNode(const Node<DIM>& rNode)
{
    TissueCell& r_cell = mrTissue.rGetCellAtNodeIndex(rNode.GetIndex());
    if(r_cell.GetCellType()!=NECROTIC)
    {
        return -mCoefficient;
    }
    else
    {
        return 0.0;
    }
}

template<unsigned DIM>
c_matrix<double,DIM,DIM> CellwiseNutrientSinkPde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& )
{
    return identity_matrix<double>(DIM);
}  

#endif /*CELLWISENUTRIENTSINKPDE_HPP_*/
