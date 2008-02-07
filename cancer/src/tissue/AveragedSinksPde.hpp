#ifndef AVERAGEDSINKSPDE_HPP_
#define AVERAGEDSINKSPDE_HPP_

#include "MeshBasedTissue.hpp"
#include "AbstractLinearEllipticPde.hpp"



/** A pde which calculates the source term by adding the number of cells
 *  in the element containing that point and scaling by the element area
 */
template<unsigned DIM>
class AveragedSinksPde : public AbstractLinearEllipticPde<DIM>
{
private:
    MeshBasedTissue<DIM>& mrTissue;
    double mCoefficient;
    std::vector<double> mCellDensityOnCoarseElements;

public:
    AveragedSinksPde(MeshBasedTissue<DIM>& rTissue, double coefficient)
        : mrTissue(rTissue),
          mCoefficient(coefficient)
    {
    }

    void SetupSourceTerms(ConformingTetrahedralMesh<DIM,DIM>& rCoarseMesh) // must be called before solve
    {
        // allocate memory
        mCellDensityOnCoarseElements.resize(rCoarseMesh.GetNumElements());
        for(unsigned elem_index=0; elem_index<mCellDensityOnCoarseElements.size(); elem_index++)
        {
            mCellDensityOnCoarseElements[elem_index]=0.0;
        } 
        //loop over cells, find which coarse element it is in, and add 1 to the mSourceTermOnCoarseElements[elem_index];
        for(typename MeshBasedTissue<DIM>::Iterator cell_iter = mrTissue.Begin();
            cell_iter != mrTissue.End();
            ++cell_iter)
        {
            const ChastePoint<DIM>& r_position_of_cell = cell_iter.rGetLocation();
            unsigned elem_index = rCoarseMesh.GetContainingElementIndex(r_position_of_cell);
            mCellDensityOnCoarseElements[elem_index] += 1.0;
        }    
        
        // then divide each entry of mSourceTermOnCoarseElements by the element's area
        for(unsigned elem_index=0; elem_index<mCellDensityOnCoarseElements.size(); elem_index++)
        {
            mCellDensityOnCoarseElements[elem_index]/= rCoarseMesh.GetElement(elem_index)->GetVolume();
        }
    
    }

    double ComputeConstantInUSourceTerm(const ChastePoint<DIM>& x)
    {
        return 0.0;
    }
    
    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<DIM>& x, Element<DIM,DIM>* pElement) // now takes in element
    {
        assert(mCellDensityOnCoarseElements.size()>0);
        return mCoefficient*mCellDensityOnCoarseElements[pElement->GetIndex()];
    }
   
    c_matrix<double,DIM,DIM> ComputeDiffusionTerm(const ChastePoint<DIM>& )
    {
        return identity_matrix<double>(DIM);
    }   
};



#endif /*AVERAGEDSINKPDE_HPP_*/
