
#include "ChemotacticForce.hpp"
#include "CellwiseDataGradient.hpp"

template<unsigned DIM>
ChemotacticForce<DIM>::ChemotacticForce()
   : AbstractForce<DIM>()
{
}

template<unsigned DIM>
ChemotacticForce<DIM>::~ChemotacticForce()
{
}

template<unsigned DIM>
double ChemotacticForce<DIM>::GetChemotacticForceMagnitude(const double concentration, const double concentrationGradientMagnitude)
{
    return concentration; // temporary force law - can be changed to something realistic
                          // without tests failing
}

template<unsigned DIM>
void ChemotacticForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                    AbstractTissue<DIM>& rTissue)
{
    CellwiseDataGradient<DIM> gradients;
    gradients.SetupGradients();

    for (typename AbstractTissue<DIM>::Iterator cell_iter = rTissue.Begin();
         cell_iter != rTissue.End();
         ++cell_iter)
    {
        // Only LABELLED cells move chemotactically
        if (cell_iter->GetMutationState() == LABELLED)
        {
            TissueCell& cell = *cell_iter;
            unsigned node_global_index = cell.GetLocationIndex();

            c_vector<double,DIM>& r_gradient = gradients.rGetGradient(cell.GetLocationIndex());
            double nutrient_concentration = CellwiseData<DIM>::Instance()->GetValue(&cell,0);
            double magnitude_of_gradient = norm_2(r_gradient);

            double force_magnitude = GetChemotacticForceMagnitude(nutrient_concentration, magnitude_of_gradient);

            // force +=  chi * gradC/|gradC|
            if (magnitude_of_gradient > 0)
            {
                rForces[node_global_index] += (force_magnitude/magnitude_of_gradient)*r_gradient;
            }
            // else Fc=0
        }
    }
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class ChemotacticForce<1>;
template class ChemotacticForce<2>;
template class ChemotacticForce<3>;
