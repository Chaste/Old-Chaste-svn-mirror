#ifndef CELLMUTATIONSTATES_HPP_
#define CELLMUTATIONSTATES_HPP_

/**
 * Possible types mutation state for TissueCells.
 */
typedef enum CellMutationState_
{
    HEALTHY,				// Wild-type cell
    APC_ONE_HIT,			// APC +/-
    APC_TWO_HIT,			// APC -/-
    BETA_CATENIN_ONE_HIT,	// Beta-catenin with a change at residue 45
    LABELLED,               // To paint a different colour but not actually mutant
    ALARCON_NORMAL,         // for use in Alarcon2004OxygenBasedCellCycleOdeSystem
    ALARCON_CANCER         // for use in Alarcon2004OxygenBasedCellCycleOdeSystem
} CellMutationState;


#endif /*CELLMUTATIONSTATES_HPP_*/
