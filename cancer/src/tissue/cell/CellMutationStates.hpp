#ifndef CELLMUTATIONSTATES_HPP_
#define CELLMUTATIONSTATES_HPP_

/**
 * Possible types mutation state for TissueCells.
 */
typedef enum CellMutationState_
{
    HEALTHY,                // Wild-type cell
    APC_ONE_HIT,            // APC +/-
    APC_TWO_HIT,            // APC -/-
    BETA_CATENIN_ONE_HIT,   // Beta-catenin with a change at residue 45
    LABELLED,               // To paint a different colour but not actually mutant
} CellMutationState;

const static unsigned NUM_CELL_MUTATION_STATES=5;

#endif /*CELLMUTATIONSTATES_HPP_*/
