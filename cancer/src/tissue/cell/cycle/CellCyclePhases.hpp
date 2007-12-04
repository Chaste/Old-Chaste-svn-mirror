#ifndef CELLCYCLEPHASES_HPP_
#define CELLCYCLEPHASES_HPP_

/**
 * Possible phases of the cell cycle.
 * 
 * When our cells 'divide' they are actually entering M phase,
 * so a cell progresses round the cell cycle in the following sequence from its birth time
 * Divide-> M -> G0/G1 -> S -> G2 -> Divide.
 * 
 * G0 is a cell which stays in the G1 phase and is not going to divide. (i.e. quiescent or differentiated.)
 */
typedef enum CellCyclePhase_
{
    G_ZERO_PHASE,
    G_ONE_PHASE,
    S_PHASE,
    G_TWO_PHASE,
    M_PHASE    
} CellCyclePhase;


#endif /*CELLCYCLEPHASES_HPP_*/
