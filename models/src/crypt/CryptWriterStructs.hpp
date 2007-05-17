#ifndef CRYPTWRITERSTRUCTS_HPP_
#define CRYPTWRITERSTRUCTS_HPP_

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "UblasCustomFunctions.hpp"

/**
 * Structure encapsulating variable identifiers for the node datawriter
 */
typedef struct node_writer_ids_t
{
    unsigned time;               /**< The simulation time */
    std::vector<unsigned> types; /**< Cell types */
    /** Cell positions */
    std::vector<c_vector<unsigned, 3> > position_id;
}
node_writer_ids_t;

/**
 * Structure encapsulating variable identifiers for the element datawriter
 */
typedef struct element_writer_ids_t
{
    unsigned time;/**< The simulation time */
    /** Node indices */
    std::vector<c_vector<unsigned, 4> > node_id;
}
element_writer_ids_t;

#endif /*CRYPTWRITERSTRUCTS_HPP_*/
