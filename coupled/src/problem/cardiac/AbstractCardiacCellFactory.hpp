#ifndef ABSTRACTCARDIACCELLFACTORY_HPP_
#define ABSTRACTCARDIACCELLFACTORY_HPP_

#include "AbstractCardiacCell.hpp"

/**
 * Class which returns cardiac cells (which may have been stimulated).
 * For use with MonodomainPde.
 * Needed because when running in parallel, the user should specify only
 * which cells are used in a mesh, not on which process they exist.
 * A concrete implementation will be required, possibly for each simulation.
 */
class AbstractCardiacCellFactory
{
    public:
    virtual AbstractCardiacCell* CreateCardiacCellForNode(int)=0;
    virtual int GetNumberOfNodes()=0;
    virtual ~AbstractCardiacCellFactory() {}
};

#endif /*ABSTRACTCARDIACCELLFACTORY_HPP_*/

