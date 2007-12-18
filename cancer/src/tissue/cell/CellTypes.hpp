#ifndef CELLTYPES_HPP_
#define CELLTYPES_HPP_

/**
 * Possible types of TissueCell.
 */
typedef enum CellType_
{
    STEM,
    TRANSIT,
    DIFFERENTIATED,
    NECROTIC // for use in tissue simulations with nutrients
} CellType;


#endif /*CELLTYPES_HPP_*/
