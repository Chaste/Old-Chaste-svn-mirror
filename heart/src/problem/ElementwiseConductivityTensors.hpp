#ifndef ELEMENTWISECONDUCTIVITYTENSORS_HPP_
#define ELEMENTWISECONDUCTIVITYTENSORS_HPP_

#include <vector>
#include "UblasCustomFunctions.hpp"

#include <fstream>
#include <iostream>

/**
 * 
 *  This class provides an abstraction for the definition of constant non-constant difussion tensors
 * associated to the different elements of the mesh.
 * 
 *  After instanciating the class any of SetFibreOrientationFile() or SetNonConstantConductivities()
 * (or both) can be called to implement fiber orientation or heterogeneous conductivity into the
 * tensors, respectively. If none of them is called a constant tensor (with constant conductivies)
 * will be generated for all the elements of the mesh.
 * 
 *  Init() should be called to actually create the tensors.
 * 
 *  Initial values for conductivity from "Laminar Arrangement of Ventricular Myocytes Influences Electrical 
 * Behavior of the Heart", Hooks et al. 2007
 *     
 */

class ElementwiseConductivityTensors
{
private:
    unsigned mNumElements;

    std::string mFibreOrientationFilename;

    bool mUseNonConstantConductivities;
    bool mUseFibreOrientation;

    // Constant Conductivities
    double mConstLongConductivity; // mS/cm
    double mConstTransConductivity; // mS/cm
    double mConstNormalConductivity; // mS/cm

    // Non-constant Conductivities
    std::vector<double> mLongitudinalConductivities; // mS/cm
    std::vector<double> mTransverseConductivities; // mS/cm
    std::vector<double> mNormalConductivities; // mS/cm

    // Container for Tensors (single or multiple)
    std::vector< c_matrix<double,3,3> > mTensors;
    
    bool mInitialised;
    
public:

    ElementwiseConductivityTensors(unsigned numElements)
        : mNumElements(numElements),          
          mUseNonConstantConductivities(false),
          mUseFibreOrientation(false),
          mConstLongConductivity(7.0),
          mConstTransConductivity(3.5),
          mConstNormalConductivity(1.75),
          mInitialised(false)           
    {      
    }
    /**
     *  Sets a file for reading the fibre orientation from.
     * 
     *  @param rFibreOrientationFilename Relative path to the file defining the fibre orientation
     */
    void SetFibreOrientationFile(const std::string &rFibreOrientationFilename);

    /**
     *  Sets a constant longitudinal and transverse conductivity for all the elements of the mesh.
     * 
     *  @param constLongConduc Longitudinal conductivity (x axis)
     *  @param constTransConduc Transverse conductivity (y axis)
     *  @param constNormalConduc Normal conductivity (z axis)
     */   
    void SetConstantConductivities(double constLongConduc, double constTransConduc, double constNormalConduc);

    /**
     *  Sets a different longitudinal and transverse conductivity for every elements of the mesh.
     * 
     *  @param rLongitudinalConductivities Vector containing longitudinal conductivities of the elements (x axis)
     *  @param rTransverseConductivities Vector containing transverse conductivities of the elements (y axis)
     *  @param rNormalConductivities Vector containing normal conductivities of the elements (z axis)
     */      
    void SetNonConstantConductivities(
    			std::vector<double> &rLongitudinalConductivities, 
    			std::vector<double> &rTransverseConductivities,
    			std::vector<double> &rNormalConductivities);
    
    /**
     *  Computes the tensors (see comments in .cpp for more details).
     */
    void Init();
    
    /**
     *  Returns the diffussion tensor of the element number "index"
     * 
     *  @param index Index of the element of the mesh   
     */    
    c_matrix<double,3,3>& operator[](const unsigned index);
};

#endif /*ELEMENTWISECONDUCTIVITYTENSORS_HPP_*/
