#ifndef ELEMENTWISECONDUCTIVITYTENSORS_HPP_
#define ELEMENTWISECONDUCTIVITYTENSORS_HPP_

#include <vector>
#include "UblasCustomFunctions.hpp"

#include <fstream>
#include <iostream>

/**
 * 
 *  This class provides an abstraction for the definition of constant/non-constant difussion tensors
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

template<unsigned SPACE_DIM>
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
    std::vector<double>* mpLongitudinalConductivities; // mS/cm
    std::vector<double>* mpTransverseConductivities; // mS/cm
    std::vector<double>* mpNormalConductivities; // mS/cm

    // Container for Tensors (single or multiple)
    std::vector< c_matrix<double,SPACE_DIM,SPACE_DIM> > mTensors;
    
    bool mInitialised;
    
public:

    ElementwiseConductivityTensors()
        : mNumElements(1),          
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
    void SetFibreOrientationFile(const std::string &rFibreOrientationFilename)
    {
        mUseFibreOrientation = true;
        mFibreOrientationFilename = rFibreOrientationFilename;       
    }

    /**
     *  Sets a constant longitudinal and transverse conductivity for all the elements of the mesh.
     * 
     *  @param constLongConduc Longitudinal conductivity (x axis)
     *  @param constTransConduc Transverse conductivity (y axis)
     *  @param constNormalConduc Normal conductivity (z axis)
     */   
    void SetConstantConductivities(double constLongConduc, double constTransConduc=-DBL_MAX, double constNormalConduc=-DBL_MAX)
    {
        mUseNonConstantConductivities=false;
        
        switch (SPACE_DIM){
            case 3:
                /*
                 *  Miguel: assert or throw an exception?
                 */
                assert(constNormalConduc != -DBL_MAX);
                mConstNormalConductivity = constNormalConduc;
            case 2:
                assert(constTransConduc != -DBL_MAX);
                mConstTransConductivity = constTransConduc;
            case 1:            
                mConstLongConductivity = constLongConduc;
        }
    }

    /**
     *  Sets a different longitudinal and transverse conductivity for every elements of the mesh.
     * 
     *  @param rLongitudinalConductivities Vector containing longitudinal conductivities of the elements (x axis)
     *  @param rTransverseConductivities Vector containing transverse conductivities of the elements (y axis)
     *  @param rNormalConductivities Vector containing normal conductivities of the elements (z axis)
     */      
    void SetNonConstantConductivities(
    			std::vector<double> *pLongitudinalConductivities, 
    			std::vector<double> *pTransverseConductivities=NULL,
    			std::vector<double> *pNormalConductivities=NULL)
    {
        mUseNonConstantConductivities=true;
        
        switch (SPACE_DIM){
            case 3:
                assert(pNormalConductivities != NULL);        
                mpNormalConductivities=pNormalConductivities;
            case 2:
                assert(pTransverseConductivities != NULL);
                mpTransverseConductivities=pTransverseConductivities;
            case 1:            
                mpLongitudinalConductivities=pLongitudinalConductivities;
        }
    }                
    
    /**
     *  Computes the tensors based in all the info set
     */
    void Init()
    {
        std::ifstream data_file;
    
        c_matrix<double, SPACE_DIM, SPACE_DIM> conductivity_matrix(zero_matrix<double>(SPACE_DIM,SPACE_DIM));
        
        switch (SPACE_DIM){
            case 3:
                conductivity_matrix(2,2) = mConstNormalConductivity;
            case 2: 
                conductivity_matrix(1,1) = mConstTransConductivity;
            case 1:
                conductivity_matrix(0,0) = mConstLongConductivity;
        }
          
        
        if (!mUseNonConstantConductivities && !mUseFibreOrientation)
        {
            // Constant tensor for every element
            mTensors.push_back(conductivity_matrix);
        }
        else
        {
            //c_matrix<double,SPACE_DIM,SPACE_DIM> orientation_matrix(identity_matrix<double>(SPACE_DIM));
            c_matrix<double,SPACE_DIM,SPACE_DIM> orientation_matrix;
            orientation_matrix = identity_matrix<double>(SPACE_DIM);
                        
            if (mUseFibreOrientation)
            {
                data_file.open(mFibreOrientationFilename.c_str());
                assert(data_file.is_open());    
                data_file >> mNumElements;    
            }
            else
            {
                mNumElements = mpLongitudinalConductivities->size();
            }
                
            // reserve() allocates all the memory at once, more efficient than relying 
            // on the automatic reallocation scheme.
            mTensors.reserve(mNumElements);
                           
            for (unsigned element_index=0; element_index<mNumElements; element_index++)
                {
                /*
                 *  For every element of the mesh we compute its tensor like (from 
                 * "Laminar Arrangement of VentricularMyocites Influences Electrical 
                 * Behavior of the Heart", Darren et al. 2007):
                 * 
                 *                         [g_f  0   0 ] [a_f']
                 *  tensor = [a_f a_l a_n] [ 0  g_l  0 ] [a_l']
                 *                         [ 0   0  g_n] [a_n']
                 * 
                 *              [x_i]
                 *  where a_i = [y_i], i={f,l,n} are read from the fibre orientation file and
                 *              [z_i]
                 * 
                 *  g_f = longitudinal conductivity (constant or element specific)
                 *  g_l = transverse conductivity (constant or element specific)
                 *  g_n = normal conductivity (constant or element specific) 
                 * 
                 */
                
                if(mUseNonConstantConductivities)
                {
                    switch (SPACE_DIM){
                        case 3:
                            conductivity_matrix(2,2) = (*mpNormalConductivities)[element_index];
                        case 2:
                            conductivity_matrix(1,1) = (*mpTransverseConductivities)[element_index];
                        case 1:    
                            conductivity_matrix(0,0) = (*mpLongitudinalConductivities)[element_index];
                    }                                       
                }      
                
                if (mUseFibreOrientation)
                {
                    for (unsigned vector_index=0; vector_index<SPACE_DIM; vector_index++)
                    {
                        for (unsigned component_index=0; component_index<SPACE_DIM; component_index++)
                        {
                            data_file >> orientation_matrix(component_index, vector_index);
                        }                        
                    }
                }                    
                
                c_matrix<double,SPACE_DIM,SPACE_DIM> temp;
                noalias(temp) = prod(orientation_matrix, conductivity_matrix);  
                mTensors.push_back( prod(temp, trans(orientation_matrix) ) );   
            }
            
            if(mUseNonConstantConductivities)
            {
                data_file.close();
            }
        }
        
        mInitialised = true;
    }
    
    
    /**
     *  Returns the diffussion tensor of the element number "index"
     * 
     *  @param index Index of the element of the mesh   
     */    
    c_matrix<double,SPACE_DIM,SPACE_DIM>& operator[](const unsigned index)
    {
        assert(mInitialised);    
        
        if (!mUseNonConstantConductivities && !mUseFibreOrientation)
        {
            return mTensors[0];    
        }
        else
        {
            assert(index < mNumElements);
            return mTensors[index];
        }
    }
    
};

#endif /*ELEMENTWISECONDUCTIVITYTENSORS_HPP_*/
