#include "ElementwiseConductivityTensors.hpp"

void ElementwiseConductivityTensors::SetFibreOrientationFile(const std::string &rFibreOrientationFilename)
{
    mUseFibreOrientation = true;
    mFibreOrientationFilename = rFibreOrientationFilename;       
}

void ElementwiseConductivityTensors::SetConstantConductivities(double constLongConduc, double constTransConduc)
{
    mUseNonConstantConductivities=false;
    mConstLongConductivity = constLongConduc;
    mConstTransConductivity = constTransConduc;
}

void ElementwiseConductivityTensors::SetNonConstantConductivities(std::vector<double> &rLongitudinalConductivities, std::vector<double> &rTransverseConductivities)
{
    assert(rLongitudinalConductivities.size()==mNumElements);
    assert(rTransverseConductivities.size()==mNumElements);
    
    mUseNonConstantConductivities=true;
    mLongitudinalConductivities=rLongitudinalConductivities;
    mTransverseConductivities=rTransverseConductivities;
}

void ElementwiseConductivityTensors::Init()
{
    std::ifstream data_file;
    
    if (!mUseNonConstantConductivities && !mUseFibreOrientation)
    {
        c_matrix<double,3,3> const_tensor(zero_matrix<double>(3,3));
        const_tensor(0,0) = mConstLongConductivity;
        const_tensor(1,1) = mConstTransConductivity;
        const_tensor(2,2) = mConstTransConductivity; 
        mTensors.push_back(const_tensor);
    }
    else
    {
        // reserve() allocates all the memory at once, more efficient than relying 
        // on the automatic reallocation scheme.
        mTensors.reserve(mNumElements);
        
        c_matrix<double, 3, 3> conductivity_matrix(zero_matrix<double>(3,3));
        c_matrix<double,3,3> orientation_matrix(identity_matrix<double>(3));
        
        if (!mUseNonConstantConductivities)
        {
            conductivity_matrix(0,0) = mConstLongConductivity;
            conductivity_matrix(1,1) = mConstTransConductivity;
            conductivity_matrix(2,2) = mConstTransConductivity;          
        }
    
        if (mUseFibreOrientation)
        {
            data_file.open(mFibreOrientationFilename.c_str());
            assert(data_file.is_open());
    
            unsigned elements_in_file;
            data_file >> elements_in_file;
    
            assert(elements_in_file==mNumElements);
        }
    
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
             *  g_l = g_n = transverse conductivity (constant or element specific)
             * 
             */
            
            if(mUseNonConstantConductivities)
            {
                conductivity_matrix(0,0) = mLongitudinalConductivities[element_index];
                conductivity_matrix(1,1) = mTransverseConductivities[element_index];
                conductivity_matrix(2,2) = mTransverseConductivities[element_index];
            }      
            
            if (mUseFibreOrientation)
            {
                // Read a_f                
                data_file >> orientation_matrix(0,0) >> orientation_matrix(1,0) >> orientation_matrix(2,0);                   
                // Read a_l
                data_file >> orientation_matrix(0,1) >> orientation_matrix(1,1) >> orientation_matrix(2,1);            
                // Read a_n
                data_file >> orientation_matrix(0,2) >> orientation_matrix(1,2) >> orientation_matrix(2,2);
            }                    
            
            c_matrix<double,3,3> temp;
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

c_matrix<double,3,3>& ElementwiseConductivityTensors::operator[](const unsigned index)
{
    assert(mInitialised);
    assert(index < mNumElements);
    
    if (!mUseNonConstantConductivities && !mUseFibreOrientation)
    {
        return mTensors[0];    
    }
    else
    {
        return mTensors[index];
    }
}
