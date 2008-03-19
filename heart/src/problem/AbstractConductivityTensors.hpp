#ifndef ABSTRACTCONDUCTIVITYTENSORS_HPP_
#define ABSTRACTCONDUCTIVITYTENSORS_HPP_

template<unsigned SPACE_DIM>
class AbstractConductivityTensors
{
protected:
    unsigned mNumElements;

    bool mUseNonConstantConductivities;
    bool mUseFibreOrientation;

    // Constant Conductivities
    c_vector<double, SPACE_DIM> mConstantConductivities; // mS/cm

    // Non-constant Conductivities
    std::vector<c_vector<double, SPACE_DIM> >* mpNonConstantConductivities; // mS/cm

    // Container for Tensors (single or multiple)
    std::vector< c_matrix<double,SPACE_DIM,SPACE_DIM> > mTensors;
    
    bool mInitialised;

    std::string mFibreOrientationFilename;
    std::ifstream mDataFile;


    void OpenFibreOrientationFile()
    {
        mDataFile.open(mFibreOrientationFilename.c_str());
        if (!mDataFile.is_open())
        {
            EXCEPTION("Wrong fibre orientation file name.");    
        }    
    }
    
    void CloseFibreOrientationFile()
    {
        mDataFile.close();        
    }

    template<class TOKENS_TYPE> 
    unsigned GetTokensAtNextLine(std::vector<TOKENS_TYPE>& tokens)
    {
        std::string line;

        bool comment_line;
        bool blank_line;
         
        do
        {
            getline(mDataFile, line);                       

            if (mDataFile.eof())
            {
                CloseFibreOrientationFile();
                EXCEPTION("Fibre orientation file contains less fibre definitions than the number of elements in the mesh");    
            }

            comment_line = (line.find('#',0) != std::string::npos);
            blank_line = (line.find_first_not_of(" \t",0) == std::string::npos);
        }
        while(comment_line || blank_line);
         
 
        std::stringstream line_stream(line);
    
        // Read all the numbers from the line
        while (!line_stream.eof())
        {
            TOKENS_TYPE item;
            line_stream >> item;
            tokens.push_back(item);
        }        
        
        return tokens.size();                
    }

    unsigned GetNumElementsFromFile()
    {
        std::vector<unsigned> tokens;
        
        if (GetTokensAtNextLine<unsigned>(tokens) != 1)
        {
            CloseFibreOrientationFile();
            EXCEPTION("First (non comment) line of the fibre orientation file should contain the number of elements of the mesh (and nothing else)");    
        }

        return tokens[0];
    }
    
public:

    AbstractConductivityTensors()
        : mNumElements(1),          
          mUseNonConstantConductivities(false),
          mUseFibreOrientation(false),
          mInitialised(false)           
    {
        double init_data[]={7.0, 3.5, 1.75};
        
        for (unsigned dim=0; dim<SPACE_DIM; dim++)
        {
             mConstantConductivities[dim] = init_data[dim];   
        }      
    }
    
    virtual ~AbstractConductivityTensors()
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
     *  Sets constant conductivities for all the elements of the mesh.
     * 
     *  @param constLongConduc Longitudinal conductivity (x axis)
     *  @param constTransConduc Transverse conductivity (y axis)
     *  @param constNormalConduc Normal conductivity (z axis)
     * 
     *  Explain problem with c_vector and SPACE_DIM
     */   
    void SetConstantConductivities(c_vector<double, 1> constantConductivities)
    {       
        if (SPACE_DIM != 1)
        {
            EXCEPTION("Wrong number of conductivities provided");
        }
        
        mUseNonConstantConductivities = false;        
        mConstantConductivities = constantConductivities;        
    }
    
    void SetConstantConductivities(c_vector<double, 2> constantConductivities)
    {
        if (SPACE_DIM != 2)
        {
            EXCEPTION("Wrong number of conductivities provided");
        }
        
        mUseNonConstantConductivities = false;        
        mConstantConductivities = constantConductivities;        
    }

    void SetConstantConductivities(c_vector<double, 3> constantConductivities)
    {
        if (SPACE_DIM != 3)
        {
            EXCEPTION("Wrong number of conductivities provided");
        }
          
        mUseNonConstantConductivities = false;        
        mConstantConductivities = constantConductivities;        
    }
    

    /**
     *  Sets a different longitudinal and transverse conductivity for every elements of the mesh.
     * 
     *  @param rLongitudinalConductivities Vector containing longitudinal conductivities of the elements (x axis)
     *  @param rTransverseConductivities Vector containing transverse conductivities of the elements (y axis)
     *  @param rNormalConductivities Vector containing normal conductivities of the elements (z axis)
     */      
    void SetNonConstantConductivities(std::vector<c_vector<double, SPACE_DIM> >* pNonConstantConductivities)
    {
        mUseNonConstantConductivities = true;
        mpNonConstantConductivities = pNonConstantConductivities;
    }                
    
    /**
     *  Computes the tensors based in all the info set
     */
    virtual void Init() throw (Exception) = 0;
        
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

#endif /*ABSTRACTCONDUCTIVITYTENSORS_HPP_*/
