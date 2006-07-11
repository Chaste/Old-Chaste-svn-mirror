#ifndef CARDIACPROBLEMINPUTFILEREADER_HPP_
#define CARDIACPROBLEMINPUTFILEREADER_HPP_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "InputFileReader.hpp"

template<int SPACE_DIM>
class CardiacProblemInputFileReader : public InputFileReader
{
private :
    bool mIsMonodomainProblem;
    c_matrix<double, SPACE_DIM, SPACE_DIM> mIntracellularConductivity;
    c_matrix<double, SPACE_DIM, SPACE_DIM> mExtracellularConductivity;
    double mSurfaceAreaToVolumeRatio;
    double mCapacitance;
    
    double mPdeTimeStep;
    double mOdeTimeStep;
    double mPrintingTimeStep;
    double mEndTime;

    std::string mMeshFilename;
    std::string mOutputDirectory;
    std::string mOutputFilenamePrefix;   
    
    
public :
    CardiacProblemInputFileReader(std::string fileName) : InputFileReader(fileName)
    {
        bool found;
        std::string isMonodomain = ReadString("CardiacModel", found);
        if(isMonodomain=="Monodomain")
        {
            mIsMonodomainProblem = true;
        }
        else if(isMonodomain=="Bidomain")
        {
            mIsMonodomainProblem = false;
        }
        else
        {
            throw Exception("CardiacProblemInputFileReader.hpp: error, expected either 'Monodomain' or 'Bidomain' after 'CardiacModel:'");
        }

        unsigned number_of_vars = SPACE_DIM*(SPACE_DIM+1)/2; //ie 6 in 3d, 3 in 2d, 1 in 1d
        std::vector<double> sigma_i_vec = ReadVector<double>("IntracellularConductivity", 
                                                             number_of_vars, 
                                                             found);
        mIntracellularConductivity = VecToMatrix(sigma_i_vec);

        if(!mIsMonodomainProblem)
        {
            std::vector<double> sigma_e_vec = ReadVector<double>("ExtracellularConductivity", 
                                                                 number_of_vars, 
                                                                 found);
            mExtracellularConductivity = VecToMatrix(sigma_e_vec);
        }

        mSurfaceAreaToVolumeRatio = ReadDouble("SurfaceAreaToVolumeRatio",found);
        mCapacitance              = ReadDouble("Capacitance",found);
    
        mPdeTimeStep              = ReadDouble("PdeTimeStep",found);
        mOdeTimeStep              = ReadDouble("OdeTimeStep",found);
        mPrintingTimeStep         = ReadDouble("PrintingTimeStep",found);
        mEndTime                  = ReadDouble("EndTime",found);

        mMeshFilename             = ReadString("MeshFilename",found);
        mOutputDirectory          = ReadString("OutputDirectory",found);
        mOutputFilenamePrefix     = ReadString("OutputFilenamePrefix",found);
    }
    
    
    // Get methods
    bool IsMonodomainProblem()            {return mIsMonodomainProblem; }
  
    c_matrix<double, SPACE_DIM, SPACE_DIM> GetIntracellularConductivity() { return mIntracellularConductivity; }
    c_matrix<double, SPACE_DIM, SPACE_DIM> GetExtracellularConductivity() { return mExtracellularConductivity; }
 
    double GetSurfaceAreaToVolumeRatio()  {return mSurfaceAreaToVolumeRatio; } 
    double GetCapacitance()               {return mCapacitance; }
    
    double GetPdeTimeStep()               {return mPdeTimeStep; }
    double GetOdeTimeStep()               {return mOdeTimeStep; }
    double GetPrintingTimeStep()          {return mPrintingTimeStep; }
    double GetEndTime()                   {return mEndTime; }

    std::string GetMeshFilename()         {return mMeshFilename; }
    std::string GetOutputDirectory()      {return mOutputDirectory; }
    std::string GetOutputFilenamePrefix() {return mOutputFilenamePrefix; }
    

    
private :
    // helper function
    c_matrix<double, SPACE_DIM, SPACE_DIM> VecToMatrix(std::vector<double> vec)
    {
        assert(vec.size()==SPACE_DIM*(SPACE_DIM+1)/2);
        c_matrix<double, SPACE_DIM, SPACE_DIM> ret;
        
        if(SPACE_DIM==1)
        {
            ret(0,0) = vec[0];
        }
        else if(SPACE_DIM==2)
        {
            ret(0,0) = vec[0];
            ret(0,1) = vec[1];
            ret(1,0) = vec[1];
            ret(1,1) = vec[2];
        }
        else if(SPACE_DIM==3)
        {
            ret(0,0) = vec[0];
            ret(0,1) = vec[1];
            ret(1,0) = vec[1];
            ret(0,2) = vec[2];
            ret(2,0) = vec[2];
            ret(1,1) = vec[3];
            ret(1,2) = vec[4];
            ret(2,1) = vec[4];
            ret(2,2) = vec[5];
        }
        else
        {
            throw Exception("CardiacProblemInputFileReader.hpp: error, dimension should be 1, 2 or 3!");
        }
        return ret;                
    }
};

#endif /*CARDIACPROBLEMINPUTFILEREADER_HPP_*/
