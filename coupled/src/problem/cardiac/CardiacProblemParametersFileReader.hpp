#ifndef CARDIACPROBLEMPARAMETERSFILEREADER_HPP_
#define CARDIACPROBLEMPARAMETERSFILEREADER_HPP_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "ParametersFileReader.hpp"
#include "GroupOfNumbersFileReader.hpp"

template<int SPACE_DIM>
class CardiacProblemParametersFileReader : public ParametersFileReader
{
private :
    bool mIsMonodomainProblem;

    c_matrix<double, SPACE_DIM, SPACE_DIM> mIntracellularConductivity;
    c_matrix<double, SPACE_DIM, SPACE_DIM> mExtracellularConductivity;
    double mSurfaceAreaToVolumeRatio;
    double mCapacitance;
    
    bool mIntracellularConductivityWasSet;
    bool mExtracellularConductivityWasSet;
    bool mSurfaceAreaToVolumeRatioWasSet;
    bool mCapacitanceWasSet;
    
    bool mThereAreFixedNodes;
    std::vector<unsigned> mFixedNodes;
    double mPdeTimeStep;
    double mOdeTimeStep;
    double mPrintingTimeStep;
    double mEndTime;

    std::string mMeshFilename;
    std::string mOutputDirectory;
    std::string mOutputFilenamePrefix;   
    
    std::vector<unsigned> mStimulatedNodes;
    double mStimulusStartTime;
    double mStimulusDuration;
    double mStimulusMagnitude;
    
    

public :
    CardiacProblemParametersFileReader(std::string fileName) : ParametersFileReader(fileName)
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
            EXCEPTION("Expected either 'Monodomain' or 'Bidomain' after 'CardiacModel:'");
        }

        unsigned number_of_vars = SPACE_DIM*(SPACE_DIM+1)/2; //ie 6 in 3d, 3 in 2d, 1 in 1d
        std::vector<double> sigma_i_vec = ReadVector<double>("IntracellularConductivity", 
                                                             number_of_vars, 
                                                             found,
                                                             false);  // don't quit if not found
        mIntracellularConductivityWasSet = found;
        if(found)
        {
            mIntracellularConductivity = VecToMatrix(sigma_i_vec);
        }


        if(!mIsMonodomainProblem)
        {
            std::vector<double> sigma_e_vec = ReadVector<double>("ExtracellularConductivity", 
                                                                 number_of_vars, 
                                                                 found,
                                                                 false);
            mExtracellularConductivityWasSet = found;
            if(found)
            {
                mExtracellularConductivity = VecToMatrix(sigma_e_vec);
            }
            
            std::string fixed_nodes_file = ReadString("FixedNodes",found);
            if(fixed_nodes_file!="NONE")
            {
                mThereAreFixedNodes = true;
                GroupOfNumbersFileReader<unsigned> numbers_reader(fixed_nodes_file);
                mFixedNodes = numbers_reader.GetData();
            }
            else
            {
                mThereAreFixedNodes = false;
            }
        }

        mSurfaceAreaToVolumeRatio = ReadDouble("SurfaceAreaToVolumeRatio",found,false);
        mSurfaceAreaToVolumeRatioWasSet = found;
        mCapacitance              = ReadDouble("Capacitance",found,false);
        mCapacitanceWasSet        = found;


        mPdeTimeStep              = ReadDouble("PdeTimeStep",found);
        mOdeTimeStep              = ReadDouble("OdeTimeStep",found);
        mPrintingTimeStep         = ReadDouble("PrintingTimeStep",found);
        mEndTime                  = ReadDouble("EndTime",found);

        mMeshFilename             = ReadString("MeshFilename",found);
        mOutputDirectory          = ReadString("OutputDirectory",found);
        mOutputFilenamePrefix     = ReadString("OutputFilenamePrefix",found);
        
        std::string stimulated_nodes_file = ReadString("StimulatedNodes",found);
        GroupOfNumbersFileReader<unsigned> group_of_numbers_reader(stimulated_nodes_file);
        mStimulatedNodes = group_of_numbers_reader.GetData();
            
        mStimulusStartTime        = ReadDouble("StimulusStartTime",found);   
        mStimulusDuration         = ReadDouble("StimulusDuration",found);   
        mStimulusMagnitude        = ReadDouble("StimulusMagnitude",found);   

        // eventually should have variables of the form
        //mOdeSolver                = ReadString("OdeSolver",found);
        //mCellModel                = ReadString("CellModel",found);
    }
    
    ///////////////////////////////////////////////////////////////////////////
    // Get methods
    ///////////////////////////////////////////////////////////////////////////
    bool IsMonodomainProblem()               { return mIsMonodomainProblem; }
  
    double GetPdeTimeStep()                  { return mPdeTimeStep; }
    double GetOdeTimeStep()                  { return mOdeTimeStep; }
    double GetPrintingTimeStep()             { return mPrintingTimeStep; }
    double GetEndTime()                      { return mEndTime; }

    std::string GetMeshFilename()            { return mMeshFilename; }
    std::string GetOutputDirectory()         { return mOutputDirectory; }
    std::string GetOutputFilenamePrefix()    { return mOutputFilenamePrefix; }
    
    bool ThereAreFixedNodes()                { return mThereAreFixedNodes; }
    std::vector<unsigned> GetFixedNodes()    { assert(mThereAreFixedNodes); return mFixedNodes; }
    
    std::vector<unsigned> GetStimulatedNodes() { return mStimulatedNodes; }
    double GetStimulusStartTime()            { return mStimulusStartTime; }
    double GetStimulusDuration()             { return mStimulusDuration; }
    double GetStimulusMagnitude()            { return mStimulusMagnitude; }
    
    bool IntracellularConductivityWasSet()   { return mIntracellularConductivityWasSet; }
    bool ExtracellularConductivityWasSet()   { return mExtracellularConductivityWasSet; }
    bool SurfaceAreaToVolumeRatioWasSet()    { return mSurfaceAreaToVolumeRatioWasSet; }
    bool CapacitanceWasSet()                 { return mCapacitanceWasSet; }
    
    c_matrix<double, SPACE_DIM, SPACE_DIM> GetIntracellularConductivity() 
    { 
        assert(mIntracellularConductivityWasSet);
        return mIntracellularConductivity; 
    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> GetExtracellularConductivity() 
    { 
        assert(mExtracellularConductivityWasSet);
        return mExtracellularConductivity; 
    }
 
    double GetSurfaceAreaToVolumeRatio()     
    { 
        assert(mSurfaceAreaToVolumeRatioWasSet); 
        return mSurfaceAreaToVolumeRatio; 
    } 

    double GetCapacitance()
    { 
        assert(mCapacitanceWasSet); 
        return mCapacitance; 
    }
    



    

private :
    ///////////////////////////////////////////////////////////////////////////
    // helper function - converts 6d vec into 3d symmetric matrix, etc
    ///////////////////////////////////////////////////////////////////////////
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
            EXCEPTION("Dimension should be 1, 2 or 3!");
        }
        return ret;                
    }
};

#endif /*CARDIACPROBLEMPARAMETERSFILEREADER_HPP_*/
