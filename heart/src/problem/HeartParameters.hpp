#ifndef HEARTPARAMETERS_HPP_
#define HEARTPARAMETERS_HPP_

#include <memory>
#include "HeartPhysiologicalParameters.hpp"

class HeartParameters // todo: change to HeartConfig?
{
public:
    /**
     * Call this method to access the global parameters holder.
     * 
     * @return a single instance of the class
     */
    static HeartParameters* Instance();
    HeartPhysiologicalParametersType* Parameters();
    void SetParametersFile(std::string fileName);
    void Destroy();
    
private:
    HeartParameters();
    /** The single instance of the class */
    static HeartParameters* mpInstance;
    HeartPhysiologicalParametersType* mpParameters;
};

    
#endif /*HEARTPARAMETERS_HPP_*/

