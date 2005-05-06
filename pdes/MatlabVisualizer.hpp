#ifndef _MATLABVISUALIZER_HPP_
#define _MATLABVISUALIZER_HPP_

/** 
 * Concrete version of the AbstractVisual class.
 * A MatlabVisualizer takes the base name and modify a little of 
 * the .node file and .out file so that all the comments are in 
 * Matlab style.
 */



#include <string>
#include <vector>

#include "AbstractVisualizer.hpp"

template<int SPACE_DIM>
class MatlabVisualizer: public AbstractVisualizer<SPACE_DIM>
{
private:
	std::string mPathBaseName; /**<path base name of the files */
	std::vector<double> mTimeSeries; /**< a vector to store the time steps which may be used as part of the file names. */
	bool mHasTimeFile; /**< a flag to indicate whether there is .time file, true if there is. */
	void CreateNodesFileForVisualization();	     
	void CreateOutputFileForVisualization();
	std::vector<std::string> GetRawDataFromFile(std::string fileName);
	
public:
	MatlabVisualizer(std::string pathBaseName);//, int dimension);
	~MatlabVisualizer();
	
	void CreateFilesForVisualization();
};

#endif //_MATLABVISUALIZER_HPP_
