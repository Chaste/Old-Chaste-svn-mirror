/**
 * AbstractVisualizer class. It writes data into correct file format.
 * i.e. For visualization with Matlab we need to modify the .node file
 * to change all "#" to "%". 
 */
#ifndef _ABSTRACTVISUALIZER_HPP_
#define _ABSTRACTVISUALIZER_HPP_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Exception.hpp"

class AbstractVisualizer 
{
protected:
	std::string mPathBaseName; /**<path base name of the files */
	int mDimension; /**< the dimension of nodes. */
	std::vector<double> mTimeSeries; /**< a vector to store the time steps which may be used as part of the file names. */
	bool mHasTimeFile; /**< a flag to indicate whether there is .time file, true if there is. */
	
	
public:
	AbstractVisualizer();
	virtual ~AbstractVisualizer();	
	virtual void CreateNodesFileForVisualization();/**<create .coord file */	     
	virtual void CreateOutputFileForVisualization();/**<create .val file which contains the output */
};

#endif //_ABSTRACTVISUALIZER_HPP_
