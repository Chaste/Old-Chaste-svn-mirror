/**
 * AbstractVisualizer class. It writes data into correct file format.
 * i.e. For visualization with Matlab we need to modify the .node file
 * to change all "#" to "%". 
 * 
 * The constructor should contain the path-base name, and the dimension of the nodes.
 * The required existing files are 
 * 1. path_base_name.node
 * 2. path_base_name.out (if the output value doesn't depend on time)
 * 
 * OR
 * 1. path_base_name.node
 * 2. path_base_name.time (containing one column of time steps)
 * 3. path_base_name.xx.out (where xx stands for time steps. eg.path_base_name.0.out, path_base_name.0.1.out)
 * 			These are a series of files. 
 * Exceptions will be thrown if anyone is missing. 
 * 
 * The output files are 
 * 1.path_base_name.coord (containing coordinates of each node. It's a 1, 2 or 3 columns file)
 * 2.path_base_name.val (containing the result of the calculation. eg. the voltage at each node.) 
 * 			The .val file contains a matrix with each row representing the value at each time step.
 * 
 * Both 
 * 		CreateNodesFileForVisualization() and 
 * 		CreateOutputFileForVisualization()
 * need to be invoked in order to generate necessary files for Matlab visualization.
 * 
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
