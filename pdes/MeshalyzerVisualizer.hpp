#ifndef _MESHALYZERVISUALIZER_HPP_
#define _MESHALYZERVISUALIZER_HPP_

#include "AbstractVisualizer.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include "TrianglesMeshReader.hpp"
#include <string>
#include <vector>

template<int SPACE_DIM>
class MeshalyzerVisualizer: public AbstractVisualizer<SPACE_DIM>
{
private:	
	std::string mPathBaseName; /**<path base name of the files */
	std::vector<double> mTimeSeries; /**< a vector to store the time steps which may be used as part of the file names. */
	void CreateNodesFileForVisualization();	     
	void CreateOutputFileForVisualization();
	std::vector<std::string> GetRawDataFromFile(std::string fileName);
	
public:
	MeshalyzerVisualizer(std::string PathBaseName);
	virtual ~MeshalyzerVisualizer();
	void CreateFilesForVisualization();
};

#endif //_MESHALYZERVISUALIZER_HPP_
