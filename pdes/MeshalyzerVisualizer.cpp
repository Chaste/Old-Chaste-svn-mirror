#ifndef _MESHALYZERVISUALIZER_CPP_
#define _MESHALYZERVISUALIZER_CPP_
#include "MeshalyzerVisualizer.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

#include "Exception.hpp"

/**
 * The constructor of MESHALYZERVisualizer.
 * Read in the time step file which has an extension Time.dat
 * There should be "<pathBaseName>.xx.out" files. <pathBaseName> is the same
 * as mPathBaseName, "xx" is the time (converted from double to string).
 */
template<int SPACE_DIM>  
MeshalyzerVisualizer<SPACE_DIM>::MeshalyzerVisualizer(std::string PathBaseName)
{
	if (SPACE_DIM != 3)
	{
		throw("MeshalyzerVisualizer can only handle 3D data.");
	}
	mPathBaseName = PathBaseName;
		
	std::stringstream time_file_name_stream;
	time_file_name_stream<<mPathBaseName<<"Time.dat";
	std::string time_file_name=time_file_name_stream.str();
	
	std::ifstream data_file(time_file_name.c_str());

	if (data_file.is_open())
	{
	
			// Read each line in turn
	
		std::string RawLineFromFile;
		getline(data_file, RawLineFromFile);
		getline(data_file, RawLineFromFile);
		
		while(data_file)
		{
		
			// Remove comments
			
			int hashLocation=RawLineFromFile.find('#',0);
			if (hashLocation >= 0)
			{
				RawLineFromFile=RawLineFromFile.substr(0,hashLocation);
			}
			
		
			// Remove blank lines
			
			int notBlankLocation=RawLineFromFile.find_first_not_of(" \t",0);
			if (notBlankLocation >= 0)
			{
				double value;
				std::stringstream line_data(RawLineFromFile);
				line_data >> value;
				mTimeSeries.push_back(value);
			}		
			
			// Move onto next line
			
			getline(data_file, RawLineFromFile);
		}	
		data_file.close();	
	}	
}

template<int SPACE_DIM> 
MeshalyzerVisualizer<SPACE_DIM>::~MeshalyzerVisualizer()
{
}

template<int SPACE_DIM> 
void MeshalyzerVisualizer<SPACE_DIM>::CreateFilesForVisualization()
{
	CreateNodesFileForVisualization();
	CreateOutputFileForVisualization();
}

/**
 * Create coordinates files and surface files(.tris) with 1, 2 or 3 columns.
 */
template<int SPACE_DIM> 
void MeshalyzerVisualizer<SPACE_DIM>::CreateNodesFileForVisualization()
{
	AbstractMeshWriter *pMeshWriter = new MeshalyzerMeshWriter(
							mPathBaseName);
	AbstractMeshReader *pImportMeshReader=new TrianglesMeshReader(
							mPathBaseName);							
	int i;
	for (i=0; i<pImportMeshReader->GetNumNodes();i++)
	{
		pMeshWriter->SetNextNode(pImportMeshReader->GetNextNode());
	}
	for (i=0; i<pImportMeshReader->GetNumElements();i++)
	{
		pMeshWriter->SetNextElement(pImportMeshReader->GetNextElement());
	}
	for (i=0; i<pImportMeshReader->GetNumBoundaryFaces();i++)
	{
		pMeshWriter->SetNextBoundaryFace(pImportMeshReader->GetNextBoundaryFace());
	}
	
	pMeshWriter->WriteFiles();
	
	//The commented lines do the same thing as well
//	//Open node file and store the lines as a vector of strings (minus the comments) 	
//	std::string node_file_name=mPathBaseName+".node";
//	std::vector<std::string> node_data=GetRawDataFromFile(node_file_name);
//
//	//Write node file for Meshalyzer - without any comments
//	std::string output_node_file_name=mPathBaseName+".pts";
//	std::ofstream output_node_file(output_node_file_name.c_str());
//	if (output_node_file.is_open())
//	{
//		double number_points;
//		//header
//		std::stringstream line_stream(node_data[0]);
//		line_stream >> number_points;
//		output_node_file << number_points<<std::endl;
//		
//		for (int i = 1; i<node_data.size(); i++)
//		{
//			std::stringstream line_stream(node_data[i]);
//			double value;
//			line_stream>>value;	//index, ignored
//			for (int dim=0; dim<SPACE_DIM; dim++)
//			{
//				line_stream>>value;
//				output_node_file << value<<" ";			
//			}
//			output_node_file <<"\n";
//	
//		}
//		output_node_file.close();
//	}else
//	{
//		throw("In MeshalyzerVisualizer, .pts file cannot be created.");
//	}
//
//	//Open node file and store the lines as a vector of strings (minus the comments) 	
//	std::string face_file_name=mPathBaseName+".face";
//	std::vector<std::string> face_data=GetRawDataFromFile(face_file_name);
//	//Write tris. file for Meshalyzer - without any comments
//	std::string output_tris_file_name=mPathBaseName+".tris";
//	std::ofstream output_tris_file(output_tris_file_name.c_str());	
//	if (output_tris_file.is_open())
//	{
//		double number_points;
//		//header
//		std::stringstream line_stream(face_data[0]);
//		line_stream >> number_points;	
//		output_node_file << number_points<<std::endl;
//				
//		for (int i = 1; i<face_data.size(); i++)
//		{
//			std::stringstream line_stream(face_data[i]);
//			double value;
//			line_stream>>value;	//index, ignored
//			for (int dim=0; dim<SPACE_DIM; dim++)
//			{
//				line_stream>>value;
//				output_tris_file << value<<" ";			
//			}
//			output_tris_file<<"0";
//			output_tris_file <<"\n";
//	
//		}
//		output_tris_file.close();
//	}
//	else
//	{
//		throw("The .tris file cannot be created.");
//	}
//	
}


template<int SPACE_DIM>
void MeshalyzerVisualizer<SPACE_DIM>::CreateOutputFileForVisualization()
{
	FILE *fw, *fr;
	char outfilename[1024], infilename[1024];
	sprintf(outfilename, "%s.tdat", mPathBaseName.c_str());
	if ( (fw = fopen(outfilename,"w"))==NULL)
	{	
		throw Exception("The file for writing cannot be created successfully");
	}
	
	int num_files = mTimeSeries.size();
	int i=0;
	do
	{
		if (num_files==0)	//time independent
		{
			sprintf(infilename,"%s.dat", mPathBaseName.c_str());
		}else
		{
			sprintf(infilename, "%s_%d.dat",mPathBaseName.c_str(), i);
		}
		if ( (fr=fopen(infilename,"r"))==NULL)
		{		
			throw Exception("The file for reading cannot be created successfully");
		}
		char buffer[1024];
		fgets(buffer, 1024, fr);	//header in the input file

		while(!feof( fr ))
		{			
			
			double value;
			int num_of_fields = fscanf(fr, "%lf", &value);

			if (num_of_fields>0)
			{
				fprintf(fw, "%E \n", value);
			}
		}
	
	}while(i++<num_files-1);
	fclose(fr);
	fclose(fw);
	

}

/**
 * Read the file indicated by fileName and return a vector of strings. 
 * It removes all comments and white lines.
 */
 template<int SPACE_DIM>
std::vector<std::string> MeshalyzerVisualizer<SPACE_DIM>::GetRawDataFromFile(std::string fileName)
{
	// Open raw data file
	
	std::vector<std::string> RawDataFromFile;
	std::ifstream dataFile(fileName.c_str());
	
	
	
	/**
	 * Checks that input file has been opened correctly. If not throws an 
	 * exception that should be caught by the user.
	 * 
	 */
	std::cout<<fileName<<"\n";
	if (!dataFile.is_open())	
	{
		std::cout<<"in the throw body\n";
		throw Exception("Could not open data file "+fileName+" .");
	}
	
	
	
	// Read each line in turn
	
	std::string RawLineFromFile;
	getline(dataFile, RawLineFromFile);
	
	while(dataFile){
		
		// Remove comments
		
		int hashLocation=RawLineFromFile.find('#',0);
		if (hashLocation >= 0)
		{
			RawLineFromFile=RawLineFromFile.substr(0,hashLocation);
		}
		
	
		// Remove blank lines
		
		int notBlankLocation=RawLineFromFile.find_first_not_of(" \t",0);
		if (notBlankLocation >= 0)
		{
			RawDataFromFile.push_back(RawLineFromFile);	
		}
		
		
		// Move onto next line
		
		getline(dataFile, RawLineFromFile);
	}


	dataFile.close(); // Closes the data file
	
	return(RawDataFromFile);
}

#endif //_MESHALYZERVISUALIZER_CPP_
