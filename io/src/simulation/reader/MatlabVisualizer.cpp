#ifndef _MATLABVISUALIZER_CPP_
#define _MATLABVISUALIZER_CPP_

#include "MatlabVisualizer.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

#include "global/src/Exception.hpp"

/**
 * The constructor of MatlabVisualizer.
 * Read in the time step file which has an extension .time
 * There should be "<pathBaseName>.xx.out" files. <pathBaseName> is the same
 * as mPathBaseName, "xx" is the time (converted from double to string).
 */
template<int SPACE_DIM> 
MatlabVisualizer<SPACE_DIM>::MatlabVisualizer(std::string outputPathBaseName, std::string inputPathBaseName)
{
	
	mOutputPathBaseName = outputPathBaseName;
	mHasTimeFile = false;
	
	if (inputPathBaseName.empty())
	{
		mInputPathBaseName = outputPathBaseName;
	} else {
		mInputPathBaseName = inputPathBaseName;	
	}
	
	std::stringstream time_file_name_stream;
	time_file_name_stream<<mOutputPathBaseName<<"Time.dat";
	std::string time_file_name=time_file_name_stream.str();
	
	std::ifstream data_file(time_file_name.c_str());

	if (data_file.is_open())
	{
		mHasTimeFile = true;
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
MatlabVisualizer<SPACE_DIM>::~MatlabVisualizer()
{
}

template<int SPACE_DIM> 
void MatlabVisualizer<SPACE_DIM>::CreateFilesForVisualization()
{
	CreateNodesFileForVisualization();
	CreateOutputFileForVisualization();
}

/**
 * Create coordinates files with 1, 2 or 3 columns.
 */
template<int SPACE_DIM> 
void MatlabVisualizer<SPACE_DIM>::CreateNodesFileForVisualization()
{
	//Open node file and store the lines as a vector of strings (minus the comments) 	
	std::string node_file_name=mInputPathBaseName+".node";
	std::vector<std::string> node_data=GetRawDataFromFile(node_file_name);

	//Write node file for Matlab - without any comments
	std::string output_node_file_name=mOutputPathBaseName+".coord";
	std::ofstream output_node_file(output_node_file_name.c_str());
	
	for (int i = 1; i<node_data.size(); i++)
	{
		std::stringstream line_stream(node_data[i]);
		double value;
		line_stream>>value;
		for (int dim=0; dim<SPACE_DIM; dim++)
		{
			line_stream>>value;
			output_node_file << value<<" ";			
		}
		output_node_file <<"\n";

	}
	output_node_file.close();	
}

/**
 * Create a .val file containing a matrix(time*coord) of output value.
 */
//template<int SPACE_DIM>
//void MatlabVisualizer<SPACE_DIM>::CreateOutputFileForVisualization()
//{
//	std::vector< std::string > in_file_names;
//	//Write output file for Matlab - without any comments
//	std::string output_value_file_name=mPathBaseName+".val";
//	std::ofstream output_value_file(output_value_file_name.c_str());
//	int num_files = mTimeSeries.size();
//	if (num_files == 0 )
//	{
//		//only one output file, no time series
//		in_file_names.push_back(mPathBaseName + ".dat");
//	}else
//	{		
//		for (int i=0; i< num_files; i++)
//		{
//			std::stringstream one_file_name_stream;
//			one_file_name_stream<<mPathBaseName<<"_"<<i<<".dat";
//			in_file_names.push_back(one_file_name_stream.str());
//		}
//	}
//	int i=0;
//	std::cout<<"num_files="<<num_files<<"\n";
//	do
//	{		
//		//Open node file and store the lines as a vector of strings (minus the comments) 	
//		std::vector<std::string> output_data=GetRawDataFromFile(in_file_names[i]);
//		std::cout<<"output_data.size="<<output_data.size()<<"\n";
//		for (int j=1 ; j<output_data.size(); j++)
//		{
//			std::stringstream line_stream(output_data[j]);
//			double value;
//			line_stream>>value;
//			output_value_file << value<<" ";			
//		}
//		output_value_file <<"\n";
//		std::cout<<"num_files, stops at i="<<i<<"\n";
//	}while(i++<num_files);
//
//	output_value_file.close();	
//}

template<int SPACE_DIM>
void MatlabVisualizer<SPACE_DIM>::CreateOutputFileForVisualization()
{
	FILE *fw, *fr;
	char outfilename[1024], infilename[1024];
	sprintf(outfilename, "%s.val", mOutputPathBaseName.c_str());
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
			sprintf(infilename,"%s.dat", mOutputPathBaseName.c_str());
		}else
		{
			sprintf(infilename, "%s_%d.dat",mOutputPathBaseName.c_str(), i);
		}
		if ( (fr=fopen(infilename,"r"))==NULL)
		{
		    std::cout << infilename << std::endl;
			throw Exception("The file for reading cannot be created successfully");
		}
		char buffer[1024];
		fgets(buffer, 1024, fr);	//header

		while(!feof( fr ))
		{			
			double value;
			int num_of_fields = fscanf(fr, "%lf", &value);

			if (num_of_fields>0)
			{
				fprintf(fw, "%lf ", value);
			}
		}
		fprintf(fw, "\n");
		fclose(fr);
	
	} while(i++<num_files-1);
	fclose(fw);
	
//	std::vector< std::string > in_file_names;
//	Write output file for Matlab - without any comments
//	std::string output_value_file_name=mPathBaseName+".val";
//	std::ofstream output_value_file(output_value_file_name.c_str());
//	int num_files = mTimeSeries.size();
//	if (num_files == 0 )
//	{
//		only one output file, no time series
//		in_file_names.push_back(mPathBaseName + ".dat");
//	}else
//	{		
//		for (int i=0; i< num_files; i++)
//		{
//			std::stringstream one_file_name_stream;
//			one_file_name_stream<<mPathBaseName<<"_"<<i<<".dat";
//			in_file_names.push_back(one_file_name_stream.str());
//		}
//	}
//	int i=0;
//	std::cout<<"num_files="<<num_files<<"\n";
//	do
//	{		
//		Open node file and store the lines as a vector of strings (minus the comments) 	
//		std::vector<std::string> output_data=GetRawDataFromFile(in_file_names[i]);
//		std::cout<<"output_data.size="<<output_data.size()<<"\n";
//		for (int j=1 ; j<output_data.size(); j++)
//		{
//			std::stringstream line_stream(output_data[j]);
//			double value;
//			line_stream>>value;
//			output_value_file << value<<" ";			
//		}
//		output_value_file <<"\n";
//		std::cout<<"num_files, stops at i="<<i<<"\n";
//	}while(i++<num_files);
//
//	output_value_file.close();	
}

/**
 * Read the file indicated by fileName and return a vector of strings. 
 * It removes all comments and white lines.
 */
 template<int SPACE_DIM>
std::vector<std::string> MatlabVisualizer<SPACE_DIM>::GetRawDataFromFile(std::string fileName)
{
	// Open raw data file
	
	std::vector<std::string> RawDataFromFile;
	std::ifstream dataFile(fileName.c_str());
	
	
	
	/**
	 * Checks that input file has been opened correctly. If not throws an 
	 * exception that should be caught by the user.
	 * 
	 */
	if (!dataFile.is_open())	
	{
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

#endif //_MATLABVISUALIZER_CPP_
