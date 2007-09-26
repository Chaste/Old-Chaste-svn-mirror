#ifndef CRYPTVORONOIDATAWRITER_HPP_
#define CRYPTVORONOIDATAWRITER_HPP_

#include "Crypt.cpp"
#include "OutputFileHandler.hpp"

/**
 *  CryptVoronoiDataWriter
 *  
 *  Write position and voronoi cell area and perimeter data for each cell.
 *  Usage: create a crypt, create this writer. Everytime timestep create the 
 *  tessellation on the crypt and call WriteData() on this class. 
 * 
 *  Alternatively call SetLogged() on a cell and pass true for the final
 *  argument of the constructor so just log ONE cell. Note only one cell is 
 *  logged even if SetLogged() is set on many cells.
 * 
 *  Output form: each line in the file looks like
 *  
 *  current_time node_index1 x1 y1 A1 P1 node_index2 x2 y2 A2 P2 ....
 *   
 *  where (x,y) is the position, A the area of the Voronoi cell and P the perimeter
 * 
 *  Note the output directory is not cleaned.
 * 
 *  Templated over dim but only meant for dim=2
 */
template<unsigned DIM>
class CryptVoronoiDataWriter
{
private:
    Crypt<DIM>& mrCrypt;
    out_stream mOutStream;
    bool mFollowLoggedCell;
    
    
public:
    CryptVoronoiDataWriter(Crypt<DIM>& rCrypt, std::string directory, std::string filename, bool followLoggedCell = false)
        :mrCrypt(rCrypt),
         mFollowLoggedCell(followLoggedCell)
    {
        assert(DIM==2);
        OutputFileHandler output_file_handler(directory,false);
        mOutStream = output_file_handler.OpenOutputFile(filename);
        
    }

    ~CryptVoronoiDataWriter()
    {
        mOutStream->close();
    }
    
    void WriteData()
    {
        (*mOutStream)<< SimulationTime::Instance()->GetDimensionalisedTime() << " ";
        for (typename Crypt<DIM>::Iterator cell_iter = mrCrypt.Begin();
             cell_iter != mrCrypt.End();
             ++cell_iter)
        {
            if((!mFollowLoggedCell) || ((mFollowLoggedCell) && (cell_iter->IsLogged())))
            {
                unsigned node_index = cell_iter.GetNode()->GetIndex();
                double x = cell_iter.rGetLocation()[0];
                double y = cell_iter.rGetLocation()[1];
            
                double cell_area = mrCrypt.rGetVoronoiTessellation().GetFace(node_index)->GetArea();
                double cell_perimeter = mrCrypt.rGetVoronoiTessellation().GetFace(node_index)->GetPerimeter();
            
                (*mOutStream)<< node_index << " " << x << " " << y << " " << cell_area << " " << cell_perimeter << " ";
                
                if(mFollowLoggedCell)
                {
                    break;
                }
            }
        }
        (*mOutStream)<< "\n";
    }
};

#endif /*CRYPTVORONOIDATAWRITER_HPP_*/
