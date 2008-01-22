#ifndef TISSUEVORONOIDATAWRITER_HPP_
#define TISSUEVORONOIDATAWRITER_HPP_

#include "Tissue.cpp"
#include "OutputFileHandler.hpp"
#include "CellTypes.hpp"

/**
 *  TissueVoronoiDataWriter
 *  
 *  Write position and voronoi cell area and perimeter data for each cell.
 *  Usage: create a tissue, create this writer. Everytime timestep create the 
 *  tessellation on the tissue and call WriteData() on this class. 
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
class TissueVoronoiDataWriter
{
private:
    Tissue<DIM>& mrTissue;
    out_stream mOutStream;
    bool mFollowLoggedCell;
    
    
public:
    TissueVoronoiDataWriter(Tissue<DIM>& rTissue, std::string directory, std::string filename, bool followLoggedCell = false)
        :mrTissue(rTissue),
         mFollowLoggedCell(followLoggedCell)
    {
        #define COVERAGE_IGNORE
        assert(DIM==2);
        #undef COVERAGE_IGNORE
        OutputFileHandler output_file_handler(directory,false);
        mOutStream = output_file_handler.OpenOutputFile(filename);
        (*mOutStream)<< std::setprecision(8);        
    }

    ~TissueVoronoiDataWriter()
    {
        mOutStream->close();
    }
    
    void WriteData()
    {
        (*mOutStream)<< SimulationTime::Instance()->GetDimensionalisedTime() << " ";
        for (typename Tissue<DIM>::Iterator cell_iter = mrTissue.Begin();
             cell_iter != mrTissue.End();
             ++cell_iter)
        {
            if((!mFollowLoggedCell) || ((mFollowLoggedCell) && (cell_iter->IsLogged())))
            {
                unsigned node_index = cell_iter.GetNode()->GetIndex();
                double x = cell_iter.rGetLocation()[0];
                double y = cell_iter.rGetLocation()[1];
            
                double cell_area = mrTissue.rGetVoronoiTessellation().GetFaceArea(node_index);
                double cell_perimeter = mrTissue.rGetVoronoiTessellation().GetFacePerimeter(node_index);
            
                (*mOutStream)<< node_index << " " << x << " " << y << " " << cell_area << " " << cell_perimeter << " ";
                
                if(mFollowLoggedCell)
                {
                    break;
                }
            }
        }
        (*mOutStream)<< "\n";
    }
    
    void WriteTissueAreas()
    {         
        (*mOutStream)<< SimulationTime::Instance()->GetDimensionalisedTime() << " ";
        
        // Don't use the Voronoi tessellation to calculate the total area
        // because it gives huge areas for boundary cells
        double total_area = mrTissue.rGetMesh().CalculateMeshVolume();    
            
        double necrotic_area = 0.0;
        
        for (typename Tissue<DIM>::Iterator cell_iter = mrTissue.Begin();
             cell_iter != mrTissue.End();
             ++cell_iter)
        {
            // Only bother calculating the cell area if it is necrotic
            if (cell_iter->GetCellType() == NECROTIC)
            {
                unsigned node_index = cell_iter.GetNode()->GetIndex();                
                double cell_area = mrTissue.rGetVoronoiTessellation().GetFace(node_index)->GetArea();
                necrotic_area += cell_area;
            }
        }       
        
        (*mOutStream) << total_area << " " << necrotic_area << "\n";
    }
};

#endif /*TISSUEVORONOIDATAWRITER_HPP_*/
