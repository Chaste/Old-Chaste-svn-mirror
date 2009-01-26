/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef ABSTRACTELASTICITYASSEMBLER_HPP_
#define ABSTRACTELASTICITYASSEMBLER_HPP_

#include "AbstractDealiiAssembler.hpp"
#include "DofVertexIterator.hpp"
#include "OutputFileHandler.hpp"
#include <cassert>



template <unsigned DIM>
class AbstractElasticityAssembler : public AbstractDealiiAssembler<DIM>
{

private:
    /*< Data structure containing the deformed position, by vertex index, in easily
     * accessable form. Only created if asked for */
    std::vector<Vector<double> > mDeformedPosition;

    /*< Data structure containing the undeformed position, by vertex index, in easily
     * accessable form. Only created if asked for */
    std::vector<Vector<double> > mUndeformedPosition;


protected:
    /*< Whether to write any output or not */
    bool mWriteOutput;
    /*< Full path of output directory, including chaste testoutput */
    std::string mOutputDirectoryFullPath;

    /*< Number of newton iterations needed to solve the problem, for interest and testing */
    unsigned mNumNewtonIterations;

public:
    /** Constructor
     *
     *  Just takes in the mesh and passes it down to AbstractAssembler
     */
    AbstractElasticityAssembler(Triangulation<DIM>* pMesh, std::string outputDirectory)
        : AbstractDealiiAssembler<DIM>(pMesh)
    {
        if (outputDirectory!="")
        {
            mWriteOutput = true;
            OutputFileHandler output_file_handler(outputDirectory, true); // clean the directory
            mOutputDirectoryFullPath = output_file_handler.GetOutputDirectoryFullPath();
        }
        else
        {
            mWriteOutput = false;
        }

        mNumNewtonIterations = UNSIGNED_UNSET;
    }

    /**
     *  Get the deformed position. rGetDeformedPosition()[i](j) is the x_i value at node j
     */
    std::vector<Vector<double> >& rGetDeformedPosition()
    {
        mDeformedPosition.resize(DIM);
        for (unsigned i=0; i<DIM; i++)
        {
            mDeformedPosition[i].reinit(this->mpMesh->n_vertices());
        }

        DofVertexIterator<DIM> vertex_iter(this->mpMesh, &this->mDofHandler);
        while (!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            Point<DIM> old_posn = vertex_iter.GetVertex();

            for (unsigned i=0; i<DIM; i++)
            {
                mDeformedPosition[i](vertex_index) =   old_posn(i)
                                                     + mCurrentSolution(vertex_iter.GetDof(i));
            }

            vertex_iter.Next();
        }

        return mDeformedPosition;
    }


    /**
     *  Get the undeformed position. rGetUndeformedPosition()[i][j] is the X_i value at node j
     *  Obviously this data is accessible from the mesh as well, this method is more useful
     *  in some situations. Note, this data structure is not set up unless this method
     *  is called.
     *
     *  Note we don't just calculate this once and store because the undeformed mesh will
     *  change if coarsening/refinement occurs
     */
    std::vector<Vector<double> >& rGetUndeformedPosition()
    {
        // initialise
        mUndeformedPosition.resize(DIM);
        for (unsigned i=0; i<DIM; i++)
        {
            mUndeformedPosition[i].reinit(this->mpMesh->n_vertices());
        }

        // populate
        TriangulationVertexIterator<DIM> vertex_iter(this->mpMesh);
        while (!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            Point<DIM> old_posn = vertex_iter.GetVertex();

            for (unsigned i=0; i<DIM; i++)
            {
                mUndeformedPosition[i](vertex_index) = vertex_iter.GetVertex()(i);
            }

            vertex_iter.Next();
        }

        return mUndeformedPosition;
    }


    /**
     *  Output current deformed position to file (or the undeformed mesh, if the
     *  second parameter is set to false)
     *  @counter A number to suffix the file. The output file will be
     *   <out_dir>/solution_<counter.[nodes/elem/undefnodes/undefelem]
     *  @writeDeformed whether to write the deformed position or the undeformed
     *   position, defaults to deformed
     */
    void WriteOutput(unsigned counter, bool writeDeformed=true)
    {
        // only write output if the flag mWriteOutput has been set
        if (!mWriteOutput)
        {
            return;
        }

        /////////////////////////////////////////////////////////////////////
        // create an node file, by looping over vertices and writing
        //   vertex_index x [y [z]]
        /////////////////////////////////////////////////////////////////////
        std::stringstream ss_nodes;
        std::stringstream ss_elem;
        ss_nodes << mOutputDirectoryFullPath << "/solution_" << counter;
        ss_elem  << mOutputDirectoryFullPath << "/solution_" << counter;;

        if(writeDeformed)
        {
            ss_nodes << ".nodes";
            ss_elem  << ".elem";
        }
        else
        {
            ss_nodes << ".undefnodes";
            ss_elem  << ".undefelem";
        }

        std::string nodes_filename = ss_nodes.str();
        std::ofstream nodes_output(nodes_filename.c_str());

        std::vector<Vector<double> >& r_deformed_position = rGetDeformedPosition();
        std::vector<Vector<double> >& r_undeformed_position = rGetUndeformedPosition();

        // loop over nodes in the mesh using the vertex iter
        // NOTE: we don't print out every all of
        // r_deformed_position[i](index) because for some values of index,
        // it will correspond to a non-active node.
        TriangulationVertexIterator<DIM> vertex_iter(this->mpMesh);
        while (!vertex_iter.ReachedEnd())
        {
            unsigned index = vertex_iter.GetVertexGlobalIndex();

            nodes_output << index << " ";
            for (unsigned i=0; i<DIM; i++)
            {
                if(writeDeformed)
                {
                    nodes_output << r_deformed_position[i](index) << " ";
                }
                else
                {
                    nodes_output << r_undeformed_position[i](index) << " ";
                }
            }
            nodes_output << "\n";
            vertex_iter.Next();
        }
        nodes_output.close();

        /////////////////////////////////////////////////////////////////////
        // create an element file, by looping over elements and writing
        //   node1 node2  .... nodeN material_id
        // where node_i is the vertex index
        /////////////////////////////////////////////////////////////////////
        std::string elem_filename = ss_elem.str();
        std::ofstream elem_output(elem_filename.c_str());

        for( typename Triangulation<DIM>::active_cell_iterator element_iter = this->mpMesh->begin_active();
             element_iter!=this->mpMesh->end();
             element_iter++)
        {
            // loop over all vertices..
            for (unsigned i=0; i<GeometryInfo<DIM>::vertices_per_cell; i++)
            {
                elem_output << element_iter->vertex_index(i) << " ";
            }
            unsigned material_id = element_iter->material_id();
            elem_output << material_id << "\n";
        }

        elem_output.close();
    }


    /** Turn writing output off (or back on). Turning back on is only valid
     *  if a non-empty directory was initially specified in the constructor
     */
    void SetWriteOutput(bool writeOutput)
    {
        // check we gave an output dir in the constructor if output is being turned
        // on
        assert(!writeOutput || (writeOutput && (mOutputDirectoryFullPath!="")));
        mWriteOutput = writeOutput;
    }

    /**
     *  Get the number of newton iterations that had been required to solve the problem
     */
    unsigned GetNumNewtonIterations()
    {
        assert(mNumNewtonIterations!=UNSIGNED_UNSET);
        return mNumNewtonIterations;
    }
};
#endif /*ABSTRACTELASTICITYASSEMBLER_HPP_*/
