/*

Copyright (C) University of Oxford, 2005-2011

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
#include "BoxCollection.hpp"
#include "PetscTools.hpp"

/////////////////////////////////////////////////////////////////////////////
// Box methods
/////////////////////////////////////////////////////////////////////////////


template<unsigned DIM>
Box<DIM>::Box(c_vector<double, 2*DIM>& rMinAndMaxValues)
{
    mMinAndMaxValues = rMinAndMaxValues;
}

template<unsigned DIM>
c_vector<double, 2*DIM>& Box<DIM>::rGetMinAndMaxValues()
{
    return mMinAndMaxValues;
}

template<unsigned DIM>
void Box<DIM>::AddNode(Node<DIM>* pNode)
{
    mNodesContained.insert(pNode);
}

template<unsigned DIM>
void Box<DIM>::RemoveNode(Node<DIM>* pNode)
{
    mNodesContained.erase(pNode);
}

template<unsigned DIM>
std::set< Node<DIM>* >& Box<DIM>::rGetNodesContained()
{
    return mNodesContained;
}

template<unsigned DIM>
void Box<DIM>::AddElement(Element<DIM,DIM>* pElement)
{
    mElementsContained.insert(pElement);
}

template<unsigned DIM>
std::set< Element<DIM,DIM>* >& Box<DIM>::rGetElementsContained()
{
    return mElementsContained;
}

/////////////////////////////////////////////////////////////////////////////
// BoxCollection methodsmpDistributedBoxStacks= new DistributedVector(PetscTools::CreateVec(stacks_vector), &factory);
/////////////////////////////////////////////////////////////////////////////


template<unsigned DIM>
BoxCollection<DIM>::BoxCollection(double boxWidth, c_vector<double, 2*DIM> domainSize)
    : mDomainSize(domainSize),
      mBoxWidth(boxWidth)
{
    /*
     * Start by calculating the number of boxes in each direction and total number of boxes.
     * Also create a helper vector of coefficients, whose first entry is 1 and whose i-th
     * entry (for i>1) is the i-th partial product of the vector mNumBoxesEachDirection. This
     * vector of coefficients will be used in the next code block to compute how many boxes
     * along in each dimension each box, given its index.
     *
     * Note: This is information about the whole domain, not just boxes local to this process.
     */
    unsigned num_boxes = 1;
    std::vector<unsigned> coefficients;
    coefficients.push_back(1);

    // Two cases to separate whether the domain size is a multiple of the box width.

    for (unsigned i=0; i<DIM; i++)
    {
		mNumBoxesEachDirection(i) = (unsigned)ceil((domainSize(2*i+1) - domainSize(2*i))/boxWidth);
		if(mNumBoxesEachDirection(i)<1)
		{
			// If domainSize is zero, need to ensure we have a box to avoid errors.
			mNumBoxesEachDirection(i)=1u;
		}
		num_boxes *= mNumBoxesEachDirection(i);
		coefficients.push_back(coefficients[i]*mNumBoxesEachDirection(i));
    }

    /*
     * Set up a PETSc distributed vector (0,1,2,....mNumBixesEachDirection(0)) to divide the
     * stacks of boxes among the processes.
     */

    std::vector<double> stacks_vector;
    for(unsigned i=0;i<mNumBoxesEachDirection(0);i++)
    {
    	stacks_vector.push_back(i);
    }

    DistributedVectorFactory factory(mNumBoxesEachDirection(0));

    mpDistributedBoxStacks= new DistributedVector(PetscTools::CreateVec(stacks_vector), &factory);

    // Assign member variable for # of procs
    mNumProcs=PetscTools::GetNumProcs();

    // Assign left and right processes
    if(PetscTools::GetMyRank()==PetscTools::GetNumProcs()-1)
    {
    	mProcRight=MPI_PROC_NULL;
    }
    else
    {
    	mProcRight=PetscTools::GetMyRank()+1;
    }

    if(PetscTools::GetMyRank()==0)
    {
    	mProcLeft=MPI_PROC_NULL;
    }
    else
    {
    	mProcLeft=PetscTools::GetMyRank()-1;
    }

    /*
     * The boxes that we set up on each process will have x-indexes in
     * the distributed vector for that process.
     *
     * The min and max values for box with index (i,j,k) will be
     * xmin = mBoxWidth*i , xmax = mBoxWidth*(i+1) etc.
     */



	switch(DIM)
	{

		case 1:
		{
		    for (DistributedVector::Iterator stack_index = mpDistributedBoxStacks->Begin();
		         stack_index != mpDistributedBoxStacks->End();
		         ++stack_index)
		    {

		    	//Calculate x-cordinates of the boxes in this stack.
				c_vector<double, 2*DIM> box_coords;
				box_coords[0]=(double)domainSize[0]+mBoxWidth*(*mpDistributedBoxStacks)[stack_index];
				box_coords[1]=(double)domainSize[0]+mBoxWidth*((*mpDistributedBoxStacks)[stack_index]+1);

				Box<DIM> new_box(box_coords);
				mBoxes.push_back(new_box);
				// Insert the local index into the box mapping.
				mBoxesMapping[(*mpDistributedBoxStacks)[stack_index]]=mBoxes.size()-1;
		    }

			break;
		}

		case 2:
		{
			c_vector<double, 2*DIM> box_coords;
			for(unsigned i=0;i<mNumBoxesEachDirection(1);i++)
			{
			    for (DistributedVector::Iterator stack_index = mpDistributedBoxStacks->Begin();
			         stack_index != mpDistributedBoxStacks->End();
			         ++stack_index)
			    {

			    	//Calculate x-cordinates of the boxes in this stack.

					box_coords[0]=(double)domainSize[0]+mBoxWidth*(*mpDistributedBoxStacks)[stack_index];
					box_coords[1]=(double)domainSize[0]+mBoxWidth*((*mpDistributedBoxStacks)[stack_index]+1);

					// Calculate y-coords
					box_coords[2]=domainSize[2]+i*mBoxWidth;
					box_coords[3]=domainSize[2]+(i+1)*mBoxWidth;

					Box<DIM> new_box(box_coords);
					mBoxes.push_back(new_box);

					// Insert the local index into the box mapping.
					mBoxesMapping[(*mpDistributedBoxStacks)[stack_index] + i*mNumBoxesEachDirection(0)]=mBoxes.size()-1;
			    }
			}
			break;
		}
		case 3:
		{
			c_vector<double, 2*DIM> box_coords;

			for(unsigned j=0;j<mNumBoxesEachDirection(2);j++)
			{
				for(unsigned i=0;i<mNumBoxesEachDirection(1);i++)
				{
					// Calculate the y-coords
					box_coords[2]=domainSize[2]+i*mBoxWidth;
					box_coords[3]=domainSize[2]+(i+1)*mBoxWidth;

						for (DistributedVector::Iterator stack_index = mpDistributedBoxStacks->Begin();
							 stack_index != mpDistributedBoxStacks->End();
							 ++stack_index)
						{

							//Calculate x-cordinates of the boxes in this stack.

							box_coords[0]=(double)domainSize[0]+mBoxWidth*(*mpDistributedBoxStacks)[stack_index];
							box_coords[1]=(double)domainSize[0]+mBoxWidth*((*mpDistributedBoxStacks)[stack_index]+1);

							// Calculate the z-coords
							box_coords[4]=domainSize[4]+j*mBoxWidth;
							box_coords[5]=domainSize[4]+(j+1)*mBoxWidth;

							Box<DIM> new_box(box_coords);
							mBoxes.push_back(new_box);

							// Insert the local index into the box mapping.
							mBoxesMapping[(*mpDistributedBoxStacks)[stack_index] + i*mNumBoxesEachDirection(0)+j*mNumBoxesEachDirection(0)*mNumBoxesEachDirection(1)]=mBoxes.size()-1;
						}
				}
			}

			break;
		}

    }

    // Check we have set up the right number of total boxes set up.
    unsigned local_boxes=mBoxes.size();
    unsigned expected_boxes=1;
    for(unsigned i=0;i<DIM;i++)
    {
    	expected_boxes*=mNumBoxesEachDirection[i];
    }
    unsigned total_boxes;

    MPI_Allreduce(&local_boxes,&total_boxes,1,MPI_UNSIGNED,MPI_SUM,PETSC_COMM_WORLD);

	assert(total_boxes==expected_boxes);
	mNumBoxes=total_boxes;

	mAreLocalBoxesSet=false;
}

template<unsigned DIM>
bool BoxCollection<DIM>::GetBoxOwnership(unsigned globalIndex)
{
	if(mBoxesMapping.find(globalIndex)==mBoxesMapping.end())
	{
		return false;
	}
	else
	{
		return true;
	}
}

template<unsigned DIM>
bool BoxCollection<DIM>::GetHaloBoxOwnership(unsigned globalIndex)
{
	if(mHaloBoxesMapping.find(globalIndex)==mHaloBoxesMapping.end())
	{
		return false;
	}
	else
	{
		return true;
	}
}

template<unsigned DIM>
unsigned BoxCollection<DIM>::CalculateGlobalIndex(c_vector<unsigned, DIM> indices)
{
    unsigned containing_box_index = 0;
    for (unsigned i=0; i<DIM; i++)
    {
        unsigned temp = 1;
        for (unsigned j=0; j<i; j++)
        {
            temp *= mNumBoxesEachDirection(j);
        }
        containing_box_index += temp*indices[i];
    }

    return containing_box_index;
}

template<unsigned DIM>
c_vector<unsigned, DIM>  BoxCollection<DIM>::CalculateCoordinateIndices(unsigned globalIndex)
{
	c_vector<unsigned, DIM> indices;

	switch(DIM)
	{
		case 1:
		{
			indices[0]=globalIndex;
			break;
		}
		case 2:
		{
			unsigned remainder=globalIndex % mNumBoxesEachDirection(0);
			indices[0]=remainder;
			indices[1]=(unsigned)(globalIndex/mNumBoxesEachDirection(0));
			break;
		}

		case 3:
		{
			unsigned remainder1=globalIndex % (mNumBoxesEachDirection(0)*mNumBoxesEachDirection(1));
			unsigned remainder2=remainder1 % mNumBoxesEachDirection(0);
			indices[0]=remainder2;
			indices[1]=((globalIndex-indices[0])/mNumBoxesEachDirection(0))%mNumBoxesEachDirection(1);
			indices[2]=((globalIndex-indices[0]-mNumBoxesEachDirection(0)*indices[1])/(mNumBoxesEachDirection(0)*mNumBoxesEachDirection(1)));
			break;
		}
	}

	return indices;
}


template<unsigned DIM>
unsigned BoxCollection<DIM>::CalculateContainingBox(Node<DIM>* pNode)
{
    // Get the location of the node
    c_vector<double, DIM> location = pNode->rGetLocation();
    return CalculateContainingBox(location);
}


template<unsigned DIM>
unsigned BoxCollection<DIM>::CalculateContainingBox(c_vector<double, DIM>& rLocation)
{
	///\todo Deal with the problem case where the width of the cell population co-incides with the MechanicsCutOffLength of the cell population and therefore the width of the boxes.

    // The node must lie inside the boundary of the box collection
    for (unsigned i=0; i<DIM; i++)
    {
        if( (rLocation[i] < mDomainSize(2*i)) || (rLocation[i] > mDomainSize(2*i+1)) )
        {
            EXCEPTION("The point provided is outside all of the boxes");
        }
    }

    // Compute the containing box index in each dimension
    c_vector<unsigned, DIM> containing_box_indices;
    for (unsigned i=0; i<DIM; i++)
    {
        containing_box_indices[i] = (unsigned) floor((rLocation[i] - mDomainSize(2*i))/mBoxWidth);
    }

    // Use these to compute the index of the containing box
    unsigned containing_box_index=CalculateGlobalIndex(containing_box_indices);

    // This index must be less than the number of boxes
    //assert(containing_box_index < mNumBoxes);

    return containing_box_index;
}

template<unsigned DIM>
Box<DIM>& BoxCollection<DIM>::rGetBox(unsigned globalIndex)
{
	// Need to get correct process to do job.
	if(this->GetBoxOwnership(globalIndex))
	{
		std::map<unsigned,unsigned>::const_iterator local_index_it=mBoxesMapping.find(globalIndex);
		return mBoxes[local_index_it->second];
	}
	else
	{
	     EXCEPTION("Box does not exist on this process! Try calling GetBoxOwnership() before rGetBox.");
	}
}

template<unsigned DIM>
unsigned BoxCollection<DIM>::GetNumBoxes()
{
    return mNumBoxes;
}

template <unsigned DIM>
unsigned BoxCollection<DIM>::SolveBoxMapping(unsigned globalIndex) const
{
    std::map<unsigned, unsigned>::const_iterator box_position = mBoxesMapping.find(globalIndex);

    if (box_position == mBoxesMapping.end())
    {
        std::stringstream message;
        message << "Requested box with global index " << globalIndex << ", which does not belong to processor " << PetscTools::GetMyRank();
        EXCEPTION(message.str().c_str());
    }
    return box_position->second;
}


template<unsigned DIM>
void BoxCollection<DIM>::SetupLocalBoxesHalfOnly()
{
	if(mAreLocalBoxesSet)
	{
		EXCEPTION("Local Boxes Are Already Set");
	}
	else
	{
		switch (DIM)
		{
			case 1:
			{
				// We only need to look for neighbours in the current and successive boxes plus some others for halos
				mLocalBoxes.clear();

				// Iterate over the global box indices
				for(std::map<unsigned, unsigned>::iterator it=mBoxesMapping.begin();
					it!=mBoxesMapping.end();
					++it)
				{
					std::set<unsigned> local_boxes;

					// Get this box's global index.
					unsigned global_index=it->first;
					unsigned local_index=it->second;

					// Insert the current box
					local_boxes.insert(global_index);

					// If we're not at the right-most box, then insert the box to the right
					if (global_index < mNumBoxesEachDirection(0)-1)
					{
						local_boxes.insert(global_index+1);
					}
					// If we're on a left process boundary and not on process 0, insert the (halo) box to the left
					if (local_index == 0 && PetscTools::GetMyRank()!=0)
					{
						local_boxes.insert(global_index-1);
					}
					mLocalBoxes.push_back(local_boxes);
				}
				break;
			}
			case 2:
			{
				// We only need to look for neighbours in the current box and half the neighbouring boxes plus some others for halos
				mLocalBoxes.clear();

				for(std::map<unsigned, unsigned>::iterator it=mBoxesMapping.begin();
					it!=mBoxesMapping.end();
					++it)
				{
					std::set<unsigned> local_boxes;

					// Get this box's global index.
					unsigned global_index=it->first;
					unsigned local_index=it->second;
					unsigned stack_width=mpDistributedBoxStacks->GetHigh()-mpDistributedBoxStacks->GetLow();

					// Insert the current box
					local_boxes.insert(global_index);

					// If we're not on the top-most row, then insert the box above
					if (global_index < mNumBoxes - mNumBoxesEachDirection(0))
					{
						local_boxes.insert(global_index + mNumBoxesEachDirection(0));

						// If we're also not on the left-most column, then insert the box above-left
						if (global_index % mNumBoxesEachDirection(0) != 0)
						{
							local_boxes.insert(global_index + mNumBoxesEachDirection(0) - 1);
						}
					}
					// If we're not on the right-most column, then insert the box to the right
					if (global_index % mNumBoxesEachDirection(0) != mNumBoxesEachDirection(0)-1)
					{
						local_boxes.insert(global_index + 1);

						// If we're also not on the top-most row, then insert the box above-right
						if (global_index < mNumBoxes - mNumBoxesEachDirection(0))
						{
							local_boxes.insert(global_index + mNumBoxesEachDirection(0) + 1);
						}
					}
					if (local_index%stack_width==0 && PetscTools::GetMyRank()!=0)
					{
						// insert box to left
						local_boxes.insert(global_index-1);

						// if we're not on the bottom row, insert below left
						if(local_index>0)
						{
							local_boxes.insert(global_index-1-mNumBoxesEachDirection(0));
						}
					}
					if(local_index%stack_width==(stack_width-1) && !PetscTools::AmTopMost())
					{
						//Add lower right box if we're not on the bottom
						if((global_index+1-mNumBoxesEachDirection(0))>0 && (global_index+1-mNumBoxesEachDirection(0))<mNumBoxes)
						{
							local_boxes.insert(global_index+1-mNumBoxesEachDirection(0));
						}
					}
					mLocalBoxes.push_back(local_boxes);
				}
				break;
			}
			case 3:
			{
				// We only need to look for neighbours in the current box and half the neighbouring boxes plus some others for halos
				mLocalBoxes.clear();
				unsigned num_boxes_xy = mNumBoxesEachDirection(0)*mNumBoxesEachDirection(1);


				for(std::map<unsigned, unsigned>::iterator it=mBoxesMapping.begin();
					it!=mBoxesMapping.end();
					++it)
				{
					std::set<unsigned> local_boxes;

					//Get global index
					unsigned global_index=it->first;
					unsigned local_index=it->second;
					unsigned stack_width=mpDistributedBoxStacks->GetHigh()-mpDistributedBoxStacks->GetLow();;

					// Insert the current box
					local_boxes.insert(global_index);

					// If we're not on the far face (y max), then insert the far box
					if (global_index % num_boxes_xy < num_boxes_xy - mNumBoxesEachDirection(0))
					{
						local_boxes.insert(global_index + mNumBoxesEachDirection(0));

						// If we're also not on the left face (x min), then insert the box to the left
						if (global_index % mNumBoxesEachDirection(0) != 0)
						{
							local_boxes.insert(global_index + mNumBoxesEachDirection(0) - 1);
						}
					}
					// If we're not on the right face (x max), then insert the box to the right
					if (global_index % mNumBoxesEachDirection(0) != mNumBoxesEachDirection(0)-1)
					{
						local_boxes.insert(global_index + 1);

						// If we're also not on the far face (y max) row, then insert the box to the far-right
						if (global_index % num_boxes_xy < num_boxes_xy - mNumBoxesEachDirection(0))
						{
							local_boxes.insert(global_index + mNumBoxesEachDirection(0) + 1);
						}
					}
					// If we're not on the top face (z max), then insert the box above
					if (global_index < mNumBoxes - num_boxes_xy)
					{
						local_boxes.insert(global_index + num_boxes_xy);

						// If we're also not on the far face (y max), then insert the above-far box
						if (global_index % num_boxes_xy < num_boxes_xy - mNumBoxesEachDirection(0))
						{
							local_boxes.insert(global_index + num_boxes_xy + mNumBoxesEachDirection(0));

							// If we're also not on the left face (x min), then insert the box to the above-left
							if (global_index % mNumBoxesEachDirection(0) != 0)
							{
								local_boxes.insert(global_index + num_boxes_xy + mNumBoxesEachDirection(0) - 1);
							}
						}
						// If we're also not on the right face (x max), then insert the box to the above-right
						if (global_index % mNumBoxesEachDirection(0) != mNumBoxesEachDirection(0)-1)
						{
							local_boxes.insert(global_index + num_boxes_xy + 1);

							// If we're also not on the far face (y max) row, then insert the box to the above-far-right
							if (global_index % num_boxes_xy < num_boxes_xy - mNumBoxesEachDirection(0))
							{
								local_boxes.insert(global_index + num_boxes_xy + mNumBoxesEachDirection(0) + 1);
							}
						}
					}
					// If we're not on the bottom face (z min), then insert the box above
					if (global_index >= num_boxes_xy)
					{
						local_boxes.insert(global_index - num_boxes_xy);

						// If we're also not on the far face (y max), then insert the below-far box
						if (global_index % num_boxes_xy < num_boxes_xy - mNumBoxesEachDirection(0))
						{
							local_boxes.insert(global_index - num_boxes_xy + mNumBoxesEachDirection(0));

							// If we're also not on the left face (x min), then insert the box to the below-left
							if (global_index % mNumBoxesEachDirection(0) != 0)
							{
								local_boxes.insert(global_index - num_boxes_xy + mNumBoxesEachDirection(0) - 1);
							}
						}
						// If we're also not on the right face (x max), then insert the box to the below-right
						if (global_index % mNumBoxesEachDirection(0) != mNumBoxesEachDirection(0)-1)
						{
							local_boxes.insert(global_index - num_boxes_xy + 1);

							// If we're also not on the far face (y max) row, then insert the box to the below-far-right
							if (global_index % num_boxes_xy < num_boxes_xy - mNumBoxesEachDirection(0))
							{
								local_boxes.insert(global_index - num_boxes_xy + mNumBoxesEachDirection(0) + 1);
							}
						}
					}
					// If we lie on a left-hand boundary with another process, and are not the lowest process
					if( local_index%stack_width==0 && !PetscTools::AmMaster())
					{
						// Include left hand halos

						// 1 to left
						local_boxes.insert(global_index-1);

						// if we're not on the bottom row (ymin)
						if(!(global_index%num_boxes_xy<mNumBoxesEachDirection(0)))
						{
							local_boxes.insert(global_index-1-mNumBoxesEachDirection(0));

							// if we're not on the front face (zmin)
							if(!(global_index<num_boxes_xy))
							{
								local_boxes.insert(global_index-1-num_boxes_xy);
								local_boxes.insert(global_index-1-num_boxes_xy-mNumBoxesEachDirection(0));
							}

							// if we're not on the back face (zmax)
							if(global_index<mNumBoxes-num_boxes_xy)
							{
								local_boxes.insert(global_index-1+num_boxes_xy);
								local_boxes.insert(global_index-1+num_boxes_xy-mNumBoxesEachDirection(0));
							}
						}
						else
						{
							// we're on the bottom row

							// if we're not on the front face (zmin)
							if(!(global_index<num_boxes_xy))
							{
								local_boxes.insert(global_index-1-num_boxes_xy);
							}

							// if we're not on the back face (zmax)
							if(global_index<mNumBoxes-num_boxes_xy)
							{
								local_boxes.insert(global_index-1+num_boxes_xy);
							}
						}
					}
					// If we lie on a right-hand boudnary with another process, and are not the top most process
					if( local_index%stack_width==(stack_width-1) && !PetscTools::AmTopMost() )
					{
						// Include right hand halos

						// if we're not on the bottom row (ymin)
						if(!(global_index%num_boxes_xy<mNumBoxesEachDirection(0)))
						{
							local_boxes.insert(global_index+1-mNumBoxesEachDirection(0));

							// if we're not on the front face (zmin)
							if(!(global_index<num_boxes_xy))
							{
								local_boxes.insert(global_index+1-num_boxes_xy-mNumBoxesEachDirection(0));
							}

							// if we're not on the back face (zmax)
							if(global_index<mNumBoxes-num_boxes_xy)
							{
								local_boxes.insert(global_index+1+num_boxes_xy-mNumBoxesEachDirection(0));
							}
						}

					}
					mLocalBoxes.push_back(local_boxes);
				}
				break;
			}
			default:
				NEVER_REACHED;
		}

		mAreLocalBoxesSet=true;
	}
}



template<unsigned DIM>
void BoxCollection<DIM>::SetupAllLocalBoxes()
{
	if(!PetscTools::IsSequential())
	{
		EXCEPTION("Using SetupAllLocalBoxes() in parallel. This will likely lead to errors. If possible use SetupLocalBoxesHalfOnly()");
	}
	if(mAreLocalBoxesSet)
	{
		EXCEPTION("Local Boxes Already Set");
	}
	else
	{
		switch (DIM)
		{
			case 1:
			{
				for(std::map<unsigned, unsigned>::iterator it=mBoxesMapping.begin();
					it!=mBoxesMapping.end();
					++it)
				{
					std::set<unsigned> local_boxes;

					// Get global index
					unsigned global_index=it->first;
					local_boxes.insert(global_index);

					// add the two neighbours
					if(global_index!=0)
					{
						local_boxes.insert(global_index-1);
					}
					if(global_index+1 != mNumBoxesEachDirection(0))
					{
						local_boxes.insert(global_index+1);
					}

					mLocalBoxes.push_back(local_boxes);
				}
				break;
			}
			case 2:
			{
				mLocalBoxes.clear();

				unsigned M = mNumBoxesEachDirection(0);
				unsigned N = mNumBoxesEachDirection(1);

				std::vector<bool> is_xmin(N*M); // far left
				std::vector<bool> is_xmax(N*M); // far right
				std::vector<bool> is_ymin(N*M); // bottom
				std::vector<bool> is_ymax(N*M); // top

				for(unsigned i=0; i<M*N; i++)
				{
					is_xmin[i] = (i%M==0);
					is_xmax[i] = ((i+1)%M==0);
					is_ymin[i] = (i%(M*N)<M);
					is_ymax[i] = (i%(M*N)>=(N-1)*M);
				}

				for(std::map<unsigned, unsigned>::iterator it=mBoxesMapping.begin();
					it!=mBoxesMapping.end();
					++it)
				{
					std::set<unsigned> local_boxes;

					//Get global index
					unsigned global_index=it->first;

					local_boxes.insert(global_index);

					// add the box to the left
					if(!is_xmin[global_index])
					{
						local_boxes.insert(global_index-1);
					}

					// add the box to the right
					if(!is_xmax[global_index])
					{
						local_boxes.insert(global_index+1);
					}

					// add the one below
					if(!is_ymin[global_index])
					{
						local_boxes.insert(global_index-M);
					}

					// add the one above
					if(!is_ymax[global_index])
					{
						local_boxes.insert(global_index+M);
					}

					// add the four corner boxes

					if( (!is_xmin[global_index]) && (!is_ymin[global_index]) )
					{
						local_boxes.insert(global_index-1-M);
					}

					if( (!is_xmin[global_index]) && (!is_ymax[global_index]) )
					{
						local_boxes.insert(global_index-1+M);
					}

					if( (!is_xmax[global_index]) && (!is_ymin[global_index]) )
					{
						local_boxes.insert(global_index+1-M);
					}

					if( (!is_xmax[global_index]) && (!is_ymax[global_index]) )
					{
						local_boxes.insert(global_index+1+M);
					}

					mLocalBoxes.push_back(local_boxes);
				}
				break;
			}
			case 3:
			{
				mLocalBoxes.clear();

				unsigned M = mNumBoxesEachDirection(0);
				unsigned N = mNumBoxesEachDirection(1);
				unsigned P = mNumBoxesEachDirection(2);

				std::vector<bool> is_xmin(N*M*P); // far left
				std::vector<bool> is_xmax(N*M*P); // far right
				std::vector<bool> is_ymin(N*M*P); // nearest
				std::vector<bool> is_ymax(N*M*P); // furthest
				std::vector<bool> is_zmin(N*M*P); // bottom layer
				std::vector<bool> is_zmax(N*M*P); // top layer

				for(unsigned i=0; i<M*N*P; i++)
				{
					is_xmin[i] = (i%M==0);
					is_xmax[i] = ((i+1)%M==0);
					is_ymin[i] = (i%(M*N)<M);
					is_ymax[i] = (i%(M*N)>=(N-1)*M);
					is_zmin[i] = (i<M*N);
					is_zmax[i] = (i>=M*N*(P-1));
				}

				for(std::map<unsigned, unsigned>::iterator it=mBoxesMapping.begin();
					it!=mBoxesMapping.end();
					++it)
				{
					std::set<unsigned> local_boxes;

					//Get global index
					unsigned global_index=it->first;

					// add itself as a local box
					local_boxes.insert(global_index);

					// now add all 26 other neighbours.....

					// add the box left
					if(!is_xmin[global_index])
					{
						local_boxes.insert(global_index-1);

						// plus some others towards the left
						if(!is_ymin[global_index])
						{
							local_boxes.insert(global_index-1-M);
						}

						if(!is_ymax[global_index])
						{
							local_boxes.insert(global_index-1+M);
						}

						if(!is_zmin[global_index])
						{
							local_boxes.insert(global_index-1-M*N);
						}

						if(!is_zmax[global_index])
						{
							local_boxes.insert(global_index-1+M*N);
						}
					}

					// add the box to the right
					if(!is_xmax[global_index])
					{
						local_boxes.insert(global_index+1);

						// plus some others towards the right
						if(!is_ymin[global_index])
						{
							local_boxes.insert(global_index+1-M);
						}

						if(!is_ymax[global_index])
						{
							local_boxes.insert(global_index+1+M);
						}

						if(!is_zmin[global_index])
						{
							local_boxes.insert(global_index+1-M*N);
						}

						if(!is_zmax[global_index])
						{
							local_boxes.insert(global_index+1+M*N);
						}
					}

					// add the boxes next along the y axis
					if(!is_ymin[global_index])
					{
						local_boxes.insert(global_index-M);

						// and more in this plane
						if(!is_zmin[global_index])
						{
							local_boxes.insert(global_index-M-M*N);
						}

						if(!is_zmax[global_index])
						{
							local_boxes.insert(global_index-M+M*N);
						}
					}

					// add the boxes next along the y axis
					if(!is_ymax[global_index])
					{
						local_boxes.insert(global_index+M);

						// and more in this plane
						if(!is_zmin[global_index])
						{
							local_boxes.insert(global_index+M-M*N);
						}

						if(!is_zmax[global_index])
						{
							local_boxes.insert(global_index+M+M*N);
						}
					}

					// add the box directly above
					if(!is_zmin[global_index])
					{
						local_boxes.insert(global_index-N*M);
					}

					// add the box directly below
					if(!is_zmax[global_index])
					{
						local_boxes.insert(global_index+N*M);
					}

					// finally, the 8 corners are left

					if( (!is_xmin[global_index]) && (!is_ymin[global_index]) && (!is_zmin[global_index]) )
					{
						local_boxes.insert(global_index-1-M-M*N);
					}

					if( (!is_xmin[global_index]) && (!is_ymin[global_index]) && (!is_zmax[global_index]) )
					{
						local_boxes.insert(global_index-1-M+M*N);
					}

					if( (!is_xmin[global_index]) && (!is_ymax[global_index]) && (!is_zmin[global_index]) )
					{
						local_boxes.insert(global_index-1+M-M*N);
					}

					if( (!is_xmin[global_index]) && (!is_ymax[global_index]) && (!is_zmax[global_index]) )
					{
						local_boxes.insert(global_index-1+M+M*N);
					}

					if( (!is_xmax[global_index]) && (!is_ymin[global_index]) && (!is_zmin[global_index]) )
					{
						local_boxes.insert(global_index+1-M-M*N);
					}

					if( (!is_xmax[global_index]) && (!is_ymin[global_index]) && (!is_zmax[global_index]) )
					{
						local_boxes.insert(global_index+1-M+M*N);
					}

					if( (!is_xmax[global_index]) && (!is_ymax[global_index]) && (!is_zmin[global_index]) )
					{
						local_boxes.insert(global_index+1+M-M*N);
					}

					if( (!is_xmax[global_index]) && (!is_ymax[global_index]) && (!is_zmax[global_index]) )
					{
						local_boxes.insert(global_index+1+M+M*N);
					}

					mLocalBoxes.push_back(local_boxes);
				}
				break;
			}
			default:
				NEVER_REACHED;
		}

		mAreLocalBoxesSet=true;
	}
}

template<unsigned DIM>
void BoxCollection<DIM>::SetupHaloBoxes()
{
	// Get top-most and bottom-most value of Distributed Box Stack.
	unsigned Hi=mpDistributedBoxStacks->GetHigh();
	unsigned Lo=mpDistributedBoxStacks->GetLow();

	c_vector<double, 2*DIM> box_coords;

	switch(DIM)
	{
		case 1:
		{
			// Different cases for Master and Top-Most Processes
			if(!PetscTools::AmTopMost())
			{
				// Right Stack
				box_coords[0]=(double)mDomainSize[0]+mBoxWidth*Hi;
				box_coords[1]=(double)mDomainSize[0]+mBoxWidth*(Hi+1);

				Box<DIM> new_box2(box_coords);
				mHaloBoxes.push_back(new_box2);
				mHaloBoxesMapping[Hi]=mHaloBoxes.size()-1;

				mHalosRight.push_back(Hi-1);
			}
			if(!PetscTools::AmMaster())
			{
				// Left Stack
				box_coords[0]=(double)mDomainSize[0]+mBoxWidth*(Lo-1);
				box_coords[1]=(double)mDomainSize[0]+mBoxWidth*(Lo);

				Box<DIM> new_box(box_coords);
				mHaloBoxes.push_back(new_box);
				mHaloBoxesMapping[Lo-1]=mHaloBoxes.size()-1;

				mHalosLeft.push_back(Lo);
			}
	    	break;
		}

		case 2:
		{
			if(!PetscTools::AmTopMost())
			{
				// Hi
				for(unsigned i=0;i<mNumBoxesEachDirection(1);i++)
				{
						//Calculate x-cordinates of the boxes in this stack.
						box_coords[0]=(double)mDomainSize[0]+mBoxWidth*Hi;
						box_coords[1]=(double)mDomainSize[0]+mBoxWidth*(Hi+1);

						// Calculate y-coords
						box_coords[2]=mDomainSize[2]+i*mBoxWidth;
						box_coords[3]=mDomainSize[2]+(i+1)*mBoxWidth;

						Box<DIM> new_box(box_coords);
						mHaloBoxes.push_back(new_box);
						mHaloBoxesMapping[Hi+i*mNumBoxesEachDirection(0)]=mHaloBoxes.size()-1;

						mHalosRight.push_back(Hi+i*mNumBoxesEachDirection(0)-1);
				}
			}
			if(!PetscTools::AmMaster())
			{
				// Lo
				for(unsigned i=0;i<mNumBoxesEachDirection(1);i++)
				{
						//Calculate x-cordinates of the boxes in this stack.
						box_coords[0]=(double)mDomainSize[0]+mBoxWidth*(Lo-1);
						box_coords[1]=(double)mDomainSize[0]+mBoxWidth*(Lo);

						// Calculate y-coords
						box_coords[2]=mDomainSize[2]+i*mBoxWidth;
						box_coords[3]=mDomainSize[2]+(i+1)*mBoxWidth;

						Box<DIM> new_box(box_coords);
						mHaloBoxes.push_back(new_box);
						mHaloBoxesMapping[Lo+i*mNumBoxesEachDirection(0)-1]=mHaloBoxes.size()-1;

						mHalosLeft.push_back(Lo+i*mNumBoxesEachDirection(0));

				}
			}
			break;
		}

		case 3:
		{
			if(!PetscTools::AmTopMost())
			{
				// Hi
				for(unsigned j=0;j<mNumBoxesEachDirection(2);j++)
				{
					for(unsigned i=0;i<mNumBoxesEachDirection(1);i++)
					{
						// Calculate the y-coords
						box_coords[2]=mDomainSize[2]+i*mBoxWidth;
						box_coords[3]=mDomainSize[2]+(i+1)*mBoxWidth;

						//Calculate x-cordinates of the boxes in this stack.
						box_coords[0]=(double)mDomainSize[0]+mBoxWidth*Hi;
						box_coords[1]=(double)mDomainSize[0]+mBoxWidth*(Hi+1);

						// Calculate the z-coords
						box_coords[4]=mDomainSize[4]+j*mBoxWidth;
						box_coords[5]=mDomainSize[4]+(j+1)*mBoxWidth;

						Box<DIM> new_box(box_coords);
						mHaloBoxes.push_back(new_box);
						mHaloBoxesMapping[Hi+i*mNumBoxesEachDirection(0)+j*mNumBoxesEachDirection(1)*mNumBoxesEachDirection(0)]=mHaloBoxes.size()-1;

						mHalosRight.push_back(Hi+i*mNumBoxesEachDirection(0)+j*mNumBoxesEachDirection(1)*mNumBoxesEachDirection(0)-1);
					}
				}
			}
			if(!PetscTools::AmMaster())
			{
				// Lo
				for(unsigned j=0;j<mNumBoxesEachDirection(2);j++)
				{
					for(unsigned i=0;i<mNumBoxesEachDirection(1);i++)
					{
						// Calculate the y-coords
						box_coords[2]=mDomainSize[2]+i*mBoxWidth;
						box_coords[3]=mDomainSize[2]+(i+1)*mBoxWidth;


						//Calculate x-cordinates of the boxes in this stack.
						box_coords[0]=(double)mDomainSize[0]+mBoxWidth*(Lo-1);
						box_coords[1]=(double)mDomainSize[0]+mBoxWidth*Lo;

						// Calculate the z-coords
						box_coords[4]=mDomainSize[4]+j*mBoxWidth;
						box_coords[5]=mDomainSize[4]+(j+1)*mBoxWidth;

						Box<DIM> new_box(box_coords);
						mHaloBoxes.push_back(new_box);
						mHaloBoxesMapping[Lo+i*mNumBoxesEachDirection(0)+j*mNumBoxesEachDirection(1)*mNumBoxesEachDirection(0)-1]=mHaloBoxes.size()-1;

						mHalosLeft.push_back(Lo+i*mNumBoxesEachDirection(0)+j*mNumBoxesEachDirection(1)*mNumBoxesEachDirection(0));
					}
				}
			}
		break;
		}

	}
}

template<unsigned DIM>
void BoxCollection<DIM>::UpdateHaloBoxes()
{

	MPI_Status status;

	// Clear halo nodes

	for(unsigned i=0;i<mHaloBoxes.size();i++)
	{
		for(typename std::set<Node<DIM>*>::iterator it=mHaloBoxes[i].rGetNodesContained().begin();
				it!=mHaloBoxes[i].rGetNodesContained().end();
				it++)
		{
			mHaloBoxes[i].RemoveNode(*it);
		}
	}

	// Send # local halo nodes left and right

	unsigned size_send_left=0;
	unsigned size_recv_left=0;
	unsigned size_send_right=0;
	unsigned size_recv_right=0;

	// Calculate number of cells to send in each direction and pass to neighbouring process.

	for(unsigned i=0;i<mHalosRight.size();i++)
	{
		size_send_right+=this->rGetBox(mHalosRight[i]).rGetNodesContained().size();
	}

	MPI_Sendrecv(	&size_send_right,
					1,
					MPI_UNSIGNED,
					mProcRight,
					123,
					&size_recv_right,
					1,
					MPI_UNSIGNED,
					mProcRight,
					123,
					PETSC_COMM_WORLD,
					&status);


	for(unsigned i=0;i<mHalosLeft.size();i++)
	{
		size_send_left+=this->rGetBox(mHalosLeft[i]).rGetNodesContained().size();
	}
	MPI_Sendrecv(	&size_send_left,
					1,
					MPI_UNSIGNED,
					mProcLeft,
					123,
					&size_recv_left,
					1,
					MPI_UNSIGNED,
					mProcLeft,
					123,
					PETSC_COMM_WORLD,
					&status);




	// Pack node location and index data

	double send_left_data[DIM*size_send_left];
	unsigned send_left_indices[size_send_left];
	double send_right_data[DIM*size_send_right];
	unsigned send_right_indices[size_send_right];


	double receive_data_from_left[DIM*size_recv_left];
	double receive_data_from_right[DIM*size_recv_right];
	unsigned indices_recvd_from_right[size_recv_right];
	unsigned indices_recvd_from_left[size_recv_left];

	unsigned node_counter=0;
	for(unsigned i=0;i<mHalosLeft.size();i++)
	{

		for(typename std::set< Node<DIM>* >::iterator it=this->rGetBox(mHalosLeft[i]).rGetNodesContained().begin();
				it!=this->rGetBox(mHalosLeft[i]).rGetNodesContained().end();
				it++)
		{
			send_left_indices[node_counter]=(*it)->GetIndex();
			for(unsigned j=0;j<DIM;j++)
			{
				send_left_data[DIM*node_counter+j]=(*it)->rGetLocation()[j];
			}
			node_counter++;
		}
	}
	node_counter=0;
	for(unsigned i=0;i<mHalosRight.size();i++)
	{

		for(typename std::set< Node<DIM>* >::iterator it=this->rGetBox(mHalosRight[i]).rGetNodesContained().begin();
				it!=this->rGetBox(mHalosRight[i]).rGetNodesContained().end();
				it++)
		{
			send_right_indices[node_counter]=(*it)->GetIndex();
			for(unsigned j=0;j<DIM;j++)
			{
				send_right_data[DIM*node_counter+j]=(*it)->rGetLocation()[j];
			}
			node_counter++;

		}
	}

	// Send node data left and right

	MPI_Sendrecv( 	send_right_data,
					DIM*size_send_right,
					MPI_DOUBLE,
					mProcRight,
					0,
					&receive_data_from_right,
					DIM*size_recv_right,
					MPI_DOUBLE,
					mProcRight,
					0,
					PETSC_COMM_WORLD,
					&status);

	MPI_Sendrecv( 	send_left_data,
					DIM*size_send_left,
					MPI_DOUBLE,
					mProcLeft,
					0,
					&receive_data_from_left,
					DIM*size_recv_left,
					MPI_DOUBLE,
					mProcLeft,
					0,
					PETSC_COMM_WORLD,
					&status);

	MPI_Sendrecv( 	send_right_indices,
					size_send_right,
					MPI_UNSIGNED,
					mProcRight,
					0,
					&indices_recvd_from_right,
					size_recv_right,
					MPI_UNSIGNED,
					mProcRight,
					0,
					PETSC_COMM_WORLD,
					&status);

	MPI_Sendrecv( 	send_left_indices,
					size_send_left,
					MPI_UNSIGNED,
					mProcLeft,
					0,
					&indices_recvd_from_left,
					size_recv_left,
					MPI_UNSIGNED,
					mProcLeft,
					0,
					PETSC_COMM_WORLD,
					&status);


	// Unpack node data & add new nodes to halo boxes.
	for(unsigned i=0;i<size_recv_right;i++)
	{
		c_vector<double, DIM> location;
		for(unsigned d=0;d<DIM;d++)
		{
			location[d]=receive_data_from_right[DIM*i+d];
		}

		//Calculate containing box
		unsigned box=this->CalculateContainingBox(location);
		mHaloBoxes[mHaloBoxesMapping[box]].AddNode(new Node<DIM>(indices_recvd_from_right[i], location));
	}
	for(unsigned i=0;i<size_recv_left;i++)
	{
		c_vector<double, DIM> location;
		for(unsigned d=0;d<DIM;d++)
		{
			location[d]=receive_data_from_left[DIM*i+d];
		}

		//Calculate containing box
		unsigned box=this->CalculateContainingBox(location);
		mHaloBoxes[mHaloBoxesMapping[box]].AddNode(new Node<DIM>(indices_recvd_from_left[i], location));
	}

}

template<unsigned DIM>
std::set<unsigned> BoxCollection<DIM>::GetLocalBoxes(unsigned globalIndex)
{

	// Make sure it lies on this process
	if(mBoxesMapping.find(globalIndex)==mBoxesMapping.end())
	{
		EXCEPTION("Called GetLocalBoxes on a boxIndex that does not belong to this process");
	}

	return mLocalBoxes[mBoxesMapping[globalIndex]];

}

template<unsigned DIM>
void BoxCollection<DIM>::CalculateNodePairs(std::vector<Node<DIM>*>& rNodes, std::set<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs)
{
    rNodePairs.clear();

    // Make sure that halo boxes are up to date

    // \todo Could remove this if called elsewhere in the code - just a safety net for now

    // Put the nodes in their appropriate boxes
    for (unsigned node_index=0; node_index<rNodes.size(); node_index++)
    {
		unsigned box_index = CalculateContainingBox(rNodes[node_index]);
		// Check that this box or its halo is owned by the process. If not, ignore it - it will be picked up by another process.
		if(this->GetBoxOwnership(box_index))
		{
			mBoxes[mBoxesMapping[box_index]].AddNode(rNodes[node_index]);
		}
    }

    this->UpdateHaloBoxes();

    // Loop over global nodes. (Global set of nodes is passed in).

    for (unsigned node_index=0; node_index<rNodes.size(); node_index++)
    {
        // Get the box containing this node
        unsigned box_index = CalculateContainingBox(rNodes[node_index]);

		// Make sure we own the box before we do anything with it.
		if(this->GetBoxOwnership(box_index))
		{
			// Get the global indices of local boxes to this box - including halo boxes.
			std::set<unsigned> local_boxes_indices = GetLocalBoxes(box_index);

			// Loop over all the local boxes
			for (std::set<unsigned>::iterator box_iter = local_boxes_indices.begin();
				 box_iter != local_boxes_indices.end();
				 box_iter++)
			{
				// Establish whether this box is owned.
				bool is_owned=this->GetBoxOwnership(*box_iter);
				// If it's not owned assume it is a halo box

				if(!is_owned)
				{
					// Check it really is a halo box
					assert(this->GetHaloBoxOwnership(*box_iter));

					// If it is not owned then two boxes cannot be the same
					assert(*box_iter != box_index);

					// Get the set of nodes contained in this box
					std::set< Node<DIM>* >& r_contained_nodes = mHaloBoxes[mHaloBoxesMapping[*box_iter]].rGetNodesContained();

					// Loop over these nodes
					for (typename std::set<Node<DIM>*>::iterator node_iter = r_contained_nodes.begin();
						 node_iter != r_contained_nodes.end();
						 ++node_iter)
					{
						// Get the index of the other node
						unsigned other_node_index = (*node_iter)->GetIndex();

						rNodePairs.insert(std::pair<Node<DIM>*, Node<DIM>*>(rNodes[node_index], rNodes[other_node_index]));
					}
				}

				// Otherwise it is owned
				else
				{
					// Get the set of nodes contained in this box
					std::set< Node<DIM>* >& r_contained_nodes= mBoxes[mBoxesMapping[*box_iter]].rGetNodesContained();

					for (typename std::set<Node<DIM>*>::iterator node_iter = r_contained_nodes.begin();
						 node_iter != r_contained_nodes.end();
						 ++node_iter)
					{
						// Get the index of the other node
						unsigned other_node_index = (*node_iter)->GetIndex();

						// If we're in the same box, then take care not to store the node pair twice
						if (*box_iter == box_index)
						{
							if (other_node_index > node_index)
							{
								rNodePairs.insert(std::pair<Node<DIM>*, Node<DIM>*>(rNodes[node_index], rNodes[other_node_index]));
							}
						}
						else
						{
							rNodePairs.insert(std::pair<Node<DIM>*, Node<DIM>*>(rNodes[node_index], rNodes[other_node_index]));
						}
					}

				}
			}
		}
    }
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class Box<1>;
template class Box<2>;
template class Box<3>;
template class BoxCollection<1>;
template class BoxCollection<2>;
template class BoxCollection<3>;
