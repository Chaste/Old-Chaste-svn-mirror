#ifndef _ELEMENT_HPP_
#define _ELEMENT_HPP_

#include "AbstractElement.cpp"
#include <set>

template <int ELEMENT_DIM, int SPACE_DIM>
class Element : public AbstractElement<ELEMENT_DIM, SPACE_DIM>
{

public:
    Element(unsigned index,
            std::vector<Node<SPACE_DIM>*> nodes,
            unsigned orderOfBasisFunctions=1): AbstractElement<ELEMENT_DIM, SPACE_DIM>(index,nodes,orderOfBasisFunctions)
    {
        this->mIsDeleted=false;
        RegisterWithNodes();
    }
    
    /***
     * Copy constructor which allows a new index to be specified
     */
    Element(const Element &element, const unsigned index)
    {
        this->mIndex=index;
        CommonConstructor(element);
        RegisterWithNodes();
    }
    
    void RegisterWithNodes()
    {
        for (unsigned i=0; i<this->mNodes.size(); i++)
        {
            this->mNodes[i]->AddElement(this->mIndex);
        }
    }
    
    void MarkAsDeleted()
    {
        this->mIsDeleted = true;
        this->mJacobianDeterminant = 0.0;
        // Update nodes in this element so they know they are not contained by us
        for (unsigned i=0; i<this->GetNumNodes(); i++)
        {
            this->mNodes[i]->RemoveElement(this->mIndex);
        }
    }
    
    /** Update node at the given index
     *  @param rIndex is an local index to which node to change
     *  @param pNode is a pointer to the replacement node
     */
    void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode)
    {
        assert(rIndex < this->mNodes.size());
        
        // Remove it from the node at this location
        this->mNodes[rIndex]->RemoveElement(this->mIndex);
        
        // Update the node at this location
        this->mNodes[rIndex] = pNode;
        
        // Add element to this node
        this->mNodes[rIndex]->AddElement(this->mIndex);
    }
    
    void ResetIndex(unsigned index)
    {
        for (unsigned i=0; i<this->GetNumNodes(); i++)
        {
            this->mNodes[i]->RemoveElement(this->mIndex);
        }
        this->mIndex=index;
        RegisterWithNodes();
    }
    
    /*Returns a vector representing the circumsphere/circumcircle
     * @returns a vector containing x_centre, y_centre,...,radius^2
    */
    c_vector<double,SPACE_DIM+1> CalculateCircumsphere()
    {
        /*Assuming that x0,y0.. is at the origin then we need to solve
         *
         * ( 2x1 2y1 2z1  ) (x)    (x1^2+y1^2+z1^2)
         * ( 2x2 2y2 2z2  ) (y)    (x2^2+y2^2+z2^2)
         * ( 2x3 2y3 2z3  ) (z)    (x3^2+y3^2+z3^2)
         * where (x,y,z) is the circumcentre
         * 
         */
        assert (ELEMENT_DIM == SPACE_DIM);
        c_vector <double, ELEMENT_DIM> rhs;
        
        for (unsigned j=0; j<ELEMENT_DIM; j++)
        {
            double squared_location=0.0;
            for (unsigned i=0; i<SPACE_DIM; i++)
            {
                //mJacobian(i,j) is the i-th component of j-th vertex (relative to vertex 0)
                squared_location += this->mJacobian(i,j)*this->mJacobian(i,j);
            }
            rhs[j]=squared_location/2.0;
        }
        
        c_vector <double, ELEMENT_DIM> centre=prod(rhs, this->mInverseJacobian);
        c_vector <double, ELEMENT_DIM+1> circum;
        double squared_radius=0.0;
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            circum[i]=centre[i] + this->GetNodeLocation(0,i);
            squared_radius += centre[i]*centre[i];
        }
        circum[SPACE_DIM]=squared_radius;
          
        return circum;
        
    }
            
    double CalculateCircumsphereVolume()
    {
        c_vector<double, SPACE_DIM+1> circum=CalculateCircumsphere();
        if (SPACE_DIM == 1)
        {
            return 2.0*sqrt(circum[SPACE_DIM]); //2*r
        }
        else if (SPACE_DIM == 2)
        {
            return M_PI*circum[SPACE_DIM]; //Pi*r^2
        }
        assert (SPACE_DIM == 3);
        return 4.0*M_PI*circum[SPACE_DIM]*sqrt(circum[SPACE_DIM])/3.0; //4*Pi*r^3/3
    }
    
    /* The quality of a triangle/tetrahedron is the ratio between the 
     * volume of the shape and the volume of its circumsphere.
     * This is normalised by dividing through by the Platonic ratio.
     */
     
    double CalculateQuality()
    {
        assert (SPACE_DIM == ELEMENT_DIM);
        if (SPACE_DIM == 1)
        {
            return 1.0;
        }
        
        c_vector<double, SPACE_DIM+1> circum=CalculateCircumsphere();
        if (SPACE_DIM == 2)
        {
            /* Want Q=(Area_Tri / Area_Cir) / (Area_Equilateral_Tri / Area_Equilateral_Cir)
             * Area_Tri = |Jacobian| /2
             * Area_Cir = Pi * r^2
             * Area_Eq_Tri = (3*sqrt(3)/4)*R^2
             * Area_Eq_Tri = Pi * R^2
             * Q= (2*|Jacobian|)/ (    bool CalculateVoronoiElement(c_vector <double, 3> first_node, c_vector <double, 3> second_node)
    {
        double x_diff_sqr = ((first_node[0] - second_node[0])*(first_node[0] - second_node[0]));
        double y_diff_sqr = ((first_node[1] - second_node[1])*(first_node[1] - second_node[1]));
        
        return ((x_diff_sqr + y_diff_sqr) > first_node[2]);
    }3*sqrt(3)*r^2)
             */
            return 2.0*this->mJacobianDeterminant/(3.0*sqrt(3)*circum[SPACE_DIM]);
        }
        assert (SPACE_DIM == 3);
       /* Want Q=(Vol_Tet / Vol_CirS) / (Vol_Plat_Tet / Vol_Plat_CirS)
         *  Vol_Tet  = |Jacobian| /6
         *  Vol_CirS = 4*Pi*r^3/3
         *  Vol_Plat_Tet  = 8*sqrt(3)*R^3/27
         *  Vol_Plat_CirS = 4*Pi*R^3/3
        * Q= 3*sqrt(3)*|Jacobian|/ (16*r^3)
         */
        
        return (3.0*sqrt(3.0)*this->mJacobianDeterminant)
                /(16.0*circum[SPACE_DIM]*sqrt(circum[SPACE_DIM]));          
    }
    
    
    c_vector<double, SPACE_DIM+1> CalculateInterpolationWeights(Point <SPACE_DIM> testPoint)
    {
        //Can only test if it's a tetrahedal mesh in 3d, triangles in 2d...
        assert (ELEMENT_DIM == SPACE_DIM);      
        
        c_vector<double, SPACE_DIM+1> weights;
        
        //Find the location with respect to node 0
        c_vector<double, SPACE_DIM> test_location=testPoint.rGetLocation()-this->GetNodeLocation(0);
        
        //Multiply by inverse Jacobian
        c_vector<double, SPACE_DIM> ans=prod(this->mInverseJacobian, test_location);
        
        //Copy 3 weights and compute the fourth weight
        weights[0]=1.0;
        for (unsigned i=1; i<=SPACE_DIM; i++)
        {
            weights[0] -= ans[i-1];
            weights[i] = ans[i-1];
        }
        return weights;
    }
    
    bool IncludesPoint(Point<SPACE_DIM> testPoint, bool strict=false)
    {
	   //Can only test if it's a tetrahedal mesh in 3d, triangles in 2d...
        assert (ELEMENT_DIM == SPACE_DIM);

        c_vector<double, SPACE_DIM+1> weights=CalculateInterpolationWeights(testPoint);
		
        //If the point is in the simplex then all the weights should be positive
        
        for (unsigned i=0;i<=SPACE_DIM;i++)
        {
            if (strict)
            {
                //Points can't be close to a face
                if (weights[i] <= DBL_EPSILON)
                {
                    return false;
                }
            }
            else
            {
                //Allow point to be close to a face
                if (weights[i] < -DBL_EPSILON)
                {
                    return false;
                }
            }
        }	    	
    	return true;
    }

    
};



#endif //_BOUNDARYELEMENT_HPP_

