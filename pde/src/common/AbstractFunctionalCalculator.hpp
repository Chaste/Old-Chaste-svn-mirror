#ifndef ABSTRACTFUNCTIONALCALCULATOR_HPP_
#define ABSTRACTFUNCTIONALCALCULATOR_HPP_

#include "LinearBasisFunction.hpp"
#include "GaussianQuadratureRule.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "GaussianQuadratureRule.hpp"
#include "ReplicatableVector.hpp"


/**
 *  This is an abstract class for computing user-defined integral-based 
 *  functionals of a solution on the mesh. The user needs to define 
 *  GetIntegrand() in the concrete class, and this class can then be used
 *  to compute the integral of f(x,u,grad_u) over the mesh, where x is 
 *  position, u is the solution at x (possibly with multiple components, for
 *  which PROBLEM_DIM>1), grad_u the gradient of u and f the integrand as 
 *  defined in the concrete class.
 * 
 *  Note linear basis functions and 2 quad points per dimension are currently
 *  hardcoded.
 */ 
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractFunctionalCalculator
{
private:
    /** replicated store of the solution vector */
    ReplicatableVector mSolutionReplicated;

    /** The integrand. Must be defined by the user */ 
    virtual double GetIntegrand(ChastePoint<SPACE_DIM> &rX,
                                c_vector<double,PROBLEM_DIM> &rU,
                                c_matrix<double,PROBLEM_DIM,SPACE_DIM> &rGradU)=0;

    /** Compute the contribution to the integral from one element */
    double CalculateOnElement(Element<ELEMENT_DIM,SPACE_DIM>& rElement)
    {
        double result_on_element = 0;
        
        GaussianQuadratureRule<ELEMENT_DIM> quad_rule(2);
        
        // \todo This assumes that the Jacobian is constant on an element.
        // This is true for linear basis functions, but not for any other type of basis function.
        const c_matrix<double, SPACE_DIM, SPACE_DIM> *p_inverse_jacobian = rElement.GetInverseJacobian();
        double jacobian_determinant = rElement.GetJacobianDeterminant();
        
        const unsigned num_nodes = rElement.GetNumNodes();
        
        // loop over Gauss points
        for (unsigned quad_index=0; quad_index < quad_rule.GetNumQuadPoints(); quad_index++)
        {
            const ChastePoint<ELEMENT_DIM>& quad_point = quad_rule.rGetQuadPoint(quad_index);
            
            c_vector<double, ELEMENT_DIM+1> phi;
            LinearBasisFunction<ELEMENT_DIM>::ComputeBasisFunctions(quad_point, phi);
            c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> grad_phi;
            LinearBasisFunction<ELEMENT_DIM>::ComputeTransformedBasisFunctionDerivatives(quad_point, *p_inverse_jacobian, grad_phi);
            
            // Location of the gauss point in the original element will be stored in x
            ChastePoint<SPACE_DIM> x(0,0,0);
            c_vector<double,PROBLEM_DIM> u = zero_vector<double>(PROBLEM_DIM);
            c_matrix<double,PROBLEM_DIM,SPACE_DIM> grad_u = zero_matrix<double>(PROBLEM_DIM,SPACE_DIM);
            
            for (unsigned i=0; i<num_nodes; i++)
            {
                const c_vector<double, SPACE_DIM>& r_node_loc = rElement.GetNode(i)->rGetLocation();
                
                // interpolate x
                x.rGetLocation() += phi(i)*r_node_loc;
                
                // interpolate u and grad u
                unsigned node_global_index = rElement.GetNodeGlobalIndex(i);
                for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
                {
                    // NOTE - following assumes that, if say there are two unknowns u and v, they
                    // are stored in the current solution vector as
                    // [U1 V1 U2 V2 ... U_n V_n]
                    unsigned index_into_vec = PROBLEM_DIM*node_global_index + index_of_unknown;
                    
                    double u_at_node = mSolutionReplicated[index_into_vec];
                    u(index_of_unknown) += phi(i)*u_at_node;
                    for (unsigned j=0; j<SPACE_DIM; j++)
                    {
                        grad_u(index_of_unknown,j) += grad_phi(j,i)*u_at_node;
                    }
                }
            }
            
            double wJ = jacobian_determinant * quad_rule.GetWeight(quad_index);
            result_on_element += GetIntegrand(x, u, grad_u) * wJ;
        }
        
        return result_on_element;
    }


public:
    virtual ~AbstractFunctionalCalculator()
    {
    }

    /** Calculate the integral over the given mesh, using the given solution
     *  vector on the mesh.
     */
    double Calculate(ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh,
                     Vec solution)
    {
        assert(solution);
        mSolutionReplicated.ReplicatePetscVector(solution);
        if (mSolutionReplicated.size() != rMesh.GetNumNodes() * PROBLEM_DIM)
        {
            EXCEPTION("The solution size does not match the mesh");
        }
        
        double result = 0;
    
        for (typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator
               iter = rMesh.GetElementIteratorBegin();
             iter != rMesh.GetElementIteratorEnd();
             ++iter)
        {
            result += CalculateOnElement(**iter);
        }
        
        return result;
    }    
};

#endif /*ABSTRACTFUNCTIONALCALCULATOR_HPP_*/
