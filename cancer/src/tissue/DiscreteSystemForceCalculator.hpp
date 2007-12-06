#ifndef DISCRETESYSTEMFORCECALCULATOR_HPP_
#define DISCRETESYSTEMFORCECALCULATOR_HPP_

#include "Meineke2001SpringSystem.hpp"
#include "OutputFileHandler.hpp"


class DiscreteSystemForceCalculator
{
    friend class TestDiscreteSystemForceCalculator;
    
private:
    
    /**
     * Spring system, used to get mesh and calculate forces between nodes.
     */ 
    Meineke2001SpringSystem<2>& mrMeinekeSpringSystem;
    
    /**
     * Small parameter, used in GetSamplingAngles().
     */ 
    double mEpsilon; 
    
    /** The file that the results of CalculateExtremalNormalForces */ 
    out_stream mpStressResultsFile;    
    
    /**
     * Given a node index, returns the set of neighbouring node indices.
     */ 
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned index)
    {
        ConformingTetrahedralMesh<2,2>& r_mesh = mrMeinekeSpringSystem.GetTissue().rGetMesh();
        Node<2>* p_node = r_mesh.GetNode(index);
        
        std::set<unsigned> neighbouring_node_indices;
        
        for (Node<2>::ContainingElementIterator it = p_node->ContainingElementsBegin();
             it != p_node->ContainingElementsEnd();
             ++it)
        {
            Element<2,2>* p_element = r_mesh.GetElement(*it);
            for(unsigned i=0; i<p_element->GetNumNodes(); i++)
            {
                unsigned node_index = p_element->GetNodeGlobalIndex(i);
                if(node_index!=index)
                {
                    neighbouring_node_indices.insert(node_index);
                }
            }
        }
        return neighbouring_node_indices;
    }
        
    /**
     * Given a node index and angle of intersecting line, returns the tangential and normal forces.
     */ 
    std::vector<double> CalculateFtAndFn(unsigned index, double theta)
    {
        ConformingTetrahedralMesh<2,2>& r_mesh = mrMeinekeSpringSystem.GetTissue().rGetMesh();
        
        std::set<unsigned> neighbouring_node_indices = GetNeighbouringNodeIndices(index);
        
        double total_tangential_force = 0.0;
        double total_normal_force = 0.0;
        
        for(std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
            iter != neighbouring_node_indices.end();
            iter++)
        {
            double alpha = r_mesh.GetAngleBetweenNodes(index, *iter);

            if(sin(alpha-theta)>0)
            {
                c_vector<double,2> force_between_nodes = mrMeinekeSpringSystem.CalculateForceBetweenNodes(index,*iter);
                //std::cout << "force between nodes = " << force_between_nodes[0] << "," <<  force_between_nodes[1] << "\n" << std::flush;
                
                c_vector<double,2> unit_vec_between_nodes;
                unit_vec_between_nodes[0] = cos(alpha);
                unit_vec_between_nodes[1] = sin(alpha);
                //std::cout << "unit vector = " << unit_vec_between_nodes[0] << "," <<  unit_vec_between_nodes[1] << "\n" << std::flush;
                double plusminus_norm_force = inner_prod(force_between_nodes,unit_vec_between_nodes);         
                
                total_tangential_force += plusminus_norm_force * cos(alpha-theta);
                total_normal_force += plusminus_norm_force * sin(alpha-theta);
            }
        }

        std::vector<double> ret(2);
        ret[0] = total_tangential_force;
        ret[1] = total_normal_force;
        
        return ret;
    }
    
    
    /**
     * Given a node index, returns a vector of sampling angles which may be used
     * by GetExtremalAngles() to find the locations of local extrema of the normal force.
     */  
    std::vector<double> GetSamplingAngles(unsigned index)
    {
        ConformingTetrahedralMesh<2,2>& r_mesh = mrMeinekeSpringSystem.GetTissue().rGetMesh();
        std::set<unsigned> neighbouring_node_indices = GetNeighbouringNodeIndices(index);
        
        std::vector<double> sampling_angles(4*neighbouring_node_indices.size());
        
        unsigned i=0;
        for(std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
            iter != neighbouring_node_indices.end();
            iter++)
        {
            double alpha = r_mesh.GetAngleBetweenNodes(index, *iter);
            double alpha_minus_epsilon = alpha - mEpsilon;
            double alpha_plus_epsilon = alpha + mEpsilon;
            double alpha_minus_epsilon_plus_pi = alpha - mEpsilon + M_PI;
            double alpha_plus_epsilon_plus_pi = alpha + mEpsilon + M_PI;
            
            if (alpha_minus_epsilon > 2*M_PI)
            {
                sampling_angles[i] = alpha_minus_epsilon - 2*M_PI;
            }
            else if (alpha_minus_epsilon < 0)
            {
                sampling_angles[i] = alpha_minus_epsilon + 2*M_PI;
            }
            else
            {
                sampling_angles[i] = alpha_minus_epsilon;
            }
            assert(sampling_angles[i] < 2*M_PI);
            assert(sampling_angles[i] >= 0);
            i++;
            
            if (alpha_plus_epsilon > 2*M_PI)
            {
                sampling_angles[i] = alpha_plus_epsilon - 2*M_PI;
            }
            else if (alpha_plus_epsilon < 0)
            {
                sampling_angles[i] = alpha_plus_epsilon + 2*M_PI;
            }
            else
            {
                sampling_angles[i] = alpha_plus_epsilon;
            }
            assert(sampling_angles[i] < 2*M_PI);
            assert(sampling_angles[i] >= 0);
            i++;
            
            if (alpha_minus_epsilon_plus_pi > 2*M_PI)
            {
                sampling_angles[i] = alpha_minus_epsilon_plus_pi - 2*M_PI;
            }
            else if (alpha_minus_epsilon_plus_pi < 0)
            {
                sampling_angles[i] = alpha_minus_epsilon_plus_pi + 2*M_PI;
            }
            else
            {
                sampling_angles[i] = alpha_minus_epsilon_plus_pi;
            }
            assert(sampling_angles[i] < 2*M_PI);
            assert(sampling_angles[i] >= 0);
            i++;
            
            if (alpha_plus_epsilon_plus_pi > 2*M_PI)
            {
                sampling_angles[i] = alpha_plus_epsilon_plus_pi - 2*M_PI;
            }
            else if (alpha_plus_epsilon_plus_pi < 0)
            {
                sampling_angles[i] = alpha_plus_epsilon_plus_pi + 2*M_PI;
            }
            else
            {
                sampling_angles[i] = alpha_plus_epsilon_plus_pi;
            }
            assert(sampling_angles[i] < 2*M_PI);
            assert(sampling_angles[i] >= 0);
            i++;            
        }
        
        sort(sampling_angles.begin(), sampling_angles.end());
        
        return sampling_angles;
    }
    
    
    /**
     * Given a node index and two sampling angles, finds the location of the root 
     * of the tangential force in the interval between the two angles.
     */ 
    double GetLocalExtremum(unsigned index, double angle1, double angle2)
    {        
        assert(angle1 < angle2);
        
        double tolerance = 1e-5;
        unsigned counter = 0;
                
        double previous_angle;        
        double current_error;     
        double current_angle = angle1; 
        
        current_error = angle2 - angle1;            
        std::vector<double> current_ft_and_fn(2);
        
        while (current_error>tolerance)
        {
            previous_angle = current_angle;
            current_ft_and_fn = CalculateFtAndFn(index,current_angle);
            current_angle -= current_ft_and_fn[0]/current_ft_and_fn[1];
            current_error = fabs(current_angle - previous_angle);
            counter++;
        }
        
        assert( current_angle>angle1 && current_angle<angle2 );            
        assert( current_error < tolerance );        
        
        if (current_angle<0.0)
        {
            current_angle += 2*M_PI;
        }
        return current_angle;        
    }
    

    /**
     * Given a vector of sampling angles, returns a vector of extremal angles, i.e. angles 
     * at which local extrema of the normal force occur.
     */ 
    std::vector<double> GetExtremalAngles(unsigned index, std::vector<double> samplingAngles)
    {
        std::vector<double> extremal_angles;
        std::vector<double> ft_and_fn(2);
        std::vector<double> tangential_force(samplingAngles.size());
        
        for (unsigned i=0; i<samplingAngles.size(); i++)
        {
            ft_and_fn = CalculateFtAndFn(index, samplingAngles[i]);
            tangential_force[i] = ft_and_fn[0];
        }         
        
        for (unsigned i=0; i<samplingAngles.size()-1; i++)
        {
            if ( ( (tangential_force[i] > 0) && (tangential_force[i+1] < 0) ) ||  
                 ( (tangential_force[i] < 0) && (tangential_force[i+1] > 0) ) )
            {
                double next_extremal_angle;
                
                // if we jump through zero, then we know exactly where the local extremum is
                if (samplingAngles[i+1] - samplingAngles[i] < 2*mEpsilon + 1e-6 )
                {
                    next_extremal_angle = (samplingAngles[i+1] + samplingAngles[i])/2.0;
                    extremal_angles.push_back(next_extremal_angle);
                }
                else
                {
                    next_extremal_angle = GetLocalExtremum(index, samplingAngles[i], samplingAngles[i+1]);
                    extremal_angles.push_back(next_extremal_angle);
                }
            }
        }
        
        // separate bit of code for the interval from the last sampling angel back round to the first...
        
        samplingAngles[samplingAngles.size()-1] -= 2*M_PI;
        
        if ( ( (tangential_force[samplingAngles.size()-1] > 0) && (tangential_force[0] < 0) ) ||  
                 ( (tangential_force[samplingAngles.size()-1] < 0) && (tangential_force[0] > 0) ) )
        {
            double next_extremal_angle;
            
            // if we jump through zero, then we know exactly where the local extremum is
            if ( samplingAngles[0] - samplingAngles[samplingAngles.size()-1] < 2*mEpsilon + 1e-6 )
            {
                next_extremal_angle = (samplingAngles[0] + samplingAngles[samplingAngles.size()-1])/2.0;
                extremal_angles.push_back(next_extremal_angle);
            }
            else
            {
                next_extremal_angle = GetLocalExtremum(index, samplingAngles[samplingAngles.size()-1], samplingAngles[0]);
                extremal_angles.push_back(next_extremal_angle);
            }
        }
        
        return extremal_angles;        
    }
    
public:
    
    DiscreteSystemForceCalculator(Meineke2001SpringSystem<2>& rMeinekeSpringSystem)
        : mrMeinekeSpringSystem(rMeinekeSpringSystem),
          mEpsilon(0.01)
    {
    }
    
    std::vector< std::vector<double> > CalculateExtremalNormalForces()
    {
        std::vector< std::vector<double> > extremal_normal_forces;
        
        ConformingTetrahedralMesh<2,2>& r_mesh = mrMeinekeSpringSystem.GetTissue().rGetMesh();
        
        std::vector<double> minimum_normal_forces(r_mesh.GetNumNodes());
        std::vector<double> maximum_normal_forces(r_mesh.GetNumNodes());
                
        for (unsigned i=0; i<r_mesh.GetNumNodes(); i++)
        {
            std::vector<double> sampling_angles = GetSamplingAngles(i);            
            std::vector<double> extremal_angles = GetExtremalAngles(i,sampling_angles);

            double minimum_normal_force_for_node_i = DBL_MAX;
            double maximum_normal_force_for_node_i = -DBL_MAX;
            
            double current_normal_force;
            
            for (unsigned j=0; j<extremal_angles.size(); j++)
            {
                current_normal_force = CalculateFtAndFn(i,extremal_angles[j])[1];
                
                if (current_normal_force > maximum_normal_force_for_node_i)
                {
                    maximum_normal_force_for_node_i = current_normal_force;
                }
                
                if (current_normal_force < minimum_normal_force_for_node_i)
                {
                    minimum_normal_force_for_node_i = current_normal_force;
                }
            }   
            
            minimum_normal_forces[i] = minimum_normal_force_for_node_i;
            maximum_normal_forces[i] = maximum_normal_force_for_node_i;
            
        } 
                
        extremal_normal_forces.push_back(minimum_normal_forces);
        extremal_normal_forces.push_back(maximum_normal_forces);
        
        return extremal_normal_forces;
    }
    
    void WriteResultsToFile(std::string simulationOutputDirectory)
    {
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        double time = p_simulation_time->GetDimensionalisedTime();
        std::ostringstream time_string;
        time_string << time;
        std::string results_directory = simulationOutputDirectory +"/results_from_time_" + time_string.str();
            
        OutputFileHandler output_file_handler2(results_directory+"/vis_results/",false); 
        mpStressResultsFile = output_file_handler2.OpenOutputFile("results.vizstress");
        
        (*mpStressResultsFile) <<  time << "\t";
        
        double global_index;
        double x;
        double y;
        double minimum;
        double maximum;

        ConformingTetrahedralMesh<2,2>& r_mesh = mrMeinekeSpringSystem.GetTissue().rGetMesh();
        
        std::vector< std::vector<double> > extremal_normal_forces = CalculateExtremalNormalForces();
               
        for (unsigned i=0; i<r_mesh.GetNumNodes(); i++)
        {
            global_index = (double) i;
            
            x = r_mesh.GetNode(i)->rGetLocation()[0];
            y = r_mesh.GetNode(i)->rGetLocation()[1];
            
            minimum = extremal_normal_forces[0][i]; 
            maximum = extremal_normal_forces[1][i];
            
            (*mpStressResultsFile) << global_index << " " << x << " " << y << " " << minimum << " " << maximum << " ";
        }
        
        (*mpStressResultsFile) << "\n";        
        mpStressResultsFile->close();
    }

};
#endif /*DISCRETESYSTEMFORCECALCULATOR_HPP_*/
