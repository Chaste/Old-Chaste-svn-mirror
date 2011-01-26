"""Copyright (C) University of Oxford, 2005-2011

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
"""

"""
This file contains various classes supporting modifications to CellML models.
"""

import pycml
from pycml import *
import validator


class ModelModificationError(ValueError):
    """Error thrown if a model modification is invalid."""
    pass

    
class ModelModifier(object):
    """Base class supporting common model modification functionality.
    
    This class contains the logic to deal with adding/deleting variables, components, equations, etc.
    and connecting things up.  It also handles re-analysing the model when modifications have been
    completed to ensure that PyCml's internal data structures are up-to-date.
    
    Instances should be created with the model to modify as a single parameter.  Once all
    modifications have been completed, the finalize method must be called to ensure later
    processing of the model (e.g. code generation) will succeed.
    """
    def __init__(self, model):
        """Constructor."""
        self.model = model
        
    def reanalyse_model(self, error_handler):
        """Re-do the model validation steps needed for further processing of the model.
        
        Checks connections, etc. and builds up the dependency graph again, then performs
        a topological sort.
        
        If any errors are found during re-validation, the error_handler will be called with the
        list.  Warnings are ignored.
        
        TODO: figure out how to determine how much of this is actually needed - InterfaceGenerator
        can probably get away with less work.
        """
        # We want to see any errors
        logging_info = validator.CellMLValidator.setup_logging(show_errors=True, show_warnings=False)
        # Re-run validation & analysis
        self.model._check_variable_mappings()
        if not self.model._cml_validation_errors:
            assignment_exprs = self.model.search_for_assignments()
            self.model._check_assigned_vars(assignment_exprs)
        if not self.model._cml_validation_errors:
            self.model._classify_variables(assignment_exprs)
            self.model._order_variables(assignment_exprs)
        if not self.model._cml_validation_errors:
            warn_on_units_errors = False
            self.model._check_dimensional_consistency(assignment_exprs,
                                                      False,
                                                      warn_on_units_errors,
                                                      False)
        if self.model._cml_validation_errors:
            error_handler(self.model._cml_validation_errors)
        # Clear up logging
        validator.CellMLValidator.cleanup_logging(logging_info)
        
    def create_new_component(self, cname):
        """Create a new component in the model, ensuring the name is unique.
        
        If a component with name cname already exists,
        underscores will be added to the component name to make it unique.
        """
        while True:
            try:
                comp = self.model.get_component_by_name(cname)
                cname += u'_'
            except KeyError:
                # Component with this name doesn't exist
                break
        # Create the component
        comp = cellml_component.create_new(self.model, cname)
        self.model._add_component(comp)
        return comp
    
    def connect_variables(self, source, target):
        """Create a connection between the given source and target variables.
        
        The variables are both specified by a pair (cname,vname).  The source
        variable must exist within the model, whereas the target might not, in
        which case it will be created.
        
        The variable names must currently be identical.  TODO: This will probably have to change.
        """
        src_cname, src_vname = source
        target_cname, target_vname = target
        assert target_vname == src_vname
        src_comp = self.model.get_component_by_name(src_cname)
        target_comp = self.model.get_component_by_name(target_cname)
        # Determine encapsulation paths from target & source to the root
        src_path = self._parent_path(src_comp)
        target_path = self._parent_path(target_comp)
        # At some point these will share a common path, even if it's just the root itself
        meeting_index = self._find_common_tail(src_path, target_path)
        # Construct path from source to target, leaving out the root (None)
        path = src_path[:meeting_index]
        if src_path[meeting_index]:
            path.append(src_path[meeting_index])
        path.extend(reversed(target_path[:meeting_index]))
        # Traverse this path, adding connections at each step
        for i, src_comp in enumerate(path[:-1]):
            target_comp = path[i+1]
            self._make_connection(src_comp, target_comp, src_vname)
    
    def _make_connection(self, src_comp, target_comp, vname):
        """Make a connection between two components using the given variable name."""
        src_var = src_comp.get_variable_by_name(vname)
        target_var = self._find_or_create_variable(target_comp.name, vname, src_var)
        # Sanity check the target variable
        if target_var.get_type() == VarTypes.Mapped:
            # It must already be mapped to src_var; we're done
            assert (target_var.get_source_variable() == src_var or
                    src_var.get_source_variable() == target_var)
#            print "Connection exists between", src_var, "and target", target_var
            return
        elif target_var.get_type() == VarTypes.Unknown:
            # We've created this variable, so should be ok, but check for gotchas
            assert not(hasattr(target_var, u'initial_value'))
            if src_comp == target_comp.parent():
                src_if = u'private'
                target_if = u'public'
            elif src_comp.parent() == target_comp:
                src_if = u'public'
                target_if = u'private'
            else:
                assert src_comp.parent() == target_comp.parent()
                src_if = u'public'
                target_if = u'public'
                # One special case: if the src_var is actually obtained from a different
                # component at this level or above, in which case we should use the real
                # source, not that given.
                if getattr(src_var, src_if + u'_interface', u'none') == u'in':
                    src_var = src_var.get_source_variable()
            # Check and set the interface attributes
#            print "Connecting source", src_var, src_if, getattr(src_var, src_if + u'_interface', u'none'),
#            print "to", target_var, target_if, getattr(target_var, target_if + u'_interface', u'none')
            assert getattr(src_var, src_if + u'_interface', u'none') != u'in'
            assert getattr(target_var, target_if + u'_interface', u'none') != u'out'
            src_var.xml_set_attribute((src_if + u'_interface', None), u'out')
            target_var.xml_set_attribute((target_if + u'_interface', None), u'in')
            # Create the connection element
            self._create_connection_element(src_var, target_var)
            # Ensure we handle a later connection attempt between these variables correctly
            target_var._set_source_variable(src_var)
        else:
            # Ouch!  The model has used the same variable name for different things...
            raise ModelModificationError("Cannot connect " + target_var.fullname() + " to source " +
                                         src_var.fullname() + " as the target has the wrong type.")
    
    def _create_connection_element(self, var1, var2):
        """Create a connection element connecting the given variables and add to the model.
        
        If there's already a connection element for the relevant pair of components,
        we just add another map_variables element to that.
        """
        xpath_template = u'(cml:map_components/@component_1 = "%s" and cml:map_components/@component_2 = "%s")'
        cn1, cn2 = var1.component.name, var2.component.name
        xpath_predicate = xpath_template % (cn1, cn2) + u' or ' + xpath_template % (cn2, cn1)
        conn = self.model.xml_xpath(u'cml:connection[%s]' % (xpath_predicate,))
        if conn:
            conn = conn[0]
            # Check if we need to swap the variables
            if conn.map_components.component_1 == cn2:
                var1, var2 = var2, var1
        else:
            conn = var1.xml_create_element(u'connection', NSS[u'cml'])
            mapc = var1.xml_create_element(u'map_components', NSS[u'cml'],
                                           attributes={u'component_1': var1.component.name,
                                                       u'component_2': var2.component.name})
            conn.xml_append(mapc)
            self.model.xml_append(conn)
        mapv = var1.xml_create_element(u'map_variables', NSS[u'cml'],
                                       attributes={u'variable_1': var1.name,
                                                   u'variable_2': var2.name})
        conn.xml_append(mapv)
    
    def _find_common_tail(self, l1, l2):
        """Find the first element at which both lists are identical from then on."""
        i = -1
        try:
            while l1[i] == l2[i]:
                i -= 1
        except IndexError:
            # One list is the tail of the other
            pass
        # i now gives the last differing element
        assert i < -1
        return i+1
        
    def _parent_path(self, comp):
        """Return a path of components from that given to the encapsulation root.
        
        The root is specified by None, since we're really dealing with a forest,
        not a tree.
        """
        path = [comp]
        while comp:
            path.append(comp.parent())
            comp = comp.parent()
        return path
    
    def _find_or_create_variable(self, cname, vname, source):
        """Find the given variable in the model, creating it if necessary.
        
        The variable will become a mapped variable with the given source.
        Hence if it is created it will have the same units.
        """
        try:
            var = self.model.get_variable_by_name(cname, vname)
        except KeyError:
            # Create it and add to model
            comp = self.model.get_component_by_name(cname)
            var = cellml_variable.create_new(comp, vname, source.units) # TODO: what if units are defined locally to source's component?
            comp._add_variable(var)
        return var


class InterfaceGenerator(ModelModifier):
    """Class for generating an interface between a CellML model and external code.
    
    This contains functionality for users to describe the interface desired by the external code, i.e.
    which variables are inputs and/or outputs, and expected units.  It will then create a new component
    within the CellML model containing these variables, and add units conversions where required.  The
    external code then only needs to interact with this new component.
    """
    pass
