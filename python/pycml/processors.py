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
        
    def finalize(self, error_handler):
        """Re-do the model validation steps needed for further processing of the model.
        
        Checks connections, etc. and builds up the dependency graph again, then performs
        a topological sort.
        
        If any errors are found during re-validation, the error_handler will be called with the
        list.  Warnings are ignored.
        
        TODO: figure out how to determine how much of this is actually needed - InterfaceGenerator
        can probably get away with less work.
        """
        self._clear_model_caches()
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
            self.model._check_dimensional_consistency(assignment_exprs,
                                                      xml_context=False,
                                                      warn_on_units_errors=self.model.get_option('warn_on_units_errors'),
                                                      check_for_units_conversions=False)
        if self.model._cml_validation_errors:
            error_handler(self.model._cml_validation_errors)
        # Clear up logging
        validator.CellMLValidator.cleanup_logging(logging_info)

    def _clear_model_caches(self):
        """
        Clear cached links in the model, since we'll need to recompute many of them
        once we've finished modifying it.  Also clears dependency information.
        """
        for comp in getattr(self.model, u'component', []):
            for math in getattr(comp, u'math', []):
                math._unset_cached_links()
        for var in self.model.get_all_variables():
            var.clear_dependency_info()
        assignment_exprs = self.model.search_for_assignments()
        for expr in assignment_exprs:
            expr.clear_dependency_info()
    
    def create_new_component(self, cname):
        """Create a new component in the model, ensuring the name is unique.
        
        If a component with name cname already exists,
        underscores will be added to the component name to make it unique.
        """
        while True:
            try:
                self.model.get_component_by_name(cname)
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
        
        The variables are both specified either by a pair (cname,vname), or as cellml_variable objects.
        The source variable must exist within the model, whereas the target might not, in
        which case it will be created.
        
        Note that in the case that both source and target exist, it might NOT be the case that
        target obtains its value from source.  They might already be connected, and source obtains
        its value from target.  Or they might both obtain their value from a common source.
        
        The variable names must currently be identical.  TODO: This will probably have to change.
        """
#        print "Connecting", target, "to source", source
        if isinstance(source, cellml_variable):
            src_cname, src_vname = source.component.name, source.name
        else:
            src_cname, src_vname = source
        if isinstance(target, cellml_variable):
            target_cname, target_vname = target.component.name, target.name
        else:
            target_cname, target_vname = target
        assert target_vname == src_vname
        src_comp = self.model.get_component_by_name(src_cname)
        target_comp = self.model.get_component_by_name(target_cname)
        if src_comp == target_comp:
            return
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
        """Make a connection between two components using the given variable name.
        
        Note that in the case that both variables already exist and are connected, the existing
        connection is allowed to flow in either direction.
        """
        src_var = src_comp.get_variable_by_name(vname)
        target_var = self._find_or_create_variable(target_comp.name, vname, src_var)
        # Sanity check the target variable
        if target_var.get_type() == VarTypes.Mapped:
            # It must already be mapped to src_var; we're done
            assert (target_var.get_source_variable() is src_var
                    or src_var.get_source_variable() is target_var)
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
    
    def _find_connection_element(self, var1, var2):
        """Find any connection element containing a connection of the given variables.
        
        Returns a pair, the first element of which is either the element or None, and the
        second of which is a boolean indicating whether the variables need to be swapped
        in order to match the order of the components in the connection.
        """
        xpath_template = u'(cml:map_components/@component_1 = "%s" and cml:map_components/@component_2 = "%s")'
        cn1, cn2 = var1.component.name, var2.component.name
        xpath_predicate = xpath_template % (cn1, cn2) + u' or ' + xpath_template % (cn2, cn1)
        conn = self.model.xml_xpath(u'cml:connection[%s]' % (xpath_predicate,))
        if conn:
            conn = conn[0]
            swap = conn.map_components.component_1 == cn2
        else:
            conn = None
            swap = False
        return conn, swap
    
    def _create_connection_element(self, var1, var2):
        """Create a connection element connecting the given variables and add to the model.
        
        If there's already a connection element for the relevant pair of components,
        we just add another map_variables element to that.
        """
        conn, swap = self._find_connection_element(var1, var2)
        if conn:
            if swap:
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
    
    def remove_connection(self, var1, var2):
        """Remove a connection between two variables.
        
        Removes the relevant map_variables element.
        If this results in an empty connection element, removes that as well.
        """
        conn, swap = self._find_connection_element(var1, var2)
        if not conn:
            raise ModelModificationError("Cannot remove non-existent connection.")
        if swap:
            var1, var2 = var2, var1
        # Find the relevant map_variables element
        mapv = conn.xml_xpath(u'cml:map_variables[@variable_1="%s" and @variable_2="%s"]'
                              % (var1.name, var2.name))
        if not mapv:
            raise ModelModificationError("Cannot remove non-existent connection.")
        conn.xml_remove_child(mapv[0])
        if not hasattr(conn, u'map_variables'):
            conn.xml_parent.xml_remove_child(conn)
    
    def remove_connections(self, var):
        """Remove all connection elements for the given variable.
        
        Removes each relevant map_variables element.
        If this results in an empty connection element, removes that as well.
        """
#        print "Removing all connections to", var
        cname, vname = var.component.name, var.name
        for conn in list(getattr(self.model, u'connection', [])):
            if cname == conn.map_components.component_1:
                vid = u'variable_1'
            elif cname == conn.map_components.component_2:
                vid = u'variable_2'
            else:
                continue
            for mapv in conn.map_variables:
                if vname == getattr(mapv, vid, ''):
                    # Found a connection
                    conn.xml_remove_child(mapv)
                    if not hasattr(conn, u'map_variables'):
                        conn.xml_parent.xml_remove_child(conn)
                    # There can't be any more matching map_variables in this connection
                    break
    
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

    def _process_operator(self, expr, operator, func, *args, **kwargs):
        """Apply func to any application of the given operator within the given tree."""
        for elt in self.model.xml_element_children(expr):
            self._process_operator(elt, operator, func, *args, **kwargs)
        if isinstance(expr, mathml_apply) and expr.operator().localName == operator:
            func(expr, *args, **kwargs)

    def _find_or_create_variable(self, cname, vname, source):
        """Find the given variable in the model, creating it if necessary.
        
        The variable will become a mapped variable with the given source.
        Hence if it is created it will have the same units.
        """
        try:
            var = self.model.get_variable_by_name(cname, vname)
        except KeyError:
            # Create it and add to model
            units = source.component.get_units_by_name(source.units)
            var = self.add_variable(cname, vname, units)
        return var
    
    def add_variable(self, comp, vname, units, **kwargs):
        """Add a new variable to the given component.
        
        Remaining arguments are as for cellml_variable.create_new.
        Returns the new variable object.
        """
        if not isinstance(comp, cellml_component):
            comp = self.model.get_component_by_name(comp)
        units = self.add_units(units)
        var = cellml_variable.create_new(comp, vname, units.name, **kwargs)
        comp._add_variable(var)
        return var
    
    def add_units(self, units):
        """Add a units definition to the model, if it doesn't already exist.
        
        If the definition isn't in the model, at whole-model level, it will be added.  If the same
        definition is already available, however, that definition should be used by preference.
        Will return the actual units object to use.
        """
#        print "add_units(",repr(units),")",
        units = self.model._get_units_obj(units)
#        print "now",repr(units),
        try:
            model_units = self.model.get_units_by_name(units.name)
        except KeyError:
            model_units = None
        if model_units:
#            print "in model", repr(model_units)
            units = model_units # TODO: Check they're the same definition!
        else:
#            print "adding with name", units.name
            self.model.add_units(units.name, units)
            self.model.xml_append(units)
            # Ensure referenced units exist
            for unit in getattr(units, u'unit', []):
                unit._set_units_element(self.add_units(unit.get_units_element()), override=True)
                unit.units = unit.get_units_element().name
        return units
    
    def add_expr_to_comp(self, comp, expr):
        """Add an expression to the mathematics in the given component.
        
        comp may be a cellml_component instance or a component name.
        """
        if not isinstance(comp, cellml_component):
            comp = self.model.get_component_by_name(comp)
        if not hasattr(comp, u'math'):
            # Create the math element
            math = comp.xml_create_element(u'math', NSS[u'm'])
            comp.xml_append(math)
        # Append this expression
        comp.math.xml_append(expr)
    
    def remove_expr(self, expr):
        """Remove an expression (ODE or assignment) from its parent."""
        assert isinstance(expr, mathml_apply)
        expr.xml_parent.safe_remove_child(expr)
        expr.xml_parent = None # Not done by Amara...
        return expr
    
    def remove_definition(self, var, keep_initial_value=False):
        """Remove any existing definition (as an equation) of the given variable.
        
        If keep_initial_value is False, then also remove any initial_value attribute.
        
        If the variable is Mapped, throw a ModelModificationError.
        """
        if var.get_type() == VarTypes.Mapped:
            raise ModelModificationError("Cannot remove the equation defining a mapped variable - remove the definition of its source instead")
        if not keep_initial_value:
            self.del_attr(var, u'initial_value')
        # Note: if this is a variable added by a protocol, then it shouldn't have
        # any dependencies set up yet, so this is a no-op.
        for dep in var.get_all_expr_dependencies():
            self.remove_expr(dep)
        # We know don't know how it will be defined
        var._set_type(VarTypes.Unknown)

    def del_attr(self, elt, localName, ns=None):
        """Delete an XML attribute from an element, if it exists."""
        for (pyname, (qname, ns_)) in elt.xml_attributes.items():
            _, name = SplitQName(qname)
            if ns_ == ns and name == localName:
                delattr(elt, pyname)


class InterfaceGenerator(ModelModifier):
    """Class for generating an interface between a CellML model and external code.
    
    This contains functionality for users to describe the interface desired by the external code, i.e.
    which variables are inputs and/or outputs, and expected units.  It will then create a new component
    within the CellML model containing these variables, and add units conversions where required.  The
    external code then only needs to interact with this new component.
    """
    def __init__(self, model, name='interface'):
        super(InterfaceGenerator, self).__init__(model)
        self._interface_component = None
        self._interface_component_name = name

    def add_input(self, var, units, annotate=True):
        """Specify a variable as an input to the model.
        
        var should be a cellml_variable object already existing in the model.
        units should be a suitable input to self._get_units_object.
        
        If adding both State and Free variables as inputs, make sure to add the Free variable first,
        otherwise you will not be able to specify units for it.
        
        Set annotate to False if you do not wish a Constant variable to be annotated as a modifiable
        parameter.
        
        The new variable added to the interface component is returned.
        """
        assert isinstance(var, cellml_variable)
        units = self._get_units_object(units)
        var = var.get_source_variable(recurse=True) # Ensure we work with source variables only
        var_name = var.name # TODO: May provide this as an optional input instead
        # Check that the variable has a suitable type to be an input
        t = var.get_type()
        if t == VarTypes.Computed:
            raise ModelModificationError("Cannot specify computed variable " + var.fullname()
                                         + " as an input")
        elif t not in [VarTypes.Constant, VarTypes.Free, VarTypes.State]:
            raise ModelModificationError("Variable " + var.fullname() + " has unexpected type " + str(t))
        # Add a new variable with desired units to the interface component
        comp = self._get_interface_component()
        newvar = self.add_variable(comp, var_name, units, id=var.cmeta_id,
                                   initial_value=getattr(var, u'initial_value', None),
                                   interfaces={u'public': u'out'})
        newvar._set_type(t)
        # Remove initial value and id from the original, if they exist
        self.del_attr(var, u'initial_value')
        self.del_attr(var, u'id', NSS['cmeta'])
        # If the original variable was a state variable, split the defining equation
        if t == VarTypes.State:
            self._split_ode(newvar, var)
        # Annotate the new variable as a parameter if the original was a constant
        if t == VarTypes.Constant and annotate:
            newvar.set_is_modifiable_parameter(True)

        # Set all variables connected to the original variable to be mapped to the new one
        vars = [v for v in self.model.get_all_variables() if v.get_source_variable(True) is var]
        # Remove old connections, including interfaces and types so creating the new connection works
        for v in vars:
            self.remove_connections(v)
            self.del_attr(v, u'public_interface')
            self.del_attr(v, u'private_interface')
            v.clear_dependency_info()
        # Create new connections
        for v in vars:
            self.connect_variables(newvar, v)
        return newvar
    
    def add_output(self, var, units, annotate=True):
        """Specify a variable as an output of the model.
        
        var should be a cellml_variable object already existing in the model.
        units should be a suitable input to self._get_units_object.
        If annotate is set to True, the new variable will be annotated as a derived quantity.
        
        TODO: Figure out what to do with metadata annotations so users get the expected behaviour.
        
        The new variable added to the interface component is returned.
        """
        assert isinstance(var, cellml_variable)
        units = self._get_units_object(units)
        var = var.get_source_variable(recurse=True)
        var_name = var.name
        comp = self._get_interface_component()
        newvar = self.add_variable(comp, var_name, units)
        self.connect_variables(var, newvar)
        if annotate:
            newvar.set_is_derived_quantity(True)
        return newvar
    
    def add_output_function(self, resultName, operator, argVars, units):
        """Add an output that's defined as a (MathML) function of existing model variables.
        
        The desired units are those of the function's result.  The function arguments will be
        imported with their units as given by the model, and the function calculated.  This result
        will then be units-converted if necessary.
        
        The new variable added to the interface component is returned.
        """
        # Add the result variable
        comp = self._get_interface_component()
        units = self._get_units_object(units)
        result_var = self.add_variable(comp, resultName, units)
        result_var.set_pe_keep(True)
        # Map the argument variables
        operands = []
        for var in argVars:
            operands.append(self.add_output(var, var.get_units(), annotate=False).name)
        # Create the new function and assign it to result_var
        expr = mathml_apply.create_new(self.model, operator, operands)
        assign = mathml_apply.create_new(self.model, u'eq', [result_var.name, expr])
        self.add_expr_to_comp(comp, assign)
        return result_var
    
    def make_var_constant(self, var, value):
        """Turn a variable into a constant."""
        self.remove_definition(var)
        var.clear_dependency_info()
        var.initial_value = unicode(str(value))
        var._set_type(VarTypes.Constant)
    
    def make_var_computed_constant(self, var, value):
        """Turn a variable into a Computed variable with constant value definition."""
        self.remove_definition(var)
        var.clear_dependency_info()
        defn = mathml_apply.create_new(self.model, u'eq', 
                                       [var.name, (unicode(str(value)), var.get_units().name)])
        self.add_expr_to_comp(var.component, defn)
        var._set_type(VarTypes.Computed)

    def finalize(self, *args, **kwargs):
        """Override finalize to also set up standard interface elements not defined individually."""
        self._add_all_odes_to_interface()
        self._transform_derivatives_on_rhs()
        super(InterfaceGenerator, self).finalize(*args, **kwargs)

    def _transform_derivatives_on_rhs(self):
        """Transform any equations with derivatives on the RHS to use the variable defining it instead.
        
        self._split_ode must have been used for all derivatives before calling this method.  This means
        that each ODE now has a variable to which the RHS is assigned.  Rather than using the derivative
        directly, which could break the dependency chain if units conversions are used for time, equations
        should refer to this new variable instead.
        """
        for expr in self.model.search_for_assignments():
            self._process_operator(list(expr.operands())[1], u'diff', self._transform_derivative_on_rhs)

    def _transform_derivative_on_rhs(self, expr):
        """Transform a derivative on the RHS of an equation to refer to the defining variable.
        
        Helper method used by self._transform_derivatives_on_rhs to do the actual transformation.
        """
        # Find the variable to use
        dep_var = expr.diff.dependent_variable.get_source_variable(recurse=True)
        indep_var = expr.diff.independent_variable.get_source_variable(recurse=True)
        ode = dep_var.get_ode_dependency(indep_var)
        rhs_var = ode.get_dependencies()[0].get_source_variable(recurse=True)
        # Ensure there's something mapped to it in this component
        self.connect_variables(rhs_var, (expr.component.name, rhs_var.name))
        # Update this expression
        parent = expr.xml_parent
        parent.xml_insert_after(expr, mathml_ci.create_new(parent, rhs_var.name))
        parent.safe_remove_child(expr)
    
    def _split_ode(self, newVar, oldVar):
        """Split an ODE definition so the derivative goes into the interface component.
        
        The RHS stays where it is, and is assigned to a new variable, which is connected to the interface
        component and assigned to the new derivative.  newVar is the new state variable in the interface
        component, and oldVar will soon be mapped to it by the caller.
        
        Any other equations in the model which use the derivative are transformed to use the new variable
        instead.
        """
        # Get the free variable in the interface component
        free_var = self.model.find_free_vars()[0]
        if free_var.component is not newVar.component:
            free_var = self.add_input(free_var, free_var.get_units())
        # Add a new variable to assign the RHS to, with units of the original derivative
        deriv_name = self._uniquify_var_name(u'd_%s_d_%s' % (oldVar.name, free_var.name), oldVar.component)
        orig_ode = oldVar.get_all_expr_dependencies()[0]
        orig_rhs_var = self.add_variable(oldVar.component, deriv_name, orig_ode.eq.lhs.get_units().extract())
        # Add an output version of this in the interface, with desired units
        desired_units = newVar.get_units().quotient(free_var.get_units())
        mapped_rhs_var = self.add_output(orig_rhs_var, desired_units, annotate=False)
        # Replace the original ODE with an assignment
        orig_rhs = orig_ode.eq.rhs
        orig_ode.safe_remove_child(orig_rhs)
        self.remove_expr(orig_ode)
        self.add_expr_to_comp(oldVar.component,
                              mathml_apply.create_new(self.model, u'eq',
                                                      [orig_rhs_var.name, orig_rhs]))
        # Create a new ODE in the interface component
        new_ode = mathml_diff.create_new(self.model, free_var.name, newVar.name, mapped_rhs_var.name)
        self.add_expr_to_comp(newVar.component, new_ode)
        new_ode.classify_variables(root=True, dependencies_only=True)
    
    def _uniquify_var_name(self, varname, comp):
        """Ensure varname is unique within the given component.
        
        Underscores will be appended to the name until it is unique.  The unique name will be returned.
        """
        while True:
            try:
                comp.get_variable_by_name(varname)
                varname += u'_'
            except:
                break
        return varname
    
    def _add_all_odes_to_interface(self):
        """All the derivatives should be considered as model outputs, and state variables as model inputs.
        
        For any that haven't been done explicitly, this method will add the corresponding state variable
        as an input, with its original units, which has the desired effect.
        """
        comp = self._get_interface_component()
        for var in self.model.find_state_vars():
            if var.component is not comp:
                self.add_input(var, var.get_units())
    
    def _get_interface_component(self):
        """Get the new component that will contain the interface.
        
        The name will be 'interface' component, unless a component with that name already exists,
        in which case underscores will be added to the component name to make it unique.
        """
        if self._interface_component is None:
            self._interface_component = self.create_new_component(unicode(self._interface_component_name))
            self.model.interface_component_name = unicode(self._interface_component_name)
            assert not self._interface_component.ignore_component_name
        return self._interface_component
    
    def _get_units_object(self, units):
        """Helper function to convert a units specification into a cellml_units object.
        
        The input can be a cellml_units object, in which case we just return it.
        However, it can also be a serialised CellML units definition, in which case it
        will be parsed to create the object.
        """
        if isinstance(units, cellml_units):
            # We're done
            pass
        else:
            units = amara_parse_cellml(unicode(units))
        assert isinstance(units, cellml_units)
        return units
