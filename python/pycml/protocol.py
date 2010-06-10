"""Copyright (C) University of Oxford, 2005-2010

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
Defines the Protocol class, which encapsulates the input & output of a
simulation protocol.
"""

import pycml
from pycml import *
import validator

class ProtocolError(ValueError):
    """Error thrown if a Protocol instance is invalid."""
    pass

class Protocol(object):
    """A class representing a simulation protocol.

     * When a protocol is initialised, it should be passed a cellml_model, and
       alter the model equations to reflect the protocol: remove ODEs for
       state variables that are set by the protocol, replace equations for
       protocol inputs with forms specified by the protocol.  Then redo the
       topological sort to ensure the equations are ordered correctly?
     * The main translator classes should then not need to do anything special
       when using a protocol, except for considering the case where there are
       no ODEs.
     * Classes need a ComputeDerivedQuantities method, which takes an (optional)
       state variable vector, and time, and computes all quantities in the model
       which are non-state-var outputs.  The protocol system, possibly via
       model annotations, will need to be able to specify what goes in here.  We
       may need the ComputedVariable functionality in OdeSystemInformation to
       name entries in the result vector.  The pe:keep behaviour for computed
       vars could change to store variables in this vector?  Although we don't
       necessarily want all ionic currents in here - the main reason they get
       annotated is so that we can generate GetIIonic.  So I think we need
       separate annotations for specifying parameters and computed vars (which
       both imply pe:keep).
     * GetIIonic may need a time parameter, since protocol inputs may depend on
       time, and the generated code wouldn't compile otherwise.
    
    A fully initialised protocol contains the following attributes:
     * inputs - a list of protocol inputs.  These may be cellml_variable instances,
       to (re)define a variable in the model, or mathml_apply instances, to add or
       modify an equation.  Once modify_model has been called, these objects will
       also exist in the model.
    """
    def __init__(self, model, multi_stage=False):
        """Create a new protocol.
        
        Eventually this will have arguments to parse a protocol definition file
        and create the protocol.  For now, however, the only option is to specify
        multi_stage as True and set up the internal data structures yourself.
        Then call self.modify_model.
        """
        self._protocol_component = None
        self.model = model
        self.inputs = []
        self.outputs = []
        if not multi_stage:
            raise NotImplemented

    def modify_model(self):
        """Actually apply protocol modifications to the model.
        
        Prior to this being called, all variable references within self.inputs must
        use full names, i.e. 'component_name,variable_name'.  Variables without a
        component part will be placed into a new 'protocol' component.
        
        This method will add the items from self.inputs into the model, replacing
        variables with the same name, and equations that assign to the same variable.
        This may involve changing a variable's type from State to Computed, or vice
        versa.
        
        After the call, all names and name references will be 'local'.  Connections
        will be created between components as needed, and units definitions added, to
        ensure a valid model.  In order for this to work, if a variable has units that
        do not already exist in the model, the object *must* have an attribute
        _cml_units referring to a suitable cellml_units instance.
        
        Finally, the protocol outputs will be used to prune the model's assignments
        list so only assignments of interest are used to generate code.
        TODO: Check interaction of this with PE.
        """
        # Add units before variables before maths so the order of inputs doesn't matter so much.
        for input in filter(lambda i: isinstance(i, cellml_units), self.inputs):
            self._add_units_to_model(input)
        for input in filter(lambda i: isinstance(i, cellml_variable), self.inputs):
            self._add_variable_to_model(input)
        for input in filter(lambda i: isinstance(i, mathml_apply), self.inputs):
            self._add_maths_to_model(input)
        self._fix_model_connections()
        self._clear_model_caches()
        self._reanalyse_model()
        self._filter_assignments()
        
    def _reanalyse_model(self):
        """Re-do the model validation steps needed for further processing of the model.
        
        Checks connections, etc. and builds up the dependency graph again, then performs
        a topological sort.
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
            raise ProtocolError("Applying protocol created an invalid model.")
        # Clear up logging
        validator.CellMLValidator.cleanup_logging(logging_info)
    
    def _fix_model_connections(self):
        """Ensure the modified model has all the necessary connections between variables.
        
        Check mathematics for ci elements that refer to variables not defined in that
        component.  These must refer to the variable by its full name (i.e. 'cname,vname').
        These variables will be renamed to use local names by this method, which will
        also create local variables mapped to the relevant source variable if needed.
        
        This needs to take account of the fact that a variable in one nested component
        may need to be connected to a variable in another nested component, and so create
        variables in the parent components to connect the whole thing up.
        """
        for expr in self.model.search_for_assignments():
            for ci_elt in self._find_ci_elts(expr):
                vname = unicode(ci_elt)
                if u',' in vname:
                    cname, vname = self._split_name(vname)
                    comp = expr.component
                    if comp.name != cname:
                        self._connect_variables((cname,vname), (comp.name,vname))
                    # Now just rename to be local
                    ci_elt._rename(vname)
        
    def _connect_variables(self, source, target):
        """Create a connection between the given source and target variables.
        
        The variables are both specified by a pair (cname,vname).  The source
        variable must exist within the model, whereas the target might not, in
        which case it will be created.
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
            #print "Connection exists between", src_var, "and target", target_var
            return
        elif target_var.get_type() == VarTypes.Unknown:
            # We've created this variable, so should be ok, but check for gotchas
            assert not(hasattr(target_var, u'initial_value'))
            public_iface = getattr(target_var, u'public_interface', u'none')
            private_iface = getattr(target_var, u'private_interface', u'none')
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
            #print "Connecting source", src_var, src_if, getattr(src_var, src_if + u'_interface', u'none'),
            #print "to", target_var, target_if, getattr(target_var, target_if + u'_interface', u'none')
            assert getattr(src_var, src_if + u'_interface', u'none') != u'in'
            assert getattr(target_var, target_if + u'_interface', u'none') != u'out'
            src_var.xml_set_attribute((src_if + u'_interface', None), u'out')
            target_var.xml_set_attribute((target_if + u'_interface', None), u'in')
            # Create the connection element
            self._create_connection_element(src_var, target_var)
            # Ensure we handle a latter connection attempt between these variables correctly
            target_var._set_source_variable(src_var)
        else:
            # Ouch!  The model has used the same variable name for different things...
            raise ProtocolError("Cannot connect " + target_var.fullname() + " to source " +
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
    
    def _get_protocol_component(self):
        """Get the protocol component in the model, creating it if necessary.
        
        New variables created just for use by the simulation protocol get put into a
        new 'protocol' component.  If a component with that name already exists,
        underscores will be added to the component name to make it unique.
        """
        if self._protocol_component is None:
            cname = u'protocol'
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
            self._protocol_component = comp
        return self._protocol_component
    
    def _split_name(self, full_name):
        """Split a full name into cname,vname, creating the component if needed.
        
        If the full_name doesn't contain a component part, the 'protocol' component
        will be used.
        """
        parts = full_name.split(',')
        if len(parts) == 2:
            cname, vname = parts
        elif len(parts) == 1:
            cname = self._get_protocol_component().name
            vname = full_name
        else:
            raise ValueError("Invalid variable name: " + full_name)
        return cname, vname
        
    def _rename_local_variables(self, expr):
        """
        Change local variable references in the given expression to refer
        explicitly to the protocol component.
        """
        for ci_elt in self._find_ci_elts(expr):
            vname = unicode(ci_elt)
            if u',' not in vname:
                # Do the rename
                cname = self._get_protocol_component().name
                full_name = cname + u',' + vname
                expr._rename(full_name)
    
    def _find_ci_elts(self, expr):
        """Get an iterator over all ci elements on the descendent-or-self axis of the given element."""
        if isinstance(expr, mathml_ci):
            yield expr
        elif hasattr(expr, 'xml_children'):
            # Recurse
            for child in expr.xml_children:
                for ci_elt in self._find_ci_elts(child):
                    yield ci_elt
    
    def _add_units_to_model(self, units):
        """Add a units definition to the model.
        
        For now, we just add to the model element without checking (TODO).
        """
        self.model.xml_append(units)
        self.model.add_units(units.name, units)
        
    def _add_variable_to_model(self, var):
        """Add or replace a variable in our model.
        
        We don't really do any checking for this case - just add the variable.
        This means that some 'possible' changes don't actually make sense, for
        instance giving a computed variable an initial value will trigger a later
        validation error, unless its definition is also changed to an ODE.
        (To change a variable to a constant, you need to replace its definition
        with a constant expression.)
        
        The one check that is made is that you don't change the interface
        definitions if replacing a variable, otherwise self._fix_model_connections
        could get confused.
        """
        cname, vname = self._split_name(var.name)
        comp = self.model.get_component_by_name(cname)
        try:
            orig_var = comp.get_variable_by_name(vname)
        except KeyError:
            orig_var = None
        if orig_var:
            # We're replacing a variable
            for iface in [u'public', u'private']:
                n = iface + u'_interface'
                assert hasattr(orig_var, n) == hasattr(var, n), "You are not allowed to change a variable's interfaces"
                if hasattr(var, n):
                    assert getattr(orig_var, n) == getattr(var, n), "You are not allowed to change a variable's interfaces"
            comp._del_variable(orig_var)
        var.name = vname
        comp._add_variable(var)

    def _add_maths_to_model(self, expr):
        """Add or replace an equation in the model.
        
        This case is more complex than variables, since we may need to change
        the type of variable assigned to, depending on the expression.
        
        Note: variable references within ci elements in the given expression
        should use full names (i.e. cname,vname).  Any local names will be
        assumed to refer to variables in the protocol component, and modified
        by self._rename_local_variables.  Later, self._fix_model_connections
        will change all references to use local names.
        """
        assert isinstance(expr, mathml_apply)
        assert expr.operator().localName == u'eq', 'Expression is not an assignment'
        # Figure out what's on the LHS of the assignment
        lhs = expr.operands().next()
        if lhs.localName == u'ci':
            # Straight assignment to variable
            cname, vname = self._split_name(unicode(lhs))
            assigned_var = self.model.get_variable_by_name(cname, vname)
            self._remove_existing_definition(assigned_var, False)
            self._add_expr_to_comp(cname, expr)
        else:
            # This had better be an ODE
            assert lhs.localName == u'apply', 'Expression is not a straight assignment or ODE'
            assert lhs.operator().localName == u'diff', 'Expression is not a straight assignment or ODE'
            dep_var = lhs.operands().next()
            assert dep_var.localName == u'ci', 'ODE is malformed'
            cname, dep_var_name = self._split_name(unicode(dep_var))
            dep_var = self.model.get_variable_by_name(cname, dep_var_name)
            self._remove_existing_definition(dep_var, True)
            self._add_expr_to_comp(cname, expr)

    def _remove_existing_definition(self, var, keep_initial_value):
        """Remove any existing definition (as an equation) of the given variable.
        
        If keep_initial_value is False, then also remove any initial_value attribute.
        
        If the variable is Mapped, throw a ProtocolError.
        """
        if var.get_type() == VarTypes.Mapped:
            raise ProtocolError("Cannot add new mathematics defining a mapped variable - change the definition of its source instead")
        if not keep_initial_value and hasattr(var, u'initial_value'):
            del var.initial_value
        # Note: if this is a variable added by the protocol, then it shouldn't have
        # any dependencies set up yet, so this is a no-op.
        for dep in var._get_all_expr_dependencies():
            assert isinstance(dep, mathml_apply)
            dep.xml_parent.xml_remove_child(dep)
            dep.xml_parent = None # Not done by Amara...
    
    def _add_expr_to_comp(self, cname, expr):
        """Add an expression to the mathematics in the given component."""
        comp = self.model.get_component_by_name(cname)
        if not hasattr(comp, u'math'):
            # Create the math element
            math = comp.xml_create_element(u'math', NSS[u'm'])
            comp.xml_append(math)
        # Append this expression
        comp.math.xml_append(expr)
    
    def _filter_assignments(self):
        """Apply protocol outputs to reduce the model size.
        
        Protocol outputs are a list of variable objects that are of interest.
        The assignments used in computing the model should be filtered so that
        only those needed for determining the outputs are used.  This has the
        potential to greatly simplify the model simulation.
        
        The only assignments should be to output variables, or nodes required
        in computing these.  Note that if one of these nodes is a state variable,
        we also require its derivative and the dependencies thereof.
        If there are no outputs specified, we leave the list unchanged.
        
        In addition, any output computed variable should be annotated as a 
        derived quantity, and any output constant annotated as a parameter, to
        ensure they are available for querying.  Other variables should have these
        annotations (and pe:keep) removed.
        """
        # Remove existing annotations
        for var in self.model.get_all_variables():
            var.set_pe_keep(False)
            var.set_is_derived_quantity(False)
            var.set_is_modifiable_parameter(False)
        # Add annotations for outputs
        for var in self.outputs:
            assert isinstance(var, cellml_variable)
            if var.get_type() == VarTypes.Constant:
                var.set_is_modifiable_parameter(True)
            elif var.get_type() == VarTypes.Computed:
                var.set_is_derived_quantity(True)
        # Filter assignments list (slightly hacky)
        if self.outputs:
            needed_nodes = self.model.calculate_extended_dependencies(self.outputs,
                                                                      state_vars_depend_on_odes=True)
            new_assignments = filter(lambda node: node in needed_nodes,
                                     self.model.get_assignments())
            self.model._cml_assignments = new_assignments
