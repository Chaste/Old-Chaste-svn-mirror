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
Defines the Protocol class, which encapsulates the input & output of a
simulation protocol.
"""

import pycml
from pycml import *
import processors
import validator

class ProtocolError(ValueError):
    """Error thrown if a Protocol instance is invalid."""
    pass

class Protocol(processors.ModelModifier):
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
     * inputs - an iterable of protocol inputs.  These may be cellml_variable instances,
       to (re)define a variable in the model, or mathml_apply instances, to add or
       modify an equation.  Once modify_model has been called, these objects will
       also exist in the model.
     * outputs - an iterable of protocol output variables.
    """
    def __init__(self, model, multi_stage=False):
        """Create a new protocol.
        
        Eventually this will have arguments to parse a protocol definition file
        and create the protocol.  For now, however, the only option is to specify
        multi_stage as True and set up the internal data structures yourself.
        Then call self.modify_model.
        """
        self._protocol_component = None
        self._units_converter = None
        self.model = model
        self.inputs = set()
        self.outputs = set()
        if not multi_stage:
            raise NotImplementedError

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
        """
        # Add units before variables before maths so the order of inputs doesn't matter so much.
        for input in filter(lambda i: isinstance(i, cellml_units), self.inputs):
            #self._check_input(input)
            self.add_units(input)
        for input in filter(lambda i: isinstance(i, cellml_variable), self.inputs):
            self._check_input(input)
            self._add_variable_to_model(input)
        for input in filter(lambda i: isinstance(i, mathml_apply), self.inputs):
            self._check_input(input)
            self._add_maths_to_model(input)
        self._fix_model_connections()
        self.finalize(self._error_handler, self._add_units_conversions)
        self._filter_assignments()
    
    def specify_as_output(self, var, units):
        """Specify the given variable within the model as a protocol output.
        
        The output is wanted in the given units, which must be added to the model if they don't exist.
        If they differ from its current units, a conversion will be needed, and hence a new version
        of this variable will be added to the protocol component, with a suitable connection.
        Otherwise, we can just record the existing variable as an output.
        """
        units = self._get_units_object(units)
        if units.equals(var.get_units()):
            output_var = var
        else:
            units = self.add_units(units)
            output_var = cellml_variable.create_new(var, var.name, units.name, id=var.cmeta_id)
            self.del_attr(var, u'id', NSS['cmeta'])
            self._get_protocol_component()._add_variable(output_var)
            self.connect_variables(var, output_var)
            self.inputs.add(output_var)
        self.outputs.add(output_var)
        return output_var
    
    def specify_as_input(self, var, units):
        """Specify the given variable within the model as a protocol input.
        
        The input is wanted in the given units, which must be added to the model if they don't exist.
        If they differ from its current units, a conversion will be needed, and hence a new version
        of this variable will be added to the protocol component, with a suitable assignment.
        
        TODO: does not units-convert the initial value.
        
        The variable that is the input must be set as a modifiable parameter, and any existing definition
        removed.
        """
        if var.get_type() == VarTypes.Mapped:
            raise ProtocolError("Cannot specify a mapped variable (%s) as an input." % var.fullname())
        # Remove any existing definition
        self.remove_definition(var, keep_initial_value=True)
        # Set up the input
        units = self._get_units_object(units)
        if units.equals(var.get_units()):
            input_var = var
            if not hasattr(var, u'initial_value'):
                var.initial_value = u'0'
        else:
            input_var = self.add_variable(self._get_protocol_component(), var.name, units, id=var.cmeta_id,
                                          initial_value=getattr(var, u'initial_value', u'0'))
            self.del_attr(var, u'id', NSS['cmeta'])
            self.del_attr(var, u'initial_value', None)
            # Set all variables connected to the original variable (including itself) to be mapped to the new one
            vars = [v for v in self.model.get_all_variables() if v.get_source_variable(True) is var]
            # Remove old connections, including interfaces and types so creating the new connection works
            for v in vars:
                self.remove_connections(v)
                self.del_attr(v, u'public_interface')
                self.del_attr(v, u'private_interface')
                v.clear_dependency_info()
            # Create new connections
            for v in vars:
                self.connect_variables(input_var, v)
        input_var._set_type(VarTypes.Constant)
        input_var._cml_ok_as_input = True
        self.inputs.add(input_var)
        return input_var
        
    def _check_input(self, input):
        """New inputs must not already exist in the model!"""
        if isinstance(input, cellml_units):
            exists = self.model.has_units(input)
        else:
            exists = self.model is getattr(input, 'xml_parent', None)
        if exists and not getattr(input, '_cml_ok_as_input', False):
            msg = "Inputs must not already exist in the model."
            msg += " (Input %s exists.)" % repr(input)
            raise ProtocolError(msg)
        
    def _error_handler(self, errors):
        """Deal with errors found when re-analysing a modified model."""
        raise ProtocolError("Applying protocol created an invalid model.")

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
                        self.connect_variables((cname,vname), (comp.name,vname))
                    # Now just rename to be local
                    ci_elt._rename(vname)
    
    def _get_protocol_component(self):
        """Get the protocol component in the model, creating it if necessary.
        
        New variables created just for use by the simulation protocol get put into a
        new 'protocol' component.  If a component with that name already exists,
        underscores will be added to the component name to make it unique.
        """
        if self._protocol_component is None:
            self._protocol_component = self.create_new_component(u'protocol')
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
    
    def _add_units_conversions(self):
        """Apply units conversions, in particular 'special' ones, to the protocol component.
        
        Also convert any equations that have been added to the model, even if they don't appear
        in the protocol component, since otherwise we might not do all necessary conversions
        between model & protocol mathematics.
        """
        converter = self.get_units_converter()
        proto_comp = self._get_protocol_component()
        converter.add_conversions_for_component(proto_comp)
        converter.convert_assignments(filter(lambda i: isinstance(i, mathml_apply), self.inputs))
        converter.finalize(self._error_handler, check_units=False)
        
    def get_units_converter(self):
        """Get the protocol's units converter object, in order to add 'special' conversions."""
        if not self._units_converter:
            warn_only = not self.model.get_option('fully_automatic') and self.model.get_option('warn_on_units_errors')
            self._units_converter = processors.UnitsConverter(self.model, warn_only)
        return self._units_converter

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
        if hasattr(var, 'xml_parent'):
            # It's already been added, e.g. by specify_as_input
            return
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
                assert getattr(orig_var, n, u'none') == getattr(var, n, u'none'), "You are not allowed to change a variable's interfaces"
            # Only keep RDF annotations if the cmeta:id is unchanged
            comp._del_variable(orig_var, keep_annotations=(orig_var.cmeta_id == var.cmeta_id))
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
            self.remove_definition(assigned_var, False)
            self.add_expr_to_comp(cname, expr)
        else:
            # This had better be an ODE
            assert lhs.localName == u'apply', 'Expression is not a straight assignment or ODE'
            assert lhs.operator().localName == u'diff', 'Expression is not a straight assignment or ODE'
            dep_var = lhs.operands().next()
            assert dep_var.localName == u'ci', 'ODE is malformed'
            cname, dep_var_name = self._split_name(unicode(dep_var))
            dep_var = self.model.get_variable_by_name(cname, dep_var_name)
            self.remove_definition(dep_var, True)
            self.add_expr_to_comp(cname, expr)
        
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
        if self.outputs:
            # Remove parts of the model that aren't needed
            needed_nodes = self.model.calculate_extended_dependencies(self.outputs,
                                                                      state_vars_depend_on_odes=True)
            needed_nodes.update([input for input in self.inputs
                                 if isinstance(input, (mathml_apply, cellml_variable))])
            for node in self.model.get_assignments()[:]:
                if node not in needed_nodes:
                    if isinstance(node, cellml_variable):
                        node.component._del_variable(node)
                    elif isinstance(node, mathml_apply):
                        node.xml_parent.xml_remove_child(node)
            # Update connection elements
            for conn in list(getattr(self.model, u'connection', [])):
                comp1 = self.model.get_component_by_name(conn.map_components.component_1)
                comp2 = self.model.get_component_by_name(conn.map_components.component_2)
                any_kept = False
                for mapv in list(conn.map_variables):
                    try:
                        var1 = comp1.get_variable_by_name(mapv.variable_1)
                        var2 = comp2.get_variable_by_name(mapv.variable_2)
                        any_kept = True
                    except KeyError:
                        # Remove connection
                        conn.xml_remove_child(mapv)
                if not any_kept:
                    self.model.xml_remove_child(conn)
            # Filter assignments list
            new_assignments = filter(lambda node: node in needed_nodes,
                                     self.model.get_assignments())
            self.model._cml_assignments = new_assignments
        # Remove existing annotations
        for var in self.model.get_all_variables():
            var.set_pe_keep(False)
            var.set_is_derived_quantity(False)
            var.set_is_modifiable_parameter(False)
        # Add annotations for inputs & outputs
        for var in [input for input in self.inputs if isinstance(input, cellml_variable)]:
            var.set_pe_keep(True)
            if var.get_type() == VarTypes.Constant:
                var.set_is_modifiable_parameter(True)
        for var in self.outputs:
            assert isinstance(var, cellml_variable)
            var.set_is_output_variable(True)
            if var.get_type() == VarTypes.Constant:
                var.set_is_modifiable_parameter(True)
            elif var.get_type() in [VarTypes.Computed, VarTypes.Mapped]:
                var.set_is_derived_quantity(True)
            else:
                assert var.get_type() in [VarTypes.State, VarTypes.Free]


def apply_protocol_file(doc, proto_file_path):
    """Apply the protocol defined in the given file to a model.
    
    Initially, this method is primarily to allow testing the protocol system.
    Hence we assume the protocol file is Python code which has a method
    apply_protocol(doc) to do the donkey work.
    """
    if proto_file_path[-3:] == '.py':
        import imp
        import os
        proto_dir = os.path.dirname(proto_file_path)
        proto_file_name = os.path.basename(proto_file_path)
        proto_module_name = os.path.splitext(proto_file_name)[0]
        (file, pathname, desc) = imp.find_module(proto_module_name, [proto_dir])
        try:
            proto = imp.load_module(proto_module_name, file, pathname, desc)
        finally:
            file.close()
        proto.apply_protocol(doc)
    elif proto_file_path[-4:] == '.xml':
        proto_xml = amara_parse(proto_file_path)
        assert hasattr(proto_xml, u'protocol')
        if hasattr(proto_xml.protocol, u'modelModification'):
            d = {'doc': doc,
                 'protocol': sys.modules[__name__],
                 '__builtins__': __builtins__}
            exec str(proto_xml.protocol.modelModification) in d
    else:
        raise ProtocolError("Unexpected protocol file extension for file: " + proto_file_path)
