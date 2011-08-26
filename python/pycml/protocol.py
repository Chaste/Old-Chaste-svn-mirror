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
    """A class representing part of a simulation protocol for the functional curation system.

    PyCml is responsible for implementing the 'model interface' section of protocols.
    See https://chaste.cs.ox.ac.uk/cgi-bin/trac.cgi/wiki/SimulationProtocolNotes#Modelinterface
    for further details.  The main user-facing methods here implement the XML elements defined
    there:
     - specify_input_variable
     - specify_output_variable
     - set_independent_variable_units
     - declare_new_variable
     - add_or_replace_equation
     - define_units_conversion_rule
    """
    def __init__(self, model, multi_stage=False, namespaces={}):
        """Create a new protocol.
        
        Eventually this will have arguments to parse a protocol definition file
        and create the protocol.  For now, however, the only option is to specify
        multi_stage as True and set up the internal data structures yourself.
        Then call self.modify_model.
        """
        self._protocol_component = None
        self._units_converter = None
        self._pending_oxmeta_assignments = []
        self._free_var_has_changed = None
        self.set_protocol_namespaces(namespaces)
        self.model = model
        self.inputs = set()
        self.outputs = set()
        if not multi_stage:
            raise NotImplementedError
    
    def specify_output_variable(self, prefixed_name, units=None):
        """Set the given variable as a protocol output, optionally in the given units.
        
        The units must be added to the model if they don't exist.  If they differ from the
        variable's original units, a conversion will be needed, and hence a new version
        of the variable will be added to the protocol component, with a suitable connection.
        Otherwise, we can just record the existing variable as an output.
        
        If units are not given, then keep the variable in its original units.
        
        The actual output variable is returned, although code shouldn't need to use it.
        """
        try:
            var = self._lookup_ontology_term(prefixed_name)
        except ValueError:
            raise ProtocolError("There is no model variable annotated with the term " + prefixed_name)
        if units is None:
            units = var.get_units()
        new_var = self.specify_as_output(var, units)
        if var.get_type() is VarTypes.Free:
            self.set_independent_variable_units(units, new_var)
        return new_var
    
    def specify_input_variable(self, prefixed_name, units=None, initial_value=None):
        """Set the given variable as a protocol input, optionally in the given units.
        
        The variable will become a constant parameter that can be set from the protocol.  If units
        are given and differ from its original units, they must be added to the model if they don't
        exist, and a new version of the variable added to the protocol component, with a suitable
        assignment.
        
        If the initial_value is given then this will overwrite the original setting if present.
        Any assignment to the variable will be removed.
        
        If the variable does not exist, this is only an error if an initial_value or units are not given.
        Otherwise we just create the variable.
        
        TODO: does not units-convert the initial value.
        """
        try:
            var = self._lookup_ontology_term(prefixed_name)
        except ValueError:
            if initial_value is None or units is None:
                raise ProtocolError("There is no model variable annotated with the term " + prefixed_name)
            else:
                var = None
        if units is None:
            units = var.get_units()
        if var is None:
            oxmeta_name = prefixed_name.split(':')[1]
            var = self.add_variable(self._get_protocol_component(), oxmeta_name, units, id=oxmeta_name)
            var.set_oxmeta_name(oxmeta_name)
            self._pending_oxmeta_assignments.append((var, oxmeta_name))
        var = self.specify_as_input(var, units)
        if initial_value:
            var.initial_value = unicode(initial_value)
        return var
    
    def set_independent_variable_units(self, units, newVar=None):
        """Set the independent variable to occur in the given units.
        
        If newVar is given, then this is being called by self.specify_output_variable, since the
        independent variable is also a protocol output, and the new version has already been
        created.
        
        Since this may mean we're called twice, we ensure the second call is a no-op.
        """
        if self._free_var_has_changed:
            assert self._free_var_has_changed.get_units().equals(units)
            assert newVar is None or newVar is self._free_var_has_changed
        else:
            t = self.model.find_free_vars()[0]
            units = self._get_units_object(units)
            if not units.equals(t.get_units()):
                # We'll need a conversion, and to convert all ODEs too
                self._free_var_has_changed = newVar or self._replace_variable(t, units)
    
    def declare_new_variable(self, name, units, initial_value=None):
        """Declare a new variable for use in the model interface.
        
        The variable will be added to the protocol component, and the units added to the model
        if they're not already present.  An assertion will be tripped if a variable with the
        given name already exists.
        """
        var = self.add_variable(self._get_protocol_component(), name, units)
        if initial_value:
            var.initial_value = unicode(initial_value)
        return var
    
    def add_or_replace_equation(self, assignment):
        """Add the given assignment equation to the model.
        
        It will replace any existing definition of the same variable.
        """
        assert isinstance(assignment, mathml_apply)
        assert assignment.operator().localName == u'eq'
        self.inputs.add(assignment)
    
    def add_units_conversion_rule(self, from_units, to_units, conv_expr):
        """Add a new biology-aware units conversion rule.
        
        The third argument must be an instance of mathml_lambda taking a single parameter.
        It will effectively be passed the RHS of an assignment expression, which has units
        dimensionally equivalent to from_units, and must return an expression with units
        dimensionally equivalent to to_units.
        
        Note that the UnitsConverter class needs a function object that will modify the
        assignment equation in-place, so that's what we create and store.
        """
        converter = self.get_units_converter()
        children = list(conv_expr.xml_element_children())
        assert len(children) == 2, "A units conversion rule must have a single bound variable: " + conv_expr.xml()
        assert children[0].localName == u'bvar', "A units conversion rule must have a single bound variable: " + conv_expr.xml()
        bvar_name = unicode(children[0].ci).strip()
        body_expr = children[1]
        # Modify variable references within the body_expr so they use fully qualified names
        # (compname, varname), except for uses of the bound variable
        try:
            self._identify_referenced_variables(body_expr, bvar_name)
        except ValueError:
            # We can't apply this rule as some required variables are missing
            print "Warning: unable to utilise units conversion rule below as required variables missing;",
            print "this may lead to later units conversion errors:"
            print "  From", from_units.description(), "to", to_units.description(), "via", conv_expr.xml()
            return
        func = lambda assignment: self._apply_conversion_rule(assignment, body_expr, bvar_name)
        converter.add_special_conversion(from_units, to_units, func)
    
    def set_protocol_namespaces(self, mapping):
        """Set the prefix->URI mapping used by the protocol file."""
        self._protocol_namespaces = mapping

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
        if self._free_var_has_changed:
            self._split_all_odes()
        self._fix_model_connections()
        self.finalize(self._error_handler, self._add_units_conversions)
        self._filter_assignments()
        # This is a bit of a hack!
        for var, oxmeta_name in self._pending_oxmeta_assignments:
            var.set_oxmeta_name(oxmeta_name)
    
    def specify_as_output(self, var, units):
        """Specify the given variable within the model as a protocol output.
        
        The output is wanted in the given units, which must be added to the model if they don't exist.
        If they differ from its current units, a conversion will be needed, and hence a new version
        of this variable will be added to the protocol component, with a suitable connection.
        Otherwise, we can just record the existing variable as an output.
        
        TODO: We can't yet specify a variable as both an output and an input but in different units.
        """
        units = self._get_units_object(units)
        if units.equals(var.get_units()):
            output_var = var
        else:
            if var in self.inputs:
                raise ProtocolError("You can't specify a variable as output and input in different units!")
            output_var = self._replace_variable(var, units)
            self.connect_variables(var, output_var)
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
            input_var = self._replace_variable(var, units)
            input_var.initial_value = getattr(var, u'initial_value', u'0')
            self.del_attr(var, u'initial_value', None)
            # Set all variables connected to the original variable (including itself) to be mapped to the new one
            self._update_connections(var, input_var)
        input_var._set_type(VarTypes.Constant)
        input_var._cml_ok_as_input = True
        self.inputs.add(input_var)
        return input_var
        
    def _replace_variable(self, var, units, allow_existing=False):
        """Replace the given variable with a version in the given units in the protocol component.
        
        Ensures that the units are added to the model if needed, and transfers the cmeta:id if
        present.  It doesn't transfer the initial_value, since this would break the output variable
        case.
        
        The new variable will be given a local name equal to the full name of the original, to avoid
        potential conflicts.  If allow_existing is False then it's an error if the variable has
        already been replaced.  If allow_existing is True, then we just reuse the existing replacement.
        """
        units = self.add_units(units)
        new_name = var.fullname(cellml=True)
        comp = self._get_protocol_component()
        try:
            existing_replacement = comp.get_variable_by_name(new_name)
        except KeyError:
            existing_replacement = None
        if existing_replacement:
            if not allow_existing:
                raise ProtocolError("Variable '%s' has already been replaced!" % new_name)
            new_var = existing_replacement
        else:
            new_var = self.add_variable(comp, new_name, units, id=var.cmeta_id)
            self.del_attr(var, u'id', NSS['cmeta'])
        return new_var
    
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
    
    def _split_all_odes(self):
        """The free variable has been units-converted, so adjust all ODEs to account for this.
        
        We copy all state variables into the protocol component, assign their RHS to a new variable,
        and create a new derivative in the protocol component to which this is assigned.
        
        Also check all equations for occurrences of derivatives on the RHS, and change them to refer
        to the new variable instead.
        
        TODO: How to deal with ODEs added by the protocol?  Will they be picked up?
        TODO: Require explicit call to set free var units even if an output?
        """
        free_var = self._free_var_has_changed
        old_free_var = self.model.find_free_vars()[0]
        self._update_connections(old_free_var, free_var)
        deriv_rhs = {}
        comp = self._get_protocol_component()
        for old_var in self.model.find_state_vars():
            if old_var.component is not comp:
                new_var = self._replace_variable(old_var, old_var.get_units(), allow_existing=True)
                new_var.initial_value = old_var.initial_value
                del old_var.initial_value
                if old_var in self.outputs:
                    self.outputs.remove(old_var)
                    self.outputs.add(new_var)
                # Add a new variable to assign the RHS to, with units of the original derivative
                deriv_name = self._uniquify_var_name(u'd_%s_d_%s' % (old_var.name, free_var.name), old_var.component)
                orig_ode = old_var.get_all_expr_dependencies()[0]
                orig_rhs_var = self.add_variable(old_var.component, deriv_name, orig_ode.eq.lhs.get_units().extract())
                deriv_rhs[new_var] = orig_rhs_var
                # Add a version of this in the protocol component, with desired units
                desired_units = new_var.get_units().quotient(free_var.get_units())
                mapped_rhs_var = self._replace_variable(orig_rhs_var, desired_units)
                self.connect_variables(orig_rhs_var, mapped_rhs_var)
                # Replace the original ODE with an assignment
                orig_rhs = orig_ode.eq.rhs
                orig_ode.safe_remove_child(orig_rhs)
                self.remove_expr(orig_ode)
                self.add_expr_to_comp(old_var.component,
                                      mathml_apply.create_new(self.model, u'eq',
                                                              [orig_rhs_var.name, orig_rhs]))
                # Create a new ODE in the interface component
                new_ode = mathml_diff.create_new(self.model, free_var.name, new_var.name, mapped_rhs_var.name)
                self.add_expr_to_comp(new_var.component, new_ode)
                new_ode.classify_variables(root=True, dependencies_only=True)
                # Update connections to the state variable
                self._update_connections(old_var, new_var)
            else:
                print "Ignoring state var", old_var
                raise NotImplementedError # TODO
        # Transform references to derivatives
        def xform(expr):
            state_var = expr.diff.dependent_variable.get_source_variable(recurse=True)
            rhs_var = deriv_rhs[state_var]
            # Ensure there's something mapped to it in this component
            rhs_var = self.connect_variables(rhs_var, (expr.component.name, rhs_var.name))
            # Update this expression
            parent = expr.xml_parent
            parent.xml_insert_after(expr, mathml_ci.create_new(parent, rhs_var.name))
            parent.safe_remove_child(expr)
        for expr in self.model.search_for_assignments():
            self._process_operator(list(expr.operands())[1], u'diff', xform)

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
        """Split a full name into cname,vname, creating the protocol component if needed.
        
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
        """Change local variable references in the given expression to refer to the protocol component.
        
        TODO: not actually used!  Remove?
        """
        for ci_elt in self._find_ci_elts(expr):
            vname = unicode(ci_elt)
            if u',' not in vname:
                # Do the rename
                cname = self._get_protocol_component().name
                full_name = cname + u',' + vname
                ci_elt._rename(full_name)
    
    def _find_ci_elts(self, expr):
        """Get an iterator over all ci elements on the descendent-or-self axis of the given element."""
        if isinstance(expr, mathml_ci):
            yield expr
        elif hasattr(expr, 'xml_children'):
            # Recurse
            for child in expr.xml_children:
                for ci_elt in self._find_ci_elts(child):
                    yield ci_elt
    
    def _identify_referenced_variables(self, expr, special_name=None):
        """Figure out which variables are referenced in the given expression, and update ci elements.
        
        The expression should contain names as used in the protocol, i.e. prefixed names giving an
        ontology name for model variables, and bare names for variables added by the protocol.
        Change each ci element to use the full "compname,varname" format.
        A ValueError is raised if any referenced variable doesn't exist.
        However, any reference to special_name is not checked and left as-is.
        Also, names already in fully qualified form are assumed to be ok and left as-is.
        """
        for ci_elt in self._find_ci_elts(expr):
            vname = unicode(ci_elt)
            if vname == special_name:
                continue
            if ',' in vname:
                continue
            if ':' in vname:
                var = self._lookup_ontology_term(vname)
            else:
                try:
                    var = self._get_protocol_component().get_variable_by_name(vname)
                except KeyError:
                    raise ValueError("The variable name '%s' has not been declared in the protocol"
                                     % vname)
            full_name = var.component.name + u',' + var.name
            ci_elt._rename(full_name)
    
    def _lookup_ontology_term(self, prefixed_name):
        """Find the variable annotated with the given term, if it exists.
        
        The term should be given in prefixed form, with the prefix appearing in the protocol's namespace
        mapping (prefix->uri, as found e.g. at elt.rootNode.xmlns_prefixes).
        
        Currently we just support the oxmeta annotations.
        
        Will throw ValueError if the variable doesn't exist in the model, or the given term is invalid.
        """
        try:
            prefix, localname = prefixed_name.split(':')
        except ValueError:
            raise ValueError("Ontology term '%s' is not a qname - it doesn't have a namespace prefix"
                             % prefixed_name)
        try:
            nsuri = self._protocol_namespaces[prefix]
        except KeyError:
            raise ValueError("The namespace prefix '%s' has not been declared" % prefix)
        if nsuri != NSS['oxmeta']:
            raise ValueError("We only support 'oxmeta' annotations at present")
        return self.model.get_variable_by_oxmeta_name(localname)
    
    def _apply_conversion_rule(self, assignment, conv_template, placeholder_name):
        """Apply a units conversion rule defined by self.add_units_conversion_rule.
        
        Modify the given assignment in-place, replacing the RHS by a copy of conv_template, except
        ci references to placeholder_name are replaced by (a copy of) the original RHS.
        
        TODO: check connection creation for used vars
        """
        rhs = assignment.eq.rhs
        assignment.safe_remove_child(rhs)
        new_rhs = mathml.clone(conv_template)
        copy_rhs = False
        for ci_elt in self._find_ci_elts(new_rhs):
            vname = unicode(ci_elt).strip()
            if vname == placeholder_name:
                # Copy the original RHS here, except if it's the first use don't bother copying
                if copy_rhs:
                    rhs = mathml.clone(rhs)
                else:
                    copy_rhs = True
                ci_elt.xml_parent.replace_child(ci_elt, rhs)
            else:
                # Ensure we have connections needed to get the variable in this component
                cname, local_name = vname.split(',')
                our_cname = assignment.component.name
                if cname != our_cname:
                    local_var = self.connect_variables((cname, local_name), (our_cname, local_name))
                    local_name = local_var.name
                ci_elt._rename(local_name)
        assignment.xml_append(new_rhs)

    def _add_units_conversions(self):
        """Apply units conversions, in particular 'special' ones, to the protocol component.
        
        Also convert any equations that have been added to the model, even if they don't appear
        in the protocol component, since otherwise we might not do all necessary conversions
        between model & protocol mathematics.
        """
        converter = self.get_units_converter()
        notifier = NotifyHandler(level=logging.WARNING)
        logging.getLogger('units-converter').addHandler(notifier)
        proto_comp = self._get_protocol_component()
        converter.add_conversions_for_component(proto_comp)
        converter.convert_assignments(filter(lambda i: isinstance(i, mathml_apply), self.inputs))
        converter.finalize(self._error_handler, check_units=False)
        notifier.flush()
        logging.getLogger('units-converter').removeHandler(notifier)
        if notifier.messages:
            raise ProtocolError("Unable to apply units conversions to the model/protocol interface")
        
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
        should use full names (i.e. cname,vname) or ontology terms.  Any local names will be
        assumed to refer to variables in the protocol component, and modified
        by self._rename_local_variables.  Later, self._fix_model_connections
        will change all references to use local names.
        """
        assert isinstance(expr, mathml_apply)
        assert expr.operator().localName == u'eq', 'Expression is not an assignment'
        self._identify_referenced_variables(expr)
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
    
    New protocols should be written in the pure XML syntax, which we parse and build up
    the protocol object ourselves.  However, legacy protocols may be Python code with a
    method apply_protocol(doc) to do the donkey work.
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
        proto_xml = amara_parse_cellml(proto_file_path)
        assert hasattr(proto_xml, u'protocol')
        proto = Protocol(doc.model, multi_stage=True, namespaces=proto_xml.xmlns_prefixes)
        proto_units = {}
        # TODO: Tidy up the units code so that we can refer to standard units without getting those in
        # the model too!
        if hasattr(proto_xml.protocol, u'units'):
            # Parse units definitions
            for defn in getattr(proto_xml.protocol.units, u'units', []):
                proto_units[defn.name] = defn
                defn.xml_parent = doc.model
        def get_units(elt, attr='units'):
            if hasattr(elt, attr):
                uname = getattr(elt, attr)
                try:
                    return proto_units[uname]
                except KeyError:
                    raise ProtocolError("Units '%s' have not been defined in the protocol" % uname)
            else:
                return None
        if hasattr(proto_xml.protocol, u'modelInterface'):
            if hasattr(proto_xml.protocol.modelInterface, u'setIndependentVariableUnits'):
                proto.set_independent_variable_units(get_units(proto_xml.protocol.modelInterface.setIndependentVariableUnits))
            for input in getattr(proto_xml.protocol.modelInterface, u'specifyInputVariable', []):
                proto.specify_input_variable(input.name, get_units(input), getattr(input, u'initial_value', None))
            for output in getattr(proto_xml.protocol.modelInterface, u'specifyOutputVariable', []):
                proto.specify_output_variable(output.name, get_units(output))
            for vardecl in getattr(proto_xml.protocol.modelInterface, u'declareNewVariable', []):
                proto.declare_new_variable(vardecl.name, get_units(vardecl), getattr(vardecl, u'initial_value', None))
            for expr in getattr(proto_xml.protocol.modelInterface, u'addOrReplaceEquation', []):
                proto.add_or_replace_equation(expr.xml_element_children().next())
            for rule in getattr(proto_xml.protocol.modelInterface, u'unitsConversionRule', []):
                proto.add_units_conversion_rule(get_units(rule, 'actualDimensions'),
                                                get_units(rule, 'desiredDimensions'),
                                                rule.xml_element_children().next())
        proto.modify_model()
    else:
        raise ProtocolError("Unexpected protocol file extension for file: " + proto_file_path)
