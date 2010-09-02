#!/usr/bin/env python

# We want 1/2==0.5
from __future__ import division

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
This part of PyCml applies various optimising transformations to CellML
models, in particular partial evaluation and the use of lookup tables.
"""

# Common CellML processing stuff
import pycml
from pycml import *  # Put contents in the local namespace as well

__version__ = "$Revision$"[11:-2]



######################################################################
#                         Partial Evaluation                         #
######################################################################

class PartialEvaluator(object):
    """Perform partial evaluation of a CellML model."""
    def _debug(self, *args):
        """Output debug info from the PE process."""
        logger = logging.getLogger('partial-evaluator')
        logger.debug(' '.join(map(str, args)))

    def _expr_lhs(self, expr):
        """Display the LHS of this expression."""
        lhs = expr.assigned_variable()
        if isinstance(lhs, cellml_variable):
            return lhs.fullname()
        else:
            return lhs[0].fullname() + u'/' + lhs[1].fullname()

    def _rename_vars(self, elt):
        """Rename variables found in ci elements in this tree to use canonical names."""
        if isinstance(elt, mathml_ci):
            var = elt.variable.get_source_variable(recurse=True)
            elt.xml_remove_child(unicode(elt))
            elt.xml_append(var.fullname(cellml=True))
            elt._cml_variable = var
            self._debug("Using canonical name", unicode(elt))
        else:
            for e in self.doc.model.xml_element_children(elt):
                self._rename_vars(e)

    def parteval(self, doc):
        """Do the partial evaluation."""
        self.doc = doc
        # BTA
        doc.model.do_binding_time_analysis()
        # Reduce/evaluate expressions
        while True:
            doc.model._pe_repeat = u'no'
            exprs = [e for e in doc.model.get_assignments()
                     if isinstance(e, mathml_apply)]
            for expr in exprs:
                if expr._get_binding_time() is BINDING_TIMES.static:
                    value = expr.evaluate()
                    self._debug("Evaluated", self._expr_lhs(expr),
                               "to", value)
                    if expr.is_ode():
                        # Replace the RHS with a <cn> element giving the value
                        rhs = expr.eq.rhs
                        new_elt = expr._eval_self()
                        expr.xml_insert_after(rhs, new_elt)
                        expr.xml_remove_child(rhs)
                    else:
                        # The variable assigned to will have its initial_value set,
                        # so we don't need the expression any more
                        expr.xml_parent.xml_remove_child(expr)
                        doc.model._remove_assignment(expr)
                    # Update variable usage counts
                    expr._update_usage_counts(list(expr.operands())[1],
                                              remove=True)
                else:
                    expr._reduce()
            if doc.model._pe_repeat == u'no':
                break
            self._debug("----- looping -----")
        del doc.model._pe_repeat
        
        # Use canonical variable names on LHS of assignments
        for expr in [e for e in doc.model.get_assignments()
                     if isinstance(e, mathml_apply)]:
            # If the assigned-to variable isn't used or kept, remove the assignment
            if isinstance(expr.eq.lhs, mathml_ci):
                var = expr.eq.lhs.variable
                if not (var.get_usage_count() or var.pe_keep):
                    doc.model._remove_assignment(expr)
                    continue
            self._rename_vars(expr.eq.lhs)

        # Tidy up kept variables, in case they aren't referenced in an eq'n.
        for comp in getattr(doc.model, u'component', []):
            for var in getattr(comp, u'variable', []):
                if var.pe_keep:
                    var._reduce()

        # Collapse into a single component
        new_comp = cellml_component.create_new(doc, u'c')
        new_comp._cml_created_by_pe = True
        old_comps = list(getattr(doc.model, u'component', []))
        doc.model._add_component(new_comp)
        # We iterate over a copy of the component list so we can
        # delete components from the model in this loop, and so the
        # new component exists in the model so we can add content to
        # it.
        for comp in old_comps:
            # Move relevant contents into new_comp
            for units in list(getattr(comp, u'units', [])):
                # Copy all <units> elements
                # TODO: Just generate the ones we need,
                # using _ensure_units_exist
                comp.xml_remove_child(units)
                new_comp.xml_append(units)
            for var in list(getattr(comp, u'variable', [])):
                # Only move used source variables
                self._debug('Variable', var.fullname(), 'usage', var.get_usage_count(),
                           'type', var.get_type(), 'kept', var.pe_keep)
                if (var.get_usage_count() and
                    var.get_type() != VarTypes.Mapped) or var.pe_keep:
                    self._debug('Moving variable', var.fullname())
                    # Remove from where it was
                    comp._del_variable(var, keep_annotations=True)
                    # Set name to canonical version
                    var.name = var.fullname(cellml=True)
                    # Place in new component
                    new_comp._add_variable(var)
            # Don't copy reactions
            for math in list(getattr(comp, u'math', [])):
                # Copy all <math> elements with content
                if math.xml_children:
                    comp.xml_remove_child(math)
                    new_comp.xml_append(math)
                    # Invalidate cached links
                    math._unset_cached_links()
            doc.model._del_component(comp)
        # Remove groups & connections
        for group in list(getattr(doc.model, u'group', [])):
            doc.model.xml_remove_child(group)
        for conn in list(getattr(doc.model, u'connection', [])):
            doc.model.xml_remove_child(conn)

        # Remove unused variable assignments from the list
        vs = [v for v in doc.model.get_assignments()
              if isinstance(v, cellml_variable)]
        for v in vs:
            if not v.xml_parent is new_comp:
                doc.model._remove_assignment(v)

        # Remove interface attributes from variables
        for v in new_comp.variable:
            for iface in [u'public', u'private']:
                try:
                    delattr(v, iface+u'_interface')
                except AttributeError:
                    pass

        # Refresh expression dependency lists
        for expr in [e for e in doc.model.get_assignments()
                     if isinstance(e, mathml_apply)]:
            expr._cml_depends_on = list(expr.vars_in(expr.eq.rhs))
            if expr.is_ode():
                # Add dependency on the independent variable
                indep_var = expr.eq.lhs.diff.independent_variable
                if not indep_var in expr._cml_depends_on:
                    expr._cml_depends_on.append(indep_var)
        return


######################################################################
#                        Lookup table analysis                       #
######################################################################

class LookupTableAnalyser(object):
    """
    Analyses & annotates a CellML model to indicate where lookup
    tables can be used.
    """

    def __init__(self):
        """Create an analyser."""
        # No model to analyse yet
        self.doc = None
        # Set default parameter values
        self.set_params()

    def var_is_membrane_potential(self, var):
        """
        Determine if the given variable represents the
        transmembrane potential.

        This method takes an instance of cellml_variable and returns
        a boolean.
        """
        return (var.name in [u'V', u'membrane__V'] and
                var.get_type(follow_maps=True) == VarTypes.State)

    def is_allowed_variable(self, var):
        """Return True iff the given variable is allowed in a lookup table.

        This method uses the config store in the document to check the
        variable object.
        """
        return self.doc._cml_config.lut_config.has_key(
            var.get_source_variable(recurse=True))

    def is_keying_var(self, var):
        """Return True iff the given variable can be used as a table key.

        Will check the config store if it exists.  If not, the variable
        name must match self.table_var.
        """
        if hasattr(self.doc, '_cml_config'):
            return self.doc._cml_config.lut_config.has_key(
                var.get_source_variable(recurse=True))
        else:
            return var.name == self.table_var

    _LT_DEFAULTS = {'table_min': u'-100.0001',
                    'table_max': u'49.9999',
                    'table_step': u'0.01',
                    'table_var': u'V'}
    def set_params(self, **kw):
        """Set parameters controlling lookup table generation.

        Keyword parameters control the lookup table settings, which are
        stored as attributes on suitable expressions.
        table_min - minimum table entry (unicode) -> lut:min
        table_max - maximum table entry (unicode) -> lut:max
        table_step - table step size (unicode) -> lut:step
        table_var - the name of the variable indexing the table (unicode)
                  -> lut:var
        """
        defaults = self._LT_DEFAULTS
        for attr in defaults:
            if kw.has_key(attr):
                setattr(self, attr, kw[attr])
            else:
                setattr(self, attr, getattr(self, attr, defaults[attr]))
        return

    def get_param(self, param_name, table_var):
        """Get the value of the lookup table parameter.

        table_var is the variable object being used to key this table.

        If the document has a config store, lookup the value there.
        If that doesn't give us a value, use that given using set_params.
        """
        try:
            val = self.doc._cml_config.lut_config[
                table_var.get_source_variable(recurse=True)][param_name]
        except AttributeError, KeyError:
            val = getattr(self, param_name)
        return val

    # One of these functions is required for a lookup table to be
    # worthwhile
    lut_expensive_funcs = frozenset(('exp', 'log', 'ln', 'root',
                                     'sin', 'cos', 'tan',
                                     'sec', 'csc', 'cot',
                                     'sinh', 'cosh', 'tanh',
                                     'sech', 'csch', 'coth',
                                     'arcsin', 'arccos', 'arctan',
                                     'arcsinh', 'arccosh', 'arctanh',
                                     'arcsec', 'arccsc', 'arccot',
                                     'arcsech', 'arccsch', 'arccoth'))

    class LUTState(object):
        """Represents the state for lookup table analysis."""
        def __init__(self):
            """Set the initial state.

            We assume at first a lookup table would not be suitable.
            """
            self.has_var = False
            self.bad_vars = set()
            self.has_func = False
            self.table_var = None

        def update(self, res):
            """Update the state with the results of a recursive call.

            res is the result of checking a sub-expression for
            suitability, and should be another instance of this class.
            """
            self.has_var = (self.has_var or res.has_var) and \
                           (not (self.table_var and res.table_var) or
                            self.table_var.get_source_variable(recurse=True) is
                             res.table_var.get_source_variable(recurse=True))
            # The second condition above specifies that the keying variables
            # must be the same if they both exist
            self.bad_vars.update(res.bad_vars)
            self.has_func = self.has_func or res.has_func
            self.table_var = self.table_var or res.table_var

        def suitable(self):
            """
            Return True iff this state indicates a suitable expression
            for replacement with a lookup table.
            """
            return (self.has_var and
                    not self.bad_vars and
                    self.has_func)

        def reason(self):
            """
            Return a unicode string describing why this state indicates the
            expression is not suitable for replacement with a lookup table.

            This can be:
             'no_var' - doesn't contain the table variable
             'bad_var <vname>' - contains a variable which isn't permitted
             'no_func' - doesn't contain an expensive function
            or a comma separated combination of the above.
            """
            r = []
            if not self.has_var:
                r.append(u'no_var')
            if not self.has_func:
                r.append(u'no_func')
            for vname in self.bad_vars:
                r.append(u'bad_var ' + vname)
            return u','.join(r)

    def analyse_for_lut(self, expr, var_checker_fn):
        """Check if the given expression can be replaced by a lookup table.

        The first argument is the expression to check; the second is a
        function which takes a variable object and returns True iff this
        variable is permitted within a lookup table expression.

        If self.annotate_failures is True then annotate <apply> and
        <piecewise> expressions which don't qualify with the reason
        why they do not.
        This can be:
         'no_var' - doesn't contain the table variable
         'bad_var <vname>' - contains a variable which isn't permitted
         'no_func' - doesn't contain an expensive function
        or a comma separated combination of the above.
        The annotation is stored as the lut:reason attribute.

        If self.annotate_outermost_only is True then only annotate the
        outermost qualifying expression, rather than also annotating
        qualifying subexpressions.
        """
        # Initialise the indicators
        state = self.LUTState()
        # Process current node
        if isinstance(expr, mathml_ci):
            # Variable reference
            if var_checker_fn(expr.variable):
                # Could be a permitted var that isn't a keying var
                if self.is_keying_var(expr.variable):
                    state.has_var = True
                    state.table_var = expr.variable
            else:
                state.bad_vars.add(expr.variable.name)
        elif isinstance(expr, mathml_piecewise):
            # Recurse into pieces & otherwise options
            if hasattr(expr, u'otherwise'):
                r = self.analyse_for_lut(child_i(expr.otherwise, 1),
                                         var_checker_fn)
                state.update(r)
            for piece in getattr(expr, u'piece', []):
                r = self.analyse_for_lut(child_i(piece, 1), var_checker_fn)
                state.update(r)
                r = self.analyse_for_lut(child_i(piece, 2), var_checker_fn)
                state.update(r)
        elif isinstance(expr, mathml_apply):
            # Check function
            if (not state.has_func and
                expr.operator().localName in self.lut_expensive_funcs):
                state.has_func = True
            # Check operands
            for operand in expr.operands():
                r = self.analyse_for_lut(operand, var_checker_fn)
                state.update(r)
            # Check qualifiers
            for qual in expr.qualifiers():
                r = self.analyse_for_lut(qual, var_checker_fn)
                state.update(r)
        else:
            # Just recurse into children
            for e in expr.xml_children:
                if getattr(e, 'nodeType', None) == Node.ELEMENT_NODE:
                    r = self.analyse_for_lut(e, var_checker_fn)
                    state.update(r)
        # Annotate the expression if appropriate
        if isinstance(expr, (mathml_apply, mathml_piecewise)):
            if state.suitable():
                if self.annotate_outermost_only:
                    # Remove annotations from (expr and) child expressions
                    self.remove_lut_annotations(expr)
                for param in ['min', 'max', 'step']:
                    expr.xml_set_attribute((u'lut:' + param, NSS['lut']),
                                           self.get_param('table_' + param,
                                                          state.table_var))
                expr.xml_set_attribute((u'lut:var', NSS['lut']),
                                       state.table_var.name)
                expr.xml_set_attribute((u'lut:possible', NSS['lut']),
                                       u'yes')
                self.doc.lookup_tables[expr] = True
            else:
                if self.annotate_failures:
                    expr.xml_set_attribute((u'lut:reason', NSS['lut']),
                                           state.reason())
        return state

    def remove_lut_annotations(self, expr, remove_reason=False):
        """Remove lookup table annotations from the given expression.

        By default this will only remove annotations from expressions
        (and sub-expressions) that can be converted to use lookup tables.
        If remove_reason is True, then the lut:reason attributes on
        non-qualifying expressions will also be removed.
        """
        # Remove from this expression
        delete_table = False
        for pyname in getattr(expr, 'xml_attributes', {}).keys():
            fullname = expr.xml_attributes[pyname]
            if fullname[1] == NSS['lut']:
                if remove_reason or fullname[0] != u'lut:reason':
                    expr.__delattr__(pyname)
                    if fullname[0] != u'lut:reason':
                        delete_table = True
        # Delete expr from list of lookup tables?
        if delete_table:
            del self.doc.lookup_tables[expr]
        # Recurse into children
        for e in expr.xml_children:
            if getattr(e, 'nodeType', None) == Node.ELEMENT_NODE:
                self.remove_lut_annotations(e, remove_reason)

    def analyse_model(self, doc,
                      annotate_failures=True,
                      annotate_outermost_only=True):
        """Analyse the given document.

        This method checks all expressions (and subexpressions)
        in the given document for whether they could be converted to
        use a lookup table, and annotates them appropriately.

        By default expressions which don't qualify will be annotated
        to indicate why; set annotate_failures to False to suppress
        this.

        Also by default only the outermost suitable expression in any
        given tree will be annotated; if you want to annotate suitable
        subexpressions of a suitable expression then pass
        annotate_outermost_only as False.
        """
        self.doc = doc
        self.annotate_failures = annotate_failures
        self.annotate_outermost_only = annotate_outermost_only
        doc.lookup_tables = {}
        # How to check for allowed variables
        if hasattr(doc, '_cml_config'):
            checker_fn = self.is_allowed_variable
        else:
            checker_fn = self.var_is_membrane_potential

        # Check all expressions
        for expr in (e for e in doc.model.get_assignments()
                     if isinstance(e, mathml_apply)):
            ops = expr.operands()
            ops.next()
            e = ops.next()
            self.analyse_for_lut(e, checker_fn)

        # Assign names (numbers) to the lookup tables found.
        # Also work out which ones can share index variables into the
        # table.
        doc.lookup_tables = doc.lookup_tables.keys()
        doc.lookup_tables.sort(cmp=element_path_cmp)
        doc.lookup_table_indexes, n = {}, 0
        for i, expr in enumerate(doc.lookup_tables):
            expr.xml_set_attribute((u'lut:table_name', NSS['lut']),
                                   unicode(i))
            comp = expr.get_component()
            var = comp.get_variable_by_name(expr.var).get_source_variable(
                recurse=True)
            key = (expr.min, expr.max, expr.step, var)
            if not doc.lookup_table_indexes.has_key(key):
                doc.lookup_table_indexes[key] = unicode(n)
                n += 1
            expr.xml_set_attribute((u'lut:table_index', NSS['lut']),
                                   doc.lookup_table_indexes[key])
        
        # Re-do dependency analysis so that an expression using lookup
        # tables only depends on the keying variable.
        for expr in (e for e in doc.model.get_assignments()
                     if isinstance(e, mathml_apply)):
            expr.classify_variables(root=True,
                                    dependencies_only=True,
                                    needs_special_treatment=self.calculate_dependencies)

    def calculate_dependencies(self, expr):
        """Determine the dependencies of an expression that might use a lookup table.

        This method is suitable for use as the needs_special_treatment function in
        mathml_apply.classify_variables.  It is used to override the default recursion
        into sub-trees.  It takes a single sub-tree as argument, and returns either
        the dependency set for that sub-tree, or None to use the default recursion.
        
        Expressions that can use a lookup table only depend on the keying variable.
        """
        if expr.getAttributeNS(NSS['lut'], u'possible', '') == u'yes':
            key_var_name = expr.getAttributeNS(NSS['lut'], u'var')
            key_var = expr.component.get_variable_by_name(key_var_name).get_source_variable(recurse=True)
            return set([key_var])
        # If not a table, use default behaviour
        return None


######################################################################
#                          Jacobian analysis                         #
######################################################################

class LinearityAnalyser(object):
    """Analyse linearity aspects of a model.

    This class performs analyses to determine which ODEs have a linear
    dependence on their state variable, discounting the presence of
    the transmembrane potential.

    This can be used to decouple the ODEs for more efficient solution,
    especially in a multi-cellular context.

    analyse_for_jacobian(doc) must be called before
    rearrange_linear_odes(doc).
    """
    LINEAR_KINDS = Enum('None', 'Linear', 'Nonlinear')

    def analyse_for_jacobian(self, doc, V=None):
        """Analyse the model for computing a symbolic Jacobian.

        Determines automatically which variables will need to be solved
        for using Newton's method, and stores their names in
        doc.model._cml_nonlinear_system_variables, as a list of variable
        objects.

        Also stores doc.model._cml_linear_vars, doc.model._cml_free_var,
        doc.model._cml_transmembrane_potential.
        """
        # TODO: Add error checking and tidy
        stvs = doc.model.find_state_vars()
        if V is None:
            Vcname, Vvname = 'membrane', 'V'
            V = doc.model.get_variable_by_name(Vcname, Vvname)
        V = V.get_source_variable(recurse=True)
        doc.model._cml_transmembrane_potential = V
        free_var = doc.model.find_free_vars()[0]
        lvs = self.find_linear_odes(stvs, V, free_var)
        # Next 3 lines for benefit of rearrange_linear_odes(doc)
        lvs.sort(key=lambda v: v.fullname())
        doc.model._cml_linear_vars = lvs
        doc.model._cml_free_var = free_var
        # Store nonlinear vars in canonical order
        nonlinear_vars = list(set(stvs) - set([V]) - set(lvs))
        nonlinear_vars.sort(key=lambda v: v.fullname())
        doc.model._cml_nonlinear_system_variables = nonlinear_vars
        # Debugging
        f = lambda var: var.fullname()
        DEBUG('linearity-analyser', 'V=', V.fullname(), '; free var=',
              free_var.fullname(), '; linear vars=', map(f, lvs),
              '; nonlinear vars=', map(f, nonlinear_vars))
        return

    def _get_rhs(self, expr):
        """Return the RHS of an assignment expression."""
        ops = expr.operands()
        ops.next()
        return ops.next()

    def _check_expr(self, expr, state_var, bad_vars):
        """The actual linearity checking function.

        Recursively determine the type of dependence expr has on
        state_var.  The presence of any members of bad_vars indicates
        a non-linear dependency.

        Return a member of the self.LINEAR_KINDS enum.
        """
        kind = self.LINEAR_KINDS
        result = None
        if isinstance(expr, mathml_ci):
            var = expr.variable.get_source_variable(recurse=True)
            if var is state_var:
                result = kind.Linear
            elif var in bad_vars:
                result = kind.Nonlinear
            elif var.get_type(follow_maps=True) == VarTypes.Computed:
                # Recurse into defining expression
                src_var = var.get_source_variable(recurse=True)
                src_expr = self._get_rhs(src_var._get_dependencies()[0])
                DEBUG('find_linear_deps', "--recurse for", src_var.name,
                      "to", src_expr)
                result = self._check_expr(src_expr, state_var, bad_vars)
            else:
                result = kind.None
            # Record the kind of this variable, for later use when
            # rearranging linear ODEs
            var._cml_linear_kind = result
        elif isinstance(expr, mathml_cn):
            result = kind.None
        elif isinstance(expr, mathml_piecewise):
            # If any conditions have a dependence, then we're
            # nonlinear.  Otherwise, all the pieces must be the same
            # (and that's what we are) or we're nonlinear.
            pieces = getattr(expr, u'piece', [])
            conds = map(lambda p: child_i(p, 2), pieces)
            chld_exprs = map(lambda p: child_i(p, 1), pieces)
            if hasattr(expr, u'otherwise'):
                chld_exprs.append(child_i(expr.otherwise, 1))
            for cond in conds:
                if self._check_expr(cond, state_var, bad_vars) != kind.None:
                    result = kind.Nonlinear
                    break
            else:
                # Conditions all OK
                for e in chld_exprs:
                    res = self._check_expr(e, state_var, bad_vars)
                    if result is not None and res != result:
                        # We have a difference
                        result = kind.Nonlinear
                        break
                    result = res
        elif isinstance(expr, mathml_apply):
            # Behaviour depends on the operator
            operator = expr.operator().localName
            operands = expr.operands()
            if operator in ['plus', 'minus']:
                # Linear if any operand linear, and none non-linear
                op_kinds = map(lambda op: self._check_expr(op, state_var,
                                                           bad_vars),
                               operands)
                result = max(op_kinds)
            elif operator == 'divide':
                # Linear iff only numerator linear
                numer = operands.next()
                denom = operands.next()
                if self._check_expr(denom, state_var, bad_vars) != kind.None:
                    result = kind.Nonlinear
                else:
                    result = self._check_expr(numer, state_var, bad_vars)
            elif operator == 'times':
                # Linear iff only 1 linear operand
                op_kinds = map(lambda op: self._check_expr(op, state_var,
                                                           bad_vars),
                               operands)
                lin, nonlin = 0, 0
                for res in op_kinds:
                    if res == kind.Linear: lin += 1
                    elif res == kind.Nonlinear: nonlin += 1
                if nonlin > 0 or lin > 1:
                    result = kind.Nonlinear
                elif lin == 1:
                    result = kind.Linear
                else:
                    result = kind.None
            else:
                # Cannot be linear; may be no dependence at all
                result = max(map(lambda op: self._check_expr(op, state_var,
                                                             bad_vars),
                                 operands))
                if result == kind.Linear:
                    result = kind.Nonlinear
        else:
            # Either a straightforward container element
            try:
                child = child_i(expr, 1)
            except ValueError:
                # Assume it's just a constant
                result = kind.None
            else:
                result = self._check_expr(child, state_var, bad_vars)
        DEBUG('find_linear_deps', "Expression", expr, "gives result", result)
        return result

    def find_linear_odes(self, state_vars, V, free_var):
        """Identify linear ODEs.

        For each ODE (except that for V), determine whether it has a
        linear dependence on the dependent variable, and thus can be
        updated directly, without using Newton's method.

        We also require it to not depend on any other state variable,
        except for V.
        """
        kind = self.LINEAR_KINDS
        candidates = set(state_vars) - set([V])
        linear_vars = []
        for var in candidates:
            ode_expr = var._get_ode_dependency(free_var)
            if self._check_expr(self._get_rhs(ode_expr), var,
                                candidates - set([var])) == kind.Linear:
                linear_vars.append(var)
        return linear_vars

    def _clone(self, expr):
        """Properly clone a MathML sub-expression."""
        return mathml.clone(expr)

    def _make_apply(self, operator, ghs, i, filter_none=True,
                    preserve=False):
        """Construct a new apply expression for g or h.

        ghs is an iterable of (g,h) pairs for operands.
        
        i indicates whether to construct g (0) or h (1).
        
        filter_none indicates the behaviour of 0 under this operator.
        If True, it's an additive zero, otherwise it's a
        multiplicative zero.
        """
        # Find g or h operands
        ghs_i = map(lambda gh: gh[i], ghs)
        if not filter_none and None in ghs_i:
            # Whole expr is None
            new_expr = None
        else:
            # Filter out None subexprs
            if operator == u'minus':
                # Do we need to retain a unary minus?
                if len(ghs_i) == 1 or ghs_i[0] is None:
                    # Original was -a or 0-a
                    retain_unary_minus = True
                else:
                    # Original was a-0 or a-b
                    retain_unary_minus = False
            else:
                # Only retain if we're told to preserve as-is
                retain_unary_minus = preserve
            ghs_i = filter(None, ghs_i)
            if ghs_i:
                if len(ghs_i) > 1 or retain_unary_minus:
                    new_expr = mathml_apply.create_new(
                        self.__expr, operator, ghs_i)
                else:
                    new_expr = self._clone(ghs_i[0])
            else:
                new_expr = None
        return new_expr

    def _transfer_lut(self, expr, gh, var):
        """Transfer lookup table annotations from expr to gh.

        gh is a pair (g, h) s.t. expr = g + h*var.

        If expr can be represented by a lookup table, then the lookup
        variable cannot be var, since if it were, then expr would have
        a non-linear dependence on var.  Hence h must be 0, since
        otherwise expr would contain a (state) variable other than the
        lookup variable, and hence not be a candidate for a table.
        Thus expr=g, so we transfer the annotations to g.
        """
        if expr.getAttributeNS(NSS['lut'], u'possible', '') != u'yes':
            return
        # Paranoia check that our reasoning is correct
        g, h = gh
        assert h is None
        # Transfer the annotations into g
        for pyname, fullname in expr.xml_attributes.iteritems():
            if fullname[1] == NSS['lut']:
                g.xml_set_attribute(fullname, getattr(expr, pyname))
        # Make sure g has a reference to its component, for use by
        # code generation.
        g._cml_component = expr.component
        return

    def _rearrange_expr(self, expr, var):
        """Rearrange an expression into the form g + h*var.

        Performs a post-order traversal of this expression's tree,
        and returns a pair (g, h)
        """
        gh = None
        if isinstance(expr, mathml_ci):
            # Variable
            ci_var = expr.variable.get_source_variable(recurse=True)
            if var is ci_var:
                gh = (None, mathml_cn.create_new(expr,
                                                 u'1', u'dimensionless'))
            else:
                if ci_var._cml_linear_kind == self.LINEAR_KINDS.None:
                    # Just put the <ci> in g, but with full name
                    gh = (mathml_ci.create_new(expr, ci_var.fullname()), None)
                    gh[0]._cml_variable = ci_var
                else:
                    # ci_var is a linear function of var, so rearrange
                    # it's definition
                    if not hasattr(ci_var, '_cml_linear_split'):
                        ci_defn = ci_var._get_dependencies()[0]
                        ci_var._cml_linear_split = self._rearrange_expr(
                            self._get_rhs(ci_defn), var)
                    gh = ci_var._cml_linear_split
        elif isinstance(expr, mathml_piecewise):
            # The tests have to move into both components of gh:
            # "if C1 then (a1,b1) elif C2 then (a2,b2) else (a0,b0)"
            # maps to "(if C1 then a1 elif C2 then a2 else a0,
            #           if C1 then b1 elif C2 then b2 else b0)"
            # Note that no test is a function of var.
            # First rearrange child expressions
            pieces = getattr(expr, u'piece', [])
            cases = map(lambda p: child_i(p, 1), pieces)
            cases_ghs = map(lambda c: self._rearrange_expr(c, var), cases)
            if hasattr(expr, u'otherwise'):
                ow_gh = self._rearrange_expr(child_i(expr.otherwise, 1), var)
            else:
                ow_gh = (None, None)
            # Now construct the new expression
            conds = map(lambda p: self._clone(child_i(p, 2)), pieces)
            def piecewise_branch(i):
                pieces_i = zip(map(lambda gh: gh[i], cases_ghs),
                               conds)
                pieces_i = filter(lambda p: p[0] is not None,
                                  pieces_i) # Remove cases that are None
                ow = ow_gh[i]
                if pieces_i:
                    new_expr = mathml_piecewise.create_new(
                        expr, pieces_i, ow)
                elif ow:
                    new_expr = ow
                else:
                    new_expr = None
                return new_expr
            gh = (piecewise_branch(0), piecewise_branch(1))
            self._transfer_lut(expr, gh, var)
        elif isinstance(expr, mathml_apply):
            # Behaviour depends on the operator
            operator = expr.operator().localName
            operands = expr.operands()
            self.__expr = expr # For self._make_apply
            if operator in ['plus', 'minus']:
                # Just split the operation into each component
                operand_ghs = map(lambda op: self._rearrange_expr(op, var),
                                  operands)
                g = self._make_apply(operator, operand_ghs, 0)
                h = self._make_apply(operator, operand_ghs, 1)
                gh = (g, h)
            elif operator == 'divide':
                # (a, b) / (c, 0) = (a/c, b/c)
                numer = self._rearrange_expr(operands.next(), var)
                denom = self._rearrange_expr(operands.next(), var)
                g = h = None
                if numer[0]:
                    g = mathml_apply.create_new(expr, operator,
                                                [numer[0], denom[0]])
                if numer[1]:
                    h = mathml_apply.create_new(expr, operator,
                                                [numer[1], denom[0]])
                gh = (g, h)
            elif operator == 'times':
                # (a1,b1)*(a2,b2) = (a1*a2, b1*a2 or a1*b2 or None)
                # Similarly for the nary case - at most one b_i is not None
                operand_ghs = map(lambda op: self._rearrange_expr(op, var),
                                  operands)
                g = self._make_apply(operator, operand_ghs, 0,
                                     filter_none=False)
                # Find non-None b_i, if any
                for i, ab in enumerate(operand_ghs):
                    if ab[1] is not None:
                        operand_ghs[i] = (ab[1], None)
                        h = self._make_apply(operator, operand_ghs, 0)
                        break
                else:
                    h = None
                gh = (g, h)
            else:
                # (a, None) op (b, None) = (a op b, None)
                operand_ghs = map(lambda op: self._rearrange_expr(op, var),
                                  operands)
                g = self._make_apply(operator, operand_ghs, 0, preserve=True)
                gh = (g, None)
            self._transfer_lut(expr, gh, var)
        else:
            # Since this expression is linear, there can't be any
            # occurrence of var in it; all possible such cases are covered
            # above.  So just clone it.
            gh = (self._clone(expr), None)
        return gh
    
    def rearrange_linear_odes(self, doc):
        """Rearrange the linear ODEs so they can be updated directly
        on solving.

        Each ODE du/dt = f(u, t) can be written in the form
            du/dt = g(t) + h(t)u.
        A backward Euler update step is then as simple as
            u_n = (u_{n-1} + g(t)dt) / (1 - h(t)dt)
        (assuming that the transmembrane potential has already been
        solved for at t_n.

        Stores the results in doc.model._cml_linear_update_exprs, a
        mapping from variable object u to pair (g, h).
        """

        odes = map(lambda v: v._get_ode_dependency(doc.model._cml_free_var),
                   doc.model._cml_linear_vars)
        result = {}
        for var, ode in itertools.izip(doc.model._cml_linear_vars, odes):
            # Do the re-arrangement for this variable.
            # We do this by a post-order traversal of its defining ODE,
            # constructing a pair (g, h) for each subexpression recursively.
            # First, get the RHS of the ODE
            rhs = self._get_rhs(ode)
            # And traverse
            result[var] = self._rearrange_expr(rhs, var)
        # Store result in model
        doc.model._cml_linear_update_exprs = result
        return result

    def show(self, d):
        """Print out a more readable report on a rearrangement,
        as given by self.rearrange_linear_odes."""
        for var, expr in d.iteritems():
            print var.fullname()
            print "G:", expr[0].xml()
            print "H:", expr[1].xml()
            print "ODE:", var._get_ode_dependency(
                doc.model._cml_free_var).xml()
            print
