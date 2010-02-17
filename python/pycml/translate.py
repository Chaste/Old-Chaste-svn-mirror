#!/usr/bin/env python

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

# Translate CellML 1.0 to computer code
# Author: Jonathan Cooper

# This first version is rather hackish - demo of principle
# TODO:
#  - Allow easy output of other languages than C++ for Chaste
#  - Make file writes cope with non-ascii unicode?

# We want 1/2==0.5
from __future__ import division

# Common CellML processing stuff
import pycml
from pycml import *  # Put contents in the local namespace as well
import validator

import optparse
import re
import time


__version__ = "$Revision$"[11:-2]

def version_comment(note_time=True):
    """Generate a version comment, with optional time info."""
    if note_time:
        t = '\non ' + time.asctime()
    else:
        t = ''
    text = """Processed by pycml - CellML Tools in Python
    (translate: %s, pycml: %s)%s""" % (__version__, pycml.__version__, t)
    return text


def debugexpr(e):
    v = None
    if isinstance(e, cellml_variable):
        v = e
    elif isinstance(e, mathml_apply):
        v = e.assigned_variable()
    if v:
        r = (v==e, v.name, v.get_usage_count())
    else:
        r = (False, '', -1)
    return r


class TranslationError(RuntimeError):
    pass

class CellMLTranslator(object):
    """
    Base class for translators from CellML to programming languages.

    Provides various methods & attributes that can be overridden to
    achieve the desired output language and style.

    Also contains a registration system for subclasses, so the
    command-line client can know what translators are available.  See
    the register method for more information.
    """
    
    translators = {}
    class NameAlreadyRegistered(ValueError):
        pass
    @classmethod
    def register(cls, subclass, name):
        """Register a new translator subclass.

        Registers the subclass `subclass' with name `name' in the
        translators class attribute of CellMLTranslator.  If the name
        given is already in use, raises NameAlreadyRegistered.
        """
        if name in cls.translators:
            raise cls.NameAlreadyRegistered(name)
        cls.translators[name] = subclass
        return

    # Methods we could add:
    # start_func(name, args, rettype)
    # end_func()

    ###########################
    # Various language tokens #
    ###########################
    STMT_END = ';'        # End of statement
    EQ_ASSIGN = ' = '     # Assignment operator
    COMMENT_START = '// ' # Start of a 1 line comment
    DOXYGEN_COMMENT_START = '//! ' # Start of a 1 line Doxygen comment

    # Variable types
    TYPE_DOUBLE = 'double '
    TYPE_CONST_DOUBLE = 'const double '

    # Special constants
    TRUE = 'true'
    FALSE = 'false'
    PI = 'M_PI'
    E = 'M_E'
    
    # Whether the target language uses a subsidiary file, such as
    # a header file in C/C++
    USES_SUBSIDIARY_FILE = False
    # Mapping from primary file extension to subsidiary file extension
    FILE_EXTENSIONS = {'cpp': 'hpp',
                       'c': 'h',
                       'cxx': 'hxx'}

    def __init__(self, add_timestamp=True):
        """Create a translator."""
        # Initially output should not be indented
        self.indent_level = 0
        # Character to indent with
        self.indent_char = ' '
        # No. of occurrences of indent_char per indent_level
        self.indent_factor = 4
        # Whether to use lookup tables where possible
        self.use_lookup_tables = True
        # Whether to add a timestamp comment to generated files
        self.add_timestamp = add_timestamp
        # Main output goes to the main file by default
        self._main_output_to_subsidiary = False

    def error(self, lines, xml=None):
        """Raise a translation error.

        lines is a list of strings describing what went wrong.
        A TranslationError with that message will be raised.

        If xml is given, it should be an element, which will be
        pretty-printed and included in the error.
        """
        if xml is not None:
            lines.extend(xml.xml(indent = u'yes',
                                 omitXmlDeclaration = u'yes').split('\n'))
        raise TranslationError('\n'.join(lines))

    def translate(self, doc, model_filename, output_filename=None,
                  subsidiary_file_name=None,
                  class_name=None, v_variable=None,
                  continuation=None,
                  lookup_method_prefix='', row_lookup_method=False,
                  lt_index_uses_floor=True, bad_tables_for_cache=False,
                  constrain_table_indices=False):
        """Generate code for the given model.

        doc should be an instance of cellml_model representing a
        valid CellML model, such as might be produced from a call
        to
        >>> valid, doc = validator.CellMLValidator().validate(
        ...     model_filename, True)

        model_filename is the filename of the input model.
        The output program will by default be generated in the same
        folder, but with a different extension.  This can be
        overridden by supplying the output_filename keyword argument.

        By default the name of the class representing the model will
        be derived from the model name.  This can be overridden by
        passing an alternative as the class_name argument.

        The variable representing the transmembrane potential should
        be passed in using the v_variable argument.

        By default this method will perform some setup and then call
            self.output_top_boilerplate()
            self.output_mathematics()
            self.output_bottom_boilerplate()
        To alter this, pass a callable as the continuation parameter;
        this will then be called instead.

        lookup_method_prefix and row_lookup_method can be used to
        customise some aspects of lookup table usage.  The former is
        used by the Chaste translator to place lookup tables within a
        singleton class, while the latter can improve cache
        performance by looking up all tables in a single call, and
        returning an array of the results.
        
        lt_index_uses_floor specifies whether to use the floor()
        function to calculate the index into the lookup tables, or
        just cast to unsigned.

        constrain_table_indices specifies whether to throw an
        exception if lookup table index variables go outside the
        bounds specified (default), or just to cap them at the bounds.
        """
        self.doc = doc
        self.model = doc.model
        # Name of the class that will represent this model
        if class_name is None:
            self.class_name = u'CML_' + doc.model.name.replace('-', '_')
        else:
            self.class_name = class_name
        # Figure out the free & state vars in this model
        self.free_vars = doc.model.find_free_vars()
        self.state_vars = doc.model.find_state_vars()
        if len(self.free_vars) > 1:
            self.error(["Model has more than one free variable; exiting.",
                        "Free vars:" + str(free_vars)])
        # If only a single component, don't add it to variable names
        self.single_component = (len(getattr(self.model, u'component', []))
                                 == 1)
        # Find the (index of the) transmembrane potential
        self.v_variable = v_variable
        self.v_variable_name = (v_variable.component.name, v_variable.name)
        for i, var in enumerate(self.state_vars):
            if var is v_variable:
                self.v_index = i
                break
        else:
            self.v_index = -1
        # Check to see if we're using lookup tables
        self.lookup_method_prefix = lookup_method_prefix
        self.row_lookup_method = row_lookup_method
        self.lt_index_uses_floor = lt_index_uses_floor
        self.constrain_table_indices = constrain_table_indices
        self.bad_tables_for_cache = bad_tables_for_cache
        self.scan_for_lookup_tables()
        if not doc.lookup_tables:
            # No tables found
            self.use_lookup_tables = False

        # Extra configuration hook
        self.final_configuration_hook()

        # Open the output file(s)
        if output_filename is None:
            output_filename = self.output_file_name(model_filename)
        if self.USES_SUBSIDIARY_FILE:
            output_filename, self.subsidiary_filename = self.subsidiary_file_name(output_filename)
        self.out = open_output_stream(output_filename)
        if self.USES_SUBSIDIARY_FILE:
            self.out2 = open_output_stream(self.subsidiary_filename)

        # Translate
        if continuation:
            continuation()
        else:
            self.output_top_boilerplate()
            self.output_mathematics()
            self.output_bottom_boilerplate()
        close_output_stream(self.out)
        if self.USES_SUBSIDIARY_FILE:
            close_output_stream(self.out2)
        return

    def final_configuration_hook(self):
        """A hook for subclasses to do some final configuration."""
        return

    def output_file_name(self, model_filename):
        """Generate a name for our output file, based on the input file."""
        return os.path.splitext(model_filename)[0] + '.cpp'
    
    def subsidiary_file_name(self, output_filename):
        """Generate a name for the subsidiary output file, based on the main one.
        
        Returns a pair (main_output_file_name, subsidiary_file_name).  This is in
        case the user specifies (e.g.) a .hpp file as the main output - we consider
        the main output to be the .cpp file.
        """
        base, ext = os.path.splitext(output_filename)
        ext = ext[1:] # Remove the '.'
        try:
            new_ext = self.FILE_EXTENSIONS[ext]
            swap = False
        except KeyError:
            swap = True
            for key, val in self.FILE_EXTENSIONS:
                if val == ext:
                    new_ext = key
                    break
            else:
                # Turn off usage of subsidiary file
                self.USES_SUBSIDIARY_FILE = False
                return output_filename, None
        subsidiary_filename = base + '.' + new_ext
        if swap:
            output_filename, subsidiary_filename = subsidiary_filename, output_filename
        return output_filename, subsidiary_filename

    def send_main_output_to_subsidiary(self, to_subsidiary=True):
        """Set subsequent main-file writes to go to the subsidiary file instead.
        
        Supplying a False argument reverts this setting.
        
        Has no effect if not using a subsidiary file.
        """
        self._main_output_to_subsidiary = to_subsidiary

    def writeln(self, *args, **kwargs):
        """Write a line to our output file.

        Takes any number of strings as input, and writes them out to file.

        Unless the keyword argument indent is given as False, then the
        output will be indented to the level set by self.set_indent().
        Setting indent_level will override this value.
        Setting indent_offset will adjust the current level temporarily.

        If nl is set to False then a newline character will not be
        appended to the output.
        
        If subsidiary=True, then the line will be written to the subsidiary
        output file instead of the main one.  An error will be raised if
        there is no subsidiary output file.
        """
        if kwargs.get('subsidiary', False) or self._main_output_to_subsidiary:
            if not self.USES_SUBSIDIARY_FILE:
                self.error('Tried to write to non-existent subsidiary file')
            else:
                target = self.out2
        else:
            target = self.out
        indent = kwargs.get('indent', True)
        nl = kwargs.get('nl', True)
        if indent:
            if kwargs.has_key('indent_level'):
                level = kwargs['indent_level']
            else:
                level = self.indent_level
            if kwargs.has_key('indent_offset'):
                level += kwargs['indent_offset']
            target.write(self.indent_char * self.indent_factor * level)
        target.write(''.join(map(str, args)))
        if nl:
            target.write('\n')

    def write(self, *args):
        """Write to our output file.

        This variant does not indent the output, or add a newline.
        """
        self.writeln(indent=False, nl=False, *args)

    def output_comment(self, *args, **kwargs):
        """Output a (multi-line) string as a comment."""
        comment = ''.join(map(str, args))
        lines = comment.split('\n')
        for line in lines:
            self.writeln(self.COMMENT_START, line, **kwargs)

    def output_doxygen(self, *args, **kwargs):
        """Output a (multi-line) string as a Doxygen comment."""
        comment = ''.join(map(str, args))
        lines = comment.split('\n')
        for line in lines:
            self.writeln(self.DOXYGEN_COMMENT_START, line, **kwargs)

    def set_indent(self, level=None, offset=None):
        """Set the indentation level for subsequent writes.

        If offset is given, adjust the level by that amount, otherwise
        set it to an absolute value.
        """
        if offset is not None:
            self.indent_level += offset
        else:
            self.indent_level = level

    def code_name(self, var, ode=False, prefix=None):
        """
        Return the full name of var in a form suitable for inclusion in a
        source file.
        
        If ode is True then return the name of the derivative of var
        instead.  We go directly to the source variable in this case,
        rather than including intermediate assignment statements as is
        done for connections.
        """
        if prefix is None:
            prefix = ['var_', 'd_dt_'][ode]
        if ode:
            var = var.get_source_variable(recurse=True)
        if self.single_component:
            name = prefix + var.name
        else:
            name = prefix + var.xml_parent.name + '__' + var.name
        return name

    @property
    def include_guard(self):
        """
        Get the include guard (for C/C++ output) for this cell model,
        based on the class name.
        """
        return self.class_name.upper() + '_HPP_'
    
    def output_top_boilerplate(self):
        """Output top boilerplate."""
        self.writeln('#ifndef _', self.include_guard, '_')
        self.writeln('#define _', self.include_guard, '_\n')
        self.output_comment('Model: ', self.model.name)
        self.output_comment(version_comment(self.add_timestamp))
        self.writeln()
        self.writeln('#include <cmath>')
        self.writeln('#include "AbstractOdeSystem.hpp"')
        self.writeln('#include "Exception.hpp"')
        self.writeln('#include "AbstractStimulusFunction.hpp"\n')
        self.writeln('class ', self.class_name, ' : public AbstractOdeSystem')
        self.writeln('{')
        self.writeln('private:')
        self.writeln('AbstractStimulusFunction *mpStimulus;\n',
                     indent_offset=1)
        self.writeln('public:')
        self.set_indent(1)
        self.writeln('const static unsigned _NB_OF_STATE_VARIABLES_ = ',
                   str(len(self.state_vars)), ';\n')
        self.writeln('//', ('-'*66))
        self.writeln('// Methods')
        self.writeln('//', ('-'*66), '\n')
        # Constructor
        self.writeln('', self.class_name,
                   '(AbstractStimulusFunction *stim)')
        self.writeln('    : AbstractOdeSystem(_NB_OF_STATE_VARIABLES_, ',
                     self.v_index, ')')
        self.open_block()
        self.writeln('mpStimulus = stim;\n')
        for var in self.state_vars:
            self.writeln('mVariableNames.push_back("', var.name, '");')
            self.writeln('mVariableUnits.push_back("', var.units, '");')
            init_val = getattr(var, u'initial_value', None)
            if init_val is None:
                init_comm = ' // Value not given in model'
                # Don't want compiler error, but shouldn't be a real number
                init_val = 'NAN'
            else:
                init_comm = ''
            self.writeln('mInitialConditions.push_back(', init_val, ');',
                       init_comm, '\n')
        if self.use_lookup_tables: self.output_lut_generation()
        self.close_block()
        # Destructor
        self.writeln('~', self.class_name, '(void)')
        self.open_block()
        if self.use_lookup_tables: self.output_lut_deletion()
        self.close_block()
        # Lookup table declarations & methods
        if self.use_lookup_tables:
            self.output_lut_declarations()
            self.output_lut_methods()
        # Evaluation function
        self.writeln('void EvaluateYDerivatives (')
        self.writeln('        double ', self.code_name(self.free_vars[0]), ',')
        self.writeln('        const std::vector<double> &rY,')
        self.writeln('        std::vector<double> &rDY)')
        self.open_block()
        self.writeln('// Inputs:')
        self.writeln('// Time units: ', self.free_vars[0].units)
        for i, var in enumerate(self.state_vars):
            self.writeln('double ', self.code_name(var),
                         ' = rY[', str(i), '];')
            self.writeln('// Units: ', var.units, '; Initial value: ',
                         getattr(var, u'initial_value', 'Unknown'))
        self.writeln()
        if self.use_lookup_tables:
            self.output_table_index_generation()
        return

    def output_mathematics(self):
        """Output the mathematics in this model."""
        self.writeln(self.COMMENT_START, 'Mathematics')
        for expr in self.model.get_assignments():
            # Check this expression is actually used; don't output if not
            var = None
            if isinstance(expr, mathml_apply) and expr.is_assignment():
                var = expr.assigned_variable()
            elif isinstance(expr, cellml_variable):
                var = expr
            if not (var and var.get_usage_count() == 0):
                self.output_assignment(expr)
        return

    def output_bottom_boilerplate(self):
        """Output bottom boilerplate"""
        self.writeln('\n')
        for i, var in enumerate(self.state_vars):
            self.writeln('rDY[', str(i), '] = ', self.code_name(var, True),
                         ';')
        self.close_block()
        self.set_indent(offset=-1)
        self.writeln('};\n')
        self.writeln('#endif')
        return

    def output_assignment(self, expr):
        """Output an assignment expression."""
        if isinstance(expr, cellml_variable):
            # This may be the assignment of a mapped variable, or a constant
            t = expr.get_type()
            if t == VarTypes.Mapped:
                # TODO: Use const double instead of double?
                self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(expr),
                             self.EQ_ASSIGN,
                             self.code_name(expr.get_source_variable()),
                             self.STMT_END)
            elif t == VarTypes.Constant:
                self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(expr),
                             self.EQ_ASSIGN, nl=False)
                self.output_number(expr.initial_value)
                self.writeln(self.STMT_END, indent=False)
        else:
            # This is a mathematical expression
            # TODO: Use const double instead of double?
            self.writeln(self.TYPE_CONST_DOUBLE, nl=False)
            opers = expr.operands()
            self.output_lhs(opers.next())
            self.write(self.EQ_ASSIGN)
            self.output_expr(opers.next(), False)
            self.writeln(self.STMT_END, indent=False)

    def output_lhs(self, expr):
        """Output the left hand side of an assignment expression."""
        if expr.localName == 'ci':
            self.output_variable(expr)
        elif expr.operator().localName == 'diff':
            self.write(self.code_name(expr.operator().dependent_variable,
                                      ode=True))

    def output_variable(self, ci_elt):
        """Output a ci element, i.e. a variable lookup."""
        self.write(self.code_name(ci_elt.variable))

    def output_expr(self, expr, paren):
        """Output the expression expr.
        
        If paren is True then the context has requested parentheses around the
        output; if expr requires them then they will be added.
        """
        if self.use_lookup_tables and self.is_lookup_table(expr):
            self.output_table_lookup(expr, paren)
        elif isinstance(expr, mathml_apply):
            self.output_apply(expr, paren)
        elif isinstance(expr, mathml_piecewise):
            self.output_piecewise(expr, paren)
        elif isinstance(expr, mathml_ci):
            self.output_variable(expr)
        elif expr.localName == u'cn':
            self.output_number(expr)
        elif expr.localName == u'degree':
            # <degree> is just a wrapper around an expression
            self.output_expr(child_i(expr, 1), paren)
        elif expr.localName == u'logbase':
            # <logbase> is just a wrapper around an expression
            self.output_expr(child_i(expr, 1), paren)
        elif expr.localName == u'true':
            self.write(self.TRUE)
        elif expr.localName == u'false':
            self.write(self.FALSE)
        elif expr.localName == u'pi':
            self.write(self.PI)
        elif expr.localName == u'exponentiale':
            self.write(self.E)
        else:
            self.error(["Unsupported expression element " + expr.localName],
                       xml=expr)

    def output_number(self, expr):
        """Output the plain number expr.
        
        We make all constants parse as doubles to avoid problems with
        integer division or numbers too large for the int type.
        
        Negative numbers will be prefixed by a space to avoid unwanted
        decrement operations.
        """
        n = self.eval_number(expr)
        num = "%.12g" % n
        if num[0] == '-':
            num = ' ' + num
        if not '.' in num and not 'e' in num:
            num = num + '.0'
        self.write(num)

    def eval_number(self, expr):
        """Evaluate a number.

        If a (unicode) string, convert to float.
        If a cn element, call its evaluate method.
        """
        if isinstance(expr, mathml_cn):
            return expr.evaluate()
        else:
            return float(unicode(expr))

    # Map from operator element names to C++ function names,
    # where the translation is straightforward.
    function_map = {'power': 'pow', 'abs': 'fabs', 'ln': 'log', 'exp': 'exp',
                    'floor': 'floor', 'ceiling': 'ceil',
                    'factorial': 'factorial', # Needs external definition
                    'not': '!',
                    'sin': 'sin', 'cos': 'cos', 'tan': 'tan',
                    'sec': '1/cos', 'csc': '1/sin', 'cot': '1/tan',
                    'sinh': 'sinh', 'cosh': 'cosh', 'tanh': 'tanh',
                    'sech': '1/cosh', 'csch': '1/sinh', 'coth': '1/tanh',
                    'arcsin': 'asin', 'arccos': 'acos', 'arctan': 'atan',
                    'arcsinh': 'asinh', 'arccosh': 'acosh', 'arctanh': 'atanh'}
    # Inverse reciprocal trig functions; these are represented by
    # key(x) = function_map[val](1/x)
    recip_trig = {'arcsec': 'arcsin', 'arccsc': 'arccos', 'arccot': 'arctan',
                  'arcsech': 'arcsinh', 'arccsch': 'arccosh', 'arccoth': 'arctanh'}
    # Operators
    nary_ops   = {'plus': '+', 'times': '*',
                  'and': '&&', 'or': '||'}
    binary_ops = {'divide': '/',
                  'xor': '!=', 'eq': '==', 'neq': '!=',
                  'geq': '>=','leq': '<=','gt': '>','lt': '<'}

    def output_apply(self, expr, paren):
        """Output an <apply> expression.
        
        paren is True if the context has requested parentheses.
        """
        op = expr.operator()
        if self.function_map.has_key(op.localName):
            self.output_function(self.function_map[op.localName],
                                 expr.operands(), paren)
        elif self.recip_trig.has_key(op.localName):
            self.output_function(
                self.function_map[self.recip_trig[op.localName]],
                expr.operands(), paren, reciprocal=True)
        elif op.localName == u'root':
            self.output_root(expr, paren)
        elif op.localName == u'log':
            self.output_log(expr, paren)
        elif self.nary_ops.has_key(op.localName):
            self.output_nary_operator(self.nary_ops[op.localName],
                                      expr.operands(), paren)
        elif self.binary_ops.has_key(op.localName):
            self.output_binary_operator(self.binary_ops[op.localName],
                                        expr.operands(), paren, expr)
        elif op.localName == u'minus':
            self.output_minus(expr, paren)
        elif op.localName == u'diff':
            # ODE occuring on the RHS
            self.write(self.code_name(op.dependent_variable, ode=True))
        else:
            # Unrecognised operator
            self.error(["Unsupported operator element " + str(op.localName)],
                       xml=expr)

    def output_function(self, func_name, args, paren, reciprocal=False):
        """Output a function call with name func_name and arguments args.
        
        Parentheses are not required so paren is ignored.
        If reciprocal is True then pass the reciprocal of each arg to
        func_name.
        """
        self.write(func_name + '(')
        comma = False
        for arg in args:
            if comma: self.write(', ')
            else: comma = True
            if reciprocal:
                self.write('1/')
                self.output_expr(arg, True)
            else:
                self.output_expr(arg, False)
        self.write(')')

    def output_nary_operator(self, operator, operands, paren):
        """Output an n-ary operator (using infix notation).
        
        If paren is True, enclose the output in parentheses.
        """
        # TODO: Optimise - to use expm1(x) for computing exp(x)-1
        self.open_paren(paren)
        op = False
        for operand in operands:
            if op: self.write(' ' + operator + ' ')
            else: op = True
            self.output_expr(operand, True)
        self.close_paren(paren)

    def output_unary_operator(self, operator, operand, paren):
        """Output a unary operator (using prefix notation)."""
        self.open_paren(paren)
        self.write(operator)
        self.output_expr(operand, True)
        self.close_paren(paren)

    def output_binary_operator(self, operator, operands, paren, expr):
        """Output a binary operator.
        
        As output_nary_operator, but checks that len(list(operands)) == 2.
        """
        operands = list(operands)
        if len(operands) != 2:
            self.error(["Binary operator" + operator +
                        "does not have 2 operands."], xml=expr)
        self.output_nary_operator(operator, operands, paren)

    special_roots = {2: 'sqrt', 3: 'cbrt'}
    def output_root(self, expr, paren):
        """Output a root taken to some degree.

        If a degree qualifier element is not provided, uses default 2.
        """
        if hasattr(expr, u'degree'):
            # A degree is given.  Compute x^(1/b)
            # TODO: Optimise for when b==2 (sqrt) or b==3 (cbrt)
            # Try to evaluate expr.degree, and if the result is a key
            # of self.special_roots, use the value as the function to call.
            self.write('pow(')
            self.output_expr(expr.operands().next(), False)
            self.write(', 1/')
            self.output_expr(expr.degree, True)
            self.write(')')
        else:
            # Compute square root
            self.output_function('sqrt', expr.operands(), paren)

    def output_log(self, expr, paren):
        """Output a logarithm to the given base, which defaults to base 10."""
        if hasattr(expr, u'logbase'):
            # A base is provided.  Use the identity log_b(x) = log(x)/log(b)
            # TODO: Optimise for log2(x)
            self.open_paren(paren)
            self.output_function('log', expr.operands(), paren)
            self.write('/log(')
            self.output_expr(expr.logbase, False)
            self.write(')')
            self.close_paren(paren)
        else:
            # Use base 10
            self.output_function('log10', expr.operands(), paren)

    def output_minus(self, expr, paren):
        """Output either a unary or binary minus.

        Which is chosen depends on the number of operands.
        """
        operands = list(expr.operands())
        if len(operands) == 1:
            self.output_unary_operator('-', operands[0], paren)
        else:
            self.output_binary_operator('-', operands, paren, expr)

    def output_piecewise(self, expr, paren):
        """Output the piecewise expression expr.

        We use a cascading ternary if expression for simplicity.
        """
        self.open_paren(paren)
        for piece in getattr(expr, u'piece', []):
            self.output_expr(child_i(piece, 2), True) # Condition
            self.write(' ? ')
            self.output_expr(child_i(piece, 1), True) # Result
            self.write(' : ')
        if hasattr(expr, u'otherwise'):
            self.output_expr(child_i(expr.otherwise, 1), True) # Default case
        self.close_paren(paren)

    def open_paren(self, paren):
        if paren: self.write('(')
    def close_paren(self, paren):
        if paren: self.write(')')

    def open_block(self, **kwargs):
        """Open a new code block and increase indent."""
        self.writeln('{', **kwargs)
        self.set_indent(offset=1)
    def close_block(self, blank_line=True, **kwargs):
        """Close a code block and decrease indent."""
        self.set_indent(offset=-1)
        self.writeln('}', **kwargs)
        if blank_line:
            self.writeln(**kwargs)
        return

    ##############################
    # Dependency related methods #
    ##############################

    # These methods allow us to calculate which equations must be
    # output in order to compute a given set of quantities.
    def calculate_extended_dependencies(self, nodes, prune=[]):
        """Calculate the extended dependencies of the given nodes.

        Recurse into the dependency graph, in order to construct a
        set, for each node in nodes, of all the nodes on which it
        depends, either directly or indirectly.

        Each node IS included in its own dependency set.

        If prune is specified, it should be a set of nodes for which
        we won't include their dependencies or the nodes themselves.
        This is useful for pruning variables required for calculating
        a stimulus if the stimulus is being provided by another
        method.

        Requires the dependency graph to be acyclic.

        Return the union of all the dependency sets.
        """
        deps = set()
        for node in nodes:
            if node in prune or (isinstance(node, mathml_apply) and
                                 isinstance(node.operator(), mathml_eq) and
                                 isinstance(node.eq.lhs, mathml_ci) and
                                 node.eq.lhs.variable in prune):
                continue
            if type(node) == tuple:
                # This is an ODE dependency, so get the defining expression
                # instead.
                ode = True
                node = node[0]._get_ode_dependency(node[1])
                free_var = node.eq.lhs.diff.independent_variable
            else:
                ode = False
            deps.add(node)
            nodedeps = node._get_dependencies()[:]
            if ode and not node._cml_ode_has_free_var_on_rhs:
                # ODEs depend on their independent variable.  However,
                # when writing out code we don't want to pull the free
                # variable in just for this, as the compiler then
                # gives us unused variable warnings.
                nodedeps.remove(free_var)
            deps.update(self.calculate_extended_dependencies(nodedeps,
                                                             prune=prune))
        return deps

    def output_equations(self, nodeset):
        """Output the mathematics described by nodeset.

        nodeset represents a subset of the assignments in the model.
        Output assignments in the order given by a topological sort,
        but only include those in nodeset.

        Since a set of assignments is given, this method does not
        check whether variables are used - it is assumed that only
        assignments that are actually wanted are given in nodeset.
        """
        for expr in (e for e in self.model.get_assignments() if e in nodeset):
            self.output_assignment(expr)
        return
    
    def get_var_units(self, varobj):
        """Get the units for the given cellml_variable."""
        assert isinstance(varobj, cellml_variable)
        uname = varobj.units
        units = varobj.component.get_units_by_name(uname)
        return units

    def varobj(self, varname):
        """Return the variable object that has code_name varname."""
        if varname[0] == '(':
            cname, vname = varname[1:-1].split(',')
            name = vname
        elif '__' in varname:
            if varname[:4] == 'var_':
                cname, vname = varname[4:].split('__')
            else:
                cname, vname = varname.split('__')
            name = cname + u'__' + vname
        else:
            return None
        if self.single_component:
            var = self.model.component.get_variable_by_name(name)
        else:
            var = self.model.get_variable_by_name(cname, vname)
        return var

    def _vars_in(self, expr):
        """Return a list of variable objects used in the given expression.

        DO include state variables."""
        res = set()
        if isinstance(expr, mathml_ci):
            varname = unicode(expr)
            varobj = self.varobj(varname.strip())
            #if varobj not in self.state_vars:
            if varobj:
                res.add(varobj)
        elif isinstance(expr, mathml_apply) and \
                 expr.operator().localName == u'diff':
            dep_varname = unicode(expr.ci)
            varobj = self.varobj(dep_varname.strip())
            res.add(varobj._get_ode_dependency(self.free_vars[0]))
        elif hasattr(expr, 'xml_children'):
            for child in expr.xml_children:
                res.update(self._vars_in(child))
        return res


    ########################
    # Lookup table methods #
    ########################

    # Lookup tables should be done in a cache- and memory-
    # efficient manner.
    #
    # Cache: Have one block of memory allocated for all tables with a
    # given index variable, such that entries are found at i*N+j,
    # where N is the no. of tables in the block, i is the index into a
    # table, and j is the table to read.  Change how lookups are done,
    # such that the lookup method is called once and returns a pointer
    # to the (i*N)'th entry.  Places where we now call the method then
    # index this pointer by j.
    #    The 'one block' part is done by default.
    #    The better lookup method is activated by --row-lookup-method.
    #
    # Memory: Extract the lookup tables into a separate class (in the
    # same .hpp file though).  This can then be made a singleton class
    # in a multi-cellular context.
    #    Chaste code generation has the option to do this, enabled by
    #    default.  Use --no-separate-lut-class to disable.

    def scan_for_lookup_tables(self):
        """Search for lookup tables used in this document.

        Store a list of suitable expressions on the document root.
        Generate a dictionary mapping tables to their index variables.
        """
        doc = self.doc
        # Get list of suitable expressions
        doc.lookup_tables = doc.xml_xpath(u"//*[@lut:possible='yes']")
        doc.lookup_tables.sort(cmp=element_path_cmp)
        # Map table keys (min, max, step, var) to an index variable
        doc.lookup_table_indexes = {}
        # Count the no. of lookup tables with each index variable
        doc.lookup_tables_num_per_index = {}
        if not doc.lookup_tables:
            # If no suitable expressions, we're done
            return
        # Search for table index variables already assigned
        table_indexes = [int(getattr(expr, u'table_index', -1))
                         for expr in doc.lookup_tables]
        tidx = max(table_indexes) + 1
        # Search for table names already assigned
        table_numbers = {}
        for expr in doc.lookup_tables:
            if hasattr(expr, u'table_name'):
                idx = expr.table_index
                table_numbers[idx] = max(table_numbers.get(idx, 0),
                                         1 + int(expr.table_name))
        # Now assign new names, and construct mapping from tables to
        # index variables
        for expr in doc.lookup_tables:
            # Get a suitable table index variable
            comp = expr.get_component()
            var = comp.get_variable_by_name(expr.var)
            var = var.get_source_variable(recurse=True)
            key = (expr.min, expr.max, expr.step, var)
            if not doc.lookup_table_indexes.has_key(key):
                if hasattr(expr, u'table_index'):
                    doc.lookup_table_indexes[key] = expr.table_index
                else:
                    doc.lookup_table_indexes[key] = unicode(tidx)
                    tidx += 1
                    expr.xml_set_attribute((u'lut:table_index', NSS['lut']),
                                           doc.lookup_table_indexes[key])
                doc.lookup_tables_num_per_index[doc.lookup_table_indexes[key]] = 0
            # Get a table name, unique over all tables with
            # this index variable
            if not hasattr(expr, u'table_name'):
                tnum = table_numbers.get(doc.lookup_table_indexes[key], 0)
                expr.xml_set_attribute((u'lut:table_name', NSS['lut']),
                                       unicode(tnum))
                table_numbers[doc.lookup_table_indexes[key]] = tnum + 1
        # Make sure each lookup table is only listed once in doc.lookup_tables,
        # so we don't get 2 tables for the same expression!
        candidates = doc.lookup_tables[:]
        doc.lookup_tables = []
        listed = set()
        for expr in candidates:
            tid = (expr.table_index, expr.table_name)
            if not tid in listed:
                listed.add(tid)
                doc.lookup_tables.append(expr)
                # Increment count for this index variable
                doc.lookup_tables_num_per_index[expr.table_index] += 1
        return

    def lut_access_code(self, table_index, table_name, i):
        """Get the code for accessing the i'th element of the given table.
        """
        name = '_lookup_table_' + str(table_index)
        if self.bad_tables_for_cache:
            name += '[' + str(table_name) + '][' + str(i) + ']'
        else:
            name += '[' + str(i) + '][' + str(table_name) + ']'
        return name

    def output_lut_generation(self):
        """Output code to generate lookup tables.

        There should be a list of suitable expressions available as
        self.doc.lookup_tables, to save having to search the whole
        model.
        """
        # Don't use table lookups to generate the tables!
        self.use_lookup_tables = False
        # Generate each table in a separate loop
        for expr in self.doc.lookup_tables:
            table_size = unicode(1+int((float(expr.max) - float(expr.min)) /
                                       float(expr.step)))
            j = expr.table_name
            var = expr.get_component().get_variable_by_name(expr.var)
            varname = self.code_name(var)
            idx = self.doc.lookup_table_indexes[
                (expr.min, expr.max, expr.step,
                 var.get_source_variable(recurse=True))]
            self.writeln('for (int i=0 ; i<', table_size, '; i++)')
            self.open_block()
            self.writeln(self.TYPE_DOUBLE, varname, self.EQ_ASSIGN, expr.min,
                         ' + i*', expr.step, self.STMT_END)
            self.writeln(self.lut_access_code(idx, j, 'i'), self.EQ_ASSIGN,
                         nl=False)
            self.output_expr(expr, False)
            self.writeln(self.STMT_END, indent=False)
            self.close_block()
        self.use_lookup_tables = True
        return

    def output_lut_deletion(self):
        """Output code to delete memory allocated for lookup tables.

        Does nothing, since memory is allocated statically.
        """
        return

    def output_lut_declarations(self, indexes_as_member=False):
        """Output declarations for the lookup tables."""
        self.output_comment('Lookup tables')
        # Allocate memory, per index variable for cache efficiency
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            min, max, step, var = key
            table_size = unicode(1+int((float(max) - float(min)) / float(step)))
            num_tables = unicode(self.doc.lookup_tables_num_per_index[idx])
            self.writeln(self.TYPE_DOUBLE,
                         self.lut_access_code(idx, num_tables, table_size),
                         self.STMT_END)
            if indexes_as_member:
                self.writeln('unsigned _table_index_', idx, ';')
                self.writeln('double _factor_', idx, ';')
        self.writeln()
        return

    def output_lut_indices(self):
        """Output declarations for the lookup table indices."""
        self.output_comment('Lookup table indices')
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            self.writeln('unsigned _table_index_', idx, self.STMT_END)
            self.writeln('double _factor_', idx, self.STMT_END)
        self.writeln()
        return
    
    def output_lut_methods(self):
        """Output the methods which look up values from lookup tables."""
        if self.row_lookup_method:
            self.output_lut_row_lookup_methods()
            return
        self.output_comment('Methods to look up values from lookup tables')
        self.output_comment('using linear interpolation')
        for expr in self.doc.lookup_tables:
            j = expr.table_name
            idx = expr.table_index
            self.writeln('inline double _lookup_', j,
                         '(unsigned i, double factor)')
            self.open_block()
            self.writeln(self.TYPE_DOUBLE, 'y1', self.EQ_ASSIGN,
                         self.lut_access_code(idx, j, 'i'), self.STMT_END)
            self.writeln(self.TYPE_DOUBLE, 'y2', self.EQ_ASSIGN,
                         self.lut_access_code(idx, j, 'i+1'), self.STMT_END)
            self.writeln('return y1 + (y2-y1)*factor;')
            self.close_block()
        self.writeln()
        return

    def output_lut_row_lookup_methods(self):
        """Write methods that return a whole row of a lookup table.

        Note: assumes that table names are numbered sequentially from 0.
        """
        self.output_comment('Row lookup methods')
        self.output_comment('using linear interpolation')
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            num_tables = unicode(self.doc.lookup_tables_num_per_index[idx])
            self.writeln('double* _lookup_', idx,
                         '_row(unsigned i, double factor)')
            self.open_block()
            self.writeln('for (unsigned j=0; j<', num_tables, '; j++)')
            self.open_block()
            self.writeln(self.TYPE_DOUBLE, 'y1', self.EQ_ASSIGN,
                         self.lut_access_code(idx, 'j', 'i'), self.STMT_END)
            self.writeln(self.TYPE_DOUBLE, 'y2', self.EQ_ASSIGN,
                         self.lut_access_code(idx, 'j', 'i+1'), self.STMT_END)
            self.writeln('_lookup_table_', idx, '_row[j]',
                         ' = y1 + (y2-y1)*factor;')
            self.close_block(blank_line=False)
            self.writeln('return _lookup_table_', idx, '_row;')
            self.close_block()
        self.writeln()
        return
    
    def output_lut_row_lookup_memory(self):
        """Output declarations for the memory used by the row lookup methods."""
        self.output_comment('Row lookup methods memory')
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            min, max, step, var = key
            table_size = unicode(1+int((float(max) - float(min)) / float(step)))
            num_tables = unicode(self.doc.lookup_tables_num_per_index[idx])
            self.writeln('double _lookup_table_', idx, '_row[', num_tables,
                         '];')
        self.writeln()
        return

    def is_lookup_table(self, expr):
        """Return True iff expr can be replaced by a lookup table.

        Uses annotations from a previous analysis."""
        return expr.getAttributeNS(NSS['lut'], u'possible', '') == u'yes'

    def output_table_lookup(self, expr, paren):
        """Output code to look up expr in the appropriate table."""
        i = expr.table_index
        if self.row_lookup_method:
            self.write('_lt_', i, '_row[', expr.table_name, ']')
        else:
            self.write(self.lookup_method_prefix, '_lookup_', expr.table_name,
                       '(_table_index_', i, ', _factor_', i, ')')
        return

    def output_table_index_generation(self, indexes_as_member=False):
        """Output code to calculate indexes into any lookup tables."""
        self.output_comment('Lookup table indexing')
        if indexes_as_member:
            index_type = ''
            factor_type = ''
        else:
            index_type = 'unsigned '
            factor_type = 'double '
        for key, i in self.doc.lookup_table_indexes.iteritems():
            min, max, step, var = key
            step_inverse = unicode(1 / float(step))
            offset = '_offset_' + i
            offset_over_step = offset + '_over_table_step'
            varname = self.code_name(var)
            self.writeln('if (', varname, '>', max,
                         ' || ', varname, '<', min, ')')
            self.open_block()
            self.writeln('#define COVERAGE_IGNORE', indent=False)
            if self.constrain_table_indices:
                self.writeln('if (', varname, '>', max, ') ', varname,
                             ' = ', max, ';')
                self.writeln('else ', varname, ' = ', max, ';')
            else:
                self.writeln('EXCEPTION(DumpState("V outside lookup table range", rY));')
            self.writeln('#undef COVERAGE_IGNORE', indent=False)
            self.close_block()
            self.writeln('double ', offset, ' = ', varname, ' - ', min, ';')
            self.writeln('double ', offset_over_step, ' = ', offset, ' * ',
                         step_inverse, ';')
            if self.lt_index_uses_floor:
                self.writeln(index_type, '_table_index_', i,
                             ' = (unsigned) floor(', offset_over_step, ');')
            else:
                self.writeln(index_type, '_table_index_', i,
                             ' = (unsigned)(', offset_over_step, ');')
            self.writeln(factor_type, '_factor_', i, ' = ', offset_over_step,
                         ' - _table_index_', i, ';')
            if self.row_lookup_method:
                self.writeln('double* _lt_', i, '_row = ',
                             self.lookup_method_prefix, '_lookup_', i,
                             '_row(_table_index_', i, ', _factor_', i, ');')
        self.writeln()
        return


class CellMLToChasteTranslator(CellMLTranslator):
    """
    As CellMLTranslator, but targets more recent Chaste style.

    Includes the ability to output a cell that can solve itself using
    backward Euler, if the appropriate analyses have been done on the
    model.  (See the -J and -j options to translate.py.)
    """

    # We want separate .cpp/.hpp files
    USES_SUBSIDIARY_FILE = True

    # Type of a reference to the state variable vector
    TYPE_VECTOR_REF = 'std::vector<double>&'
    
    def writeln_hpp(self, *args, **kwargs):
        """Convenience wrapper for writing to the header file."""
        kwargs['subsidiary'] = True
        self.writeln(*args, **kwargs)

    def translate(self, *args, **kwargs):
        """Generate code for the given model."""
        our_kwargs = {'use_chaste_stimulus': False,
                      'separate_lut_class': True,
                      'convert_interfaces': False,
                      'kept_vars_as_members': True,
                      'use_metadata': False,
                      'dynamically_loadable': False
                      }
        for key, default in our_kwargs.iteritems():
            setattr(self, key, kwargs.get(key, default))
            if key in kwargs:
                del kwargs[key]
        # Last method's access specification
        self._last_method_access = 'private'
        return super(CellMLToChasteTranslator, self).translate(*args, **kwargs)

    def final_configuration_hook(self):
        """Set the LT method prefix (requires self.class_name to be set)."""
        if self.separate_lut_class:
            self.lt_class_name = self.class_name + '_LookupTables'
            self.lookup_method_prefix = self.lt_class_name + '::Instance()->'
        return super(CellMLToChasteTranslator, self).final_configuration_hook()
    
    def ionic_current_units_conversion_factor(self, units):
        """
        Check whether the units of the transmembrane current are as
        expected by Chaste.  If they are dimensionally consistent,
        return the required conversion factor (could be 1).  If not,
        print a warning message to stderr.

        To go from model -> Chaste, multiply by the conversion factor.
        """
        chaste_units = cellml_units.create_new(
            self.model, 'uA_per_cm2',
            [{'units': 'ampere', 'prefix': 'micro'},
             {'units': 'metre', 'prefix': 'centi', 'exponent': '-2'}])
        if not chaste_units.dimensionally_equivalent(units):
            print >>sys.stderr, "Units of the ionic current are not in the " \
                "dimensions expected by Chaste (uA/cm^2).  Please add " \
                "a suitable conversion manually - we cannot do this " \
                "reliably automatically."
            return 1
        return (units.get_multiplicative_factor() /
                chaste_units.get_multiplicative_factor())
    
    def check_time_units(self):
        """Check if units conversions at the interfaces are required (#621).
        
        Sets self.conversion_factor to None if conversions are not required.
        If they are, to go from model -> Chaste, divide by the conversion factor.
        To go from Chaste -> model, multiply.
        If dealing with a derivative, do the opposite.
        
        TODO: Change the way we approach this ticket: have a method, called here,
        which adds conversion mathematics to the model itself, so we don't
        have to check in each place code is generated?  Would have to add
        new variables which Chaste reads/writes.  The function
        mathml_units_mixin._add_units_conversion could then be used.
        An alternative plan is to put the conversions in
        Compute(ExceptVoltage), but these may not be virtual.
        """
        if self.convert_interfaces:
            time_units = self.get_var_units(self.free_vars[0])
            # Create a reference 'millisecond' units definition
            ms = cellml_units.create_new(
                self.model, 'milliseconds',
                [{'units': 'second', 'prefix': 'milli'}])
            # Compare
            if not ms.dimensionally_equivalent(time_units):
                # Oops!
                raise TranslationError('Time does not have dimensions of time')
            elif ms != time_units:
                # We'll need conversions
                DEBUG('chaste-translator', 'Will convert time units')
                self.conversion_factor = (ms.get_multiplicative_factor()/
                                          time_units.get_multiplicative_factor())
            else:
                self.conversion_factor = None
            # Now to go from model -> Chaste, divide by the conversion factor.
            # To go from Chaste -> model, multiply.
            # If dealing with a derivative, do the opposite.
        else:
            self.conversion_factor = None
        return
    
    def output_includes(self, base_class=None):
        """Output the start of each output file.
        
        As well as the #include lines, it also outputs the include guard for
        the .hpp file, and doxygen comment.
        
        If base_class is not None (and self.use_backward_euler isn't set)
        then includes that class' header instead of AbstractCardiacCell.
        
        Reads self.include_serialization and self.use_backward_euler.
        Sets self.base_class_name.
        """
        self.writeln_hpp('#ifndef ', self.include_guard)
        self.writeln_hpp('#define ', self.include_guard, '\n')
        for sub in [False, True]:
            self.output_doxygen('@file\n\n',
                                'This source file was generated from CellML.\n\n',
                                'Model: ', self.model.name, '\n\n',
                                version_comment(self.add_timestamp),
                                '\n\n<autogenerated>',
                                subsidiary=sub)
            self.writeln(subsidiary=sub)
        # .cpp should include .hpp
        self.writeln('#include "', os.path.basename(self.subsidiary_filename), '"')
        if self.include_serialization:
            self.writeln_hpp('#include <boost/serialization/access.hpp>')
            self.writeln_hpp('#include <boost/serialization/base_object.hpp>')
        self.writeln('#include <cmath>')
        self.writeln('#include <cassert>')
        if self.use_backward_euler:
            self.writeln_hpp('#include "AbstractBackwardEulerCardiacCell.hpp"')
            self.writeln('#include "CardiacNewtonSolver.hpp"')
            self.base_class_name = 'AbstractBackwardEulerCardiacCell<' + \
                str(self.nonlinear_system_size) + '>'
        elif base_class:
            self.base_class_name = base_class
            self.writeln_hpp('#include "' + self.base_class_name + '.hpp"')
        else:
            self.base_class_name = 'AbstractCardiacCell'
            self.writeln_hpp('#include "' + self.base_class_name + '.hpp"')
        if self.use_metadata:
            self.writeln_hpp('#include "AbstractCardiacCellWithModifiers.hpp"')
            # Modify the base class name
            self.base_class_name = 'AbstractCardiacCellWithModifiers<' + self.base_class_name + ' >'
        self.writeln('#include "Exception.hpp"')
        self.writeln('#include "OdeSystemInformation.hpp"')
        self.writeln_hpp('#include "AbstractStimulusFunction.hpp"')
        self.writeln()
        self.writeln_hpp()
        
    def set_access(self, access):
        """Set the access specification for subsequent output.
        
        We keep track of the last access set, either via this method or
        output_method_start, and only output a new declaration to the
        header file if it changes.
        """
        if access != self._last_method_access:
            self._last_method_access = access
            self.writeln_hpp()
            self.writeln_hpp(access, ':', indent_offset=-1)

    def output_method_start(self, method_name, args, ret_type, access=None):
        """Output the start of a method declaration/definition.
        
        Will write to both the .hpp and .cpp file.
        
        We keep track of the access of the last method, and only output a new
        declaration to the header file if it changes.  The default is to use
        the same access specification as last time.
        """
        if access:
            self.set_access(access)
        if ret_type:
            if ret_type[-1] != ' ':
                ret_type = ret_type + ' '
        else:
            ret_type = ''
        args_string = ', '.join(args)
        self.writeln_hpp(ret_type, method_name, '(', args_string, ')', self.STMT_END)
        self.writeln(ret_type, self.class_name, '::', method_name, '(', args_string, ')')

    def output_cell_parameters(self):
        """Output declarations, set & get methods for cell parameters.
        
        "Parameters" are those variables annotated with pe:keep (see #666).
        They can be real parameters, which have both set & get methods,
        or computed variables, which just have get methods.
        
        Sets self.cell_parameters to be a list of the cellml_variable objects
        representing parameters.
        
        If we're using metadata (self.use_metadata = True), then also
        collect variables annotated with oxmeta:name into self.metadata_vars.
        Only constants and state variables are included.
        """
        # Find annotated variables
        if self.kept_vars_as_members:
            kept_vars = self.doc.xml_xpath(u'/*/*/cml:variable[@pe:keep]')
        else:
            kept_vars = []
        self.cell_parameters = kept_vars
        # create set of all annotated variables         
        if self.use_metadata:
            vars = self.doc.xml_xpath(u'/*/*/cml:variable[@oxmeta:name]')
            self.metadata_vars = set([v for v in vars if v.get_type() == VarTypes.Constant
                                      or v in self.state_vars])
        else:
            self.metadata_vars = set([])

        # Generate member variable declarations
        self.set_access('private')
        if kept_vars or self.metadata_vars:
            self.output_comment('\nSettable parameters and readable variables\n',
                                subsidiary=True)
        for var in kept_vars:
            self.writeln_hpp(self.TYPE_DOUBLE, self.code_name(var), self.STMT_END)
        # Generate Set & Get methods
        self.set_access('public')
        for var in kept_vars:
            if var.get_type() == VarTypes.Constant:
                # Generate a Set method
                self.output_method_start('Set_' + self.code_name(var, prefix=''),
                                         [self.TYPE_CONST_DOUBLE + 'value'],
                                         'void')
                self.open_block()
                self.writeln(self.code_name(var), self.EQ_ASSIGN, 'value;')
                self.close_block()
            # Generate Get method
            self.output_method_start('Get_' + self.code_name(var, prefix=''), [], self.TYPE_DOUBLE)
            self.open_block()
            self.writeln('return ', self.code_name(var), ';')
            self.close_block()
        
        # Methods associated with oxmeta annotated variables
        if self.use_metadata:
            for var in self.metadata_vars:
                if var.get_type() == VarTypes.Constant:
                    self.output_method_start('Get_' + var.oxmeta_name + '_constant', [], self.TYPE_DOUBLE)
                    self.open_block()
                    self.output_comment('Constant value given in CellML')
                    self.writeln('return ', var.initial_value, self.STMT_END)
                    self.close_block()
                    self.writeln()

    def output_top_boilerplate(self):
        """Output top boilerplate.
        
        This method outputs the constructor and destructor of the cell
        class, and also lookup table declarations and lookup methods.
        It also outputs a blank VerifyStateVariables method.
        """
        self.include_serialization = not (self.use_metadata or self.dynamically_loadable) # TODO: Implement
        self.check_time_units()
        # Check if we're generating a Backward Euler model
        if hasattr(self.model, u'solver_info') and \
               hasattr(self.model.solver_info, u'jacobian'):
            self.use_backward_euler = True
            # Find the size of the nonlinear system
            num_linear_odes = len(self.model.solver_info.xml_xpath(
                u'solver:linear_odes/m:math/m:apply'))
            self.nonlinear_system_size = len(self.state_vars) - 1 - num_linear_odes
            nonlinear_entries = self.model.solver_info.xml_xpath(
                u'solver:jacobian/solver:entry/@var_j')
            self.nonlinear_system_vars = \
                map(unicode, nonlinear_entries[:self.nonlinear_system_size])
        else:
            self.use_backward_euler = False
        # Start output
        self.output_includes()
        
        if self.use_backward_euler:
            # Keep the same signature as forward cell models, but note that the solver
            # isn't used
            solver1 = 'boost::shared_ptr<AbstractIvpOdeSolver> /* unused; should be empty */'
            solver2 = ''
            #solver1 = solver2 = ''
        else:
            solver1 = 'boost::shared_ptr<AbstractIvpOdeSolver> pSolver'
            solver2 = 'pSolver'

        if self.use_lookup_tables and self.separate_lut_class:
            self.output_lut_class()

        # Cell model class
        self.writeln_hpp('class ', self.class_name,
                         ' : public ', self.base_class_name)
        self.open_block(subsidiary=True)
        # Serialization
        if self.include_serialization:
            self.writeln_hpp('friend class boost::serialization::access;')
            self.writeln_hpp('template<class Archive>')
            self.writeln_hpp('void serialize(Archive & archive, const unsigned int version)')
            self.open_block(subsidiary=True)
            self.writeln_hpp('archive & boost::serialization::base_object<', self.base_class_name,
                             ' >(*this);')
            self.close_block(subsidiary=True)
        # Parameter declarations, and set & get methods (#666)
        param_vars = self.output_cell_parameters()
        # Constructor
        self.set_access('public')
        self.output_constructor([solver1, 'boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus'],
                                [solver2, len(self.state_vars), self.v_index, 'pIntracellularStimulus'])
        # Destructor
        self.output_method_start('~'+self.class_name, [], '')
        self.open_block()
        self.close_block()
        # Lookup table declarations & methods
        if self.use_lookup_tables:
            self.send_main_output_to_subsidiary()
            if self.separate_lut_class:
                if self.use_backward_euler:
                    self.output_lut_indices()
            else:
                self.output_lut_declarations()
                self.output_lut_row_lookup_memory()
                self.output_lut_methods()
            self.send_main_output_to_subsidiary(False)
        # Verify state variables method; empty at present
        self.output_method_start('VerifyStateVariables', [], 'void')
        self.writeln('{}\n')
        return
    
    def output_constructor(self, params, base_class_params):
        """Output a cell constructor.
        
        params is a list of constructor parameters, entries of which should be strings
        including both type and parameter name, which will be included verbatim in the
        generated code.
        
        base_class_params is a list of parameters to be supplied to the base class
        constructor.  Entries will be converted to strings.
        """
        self.output_method_start(self.class_name, params, '', access='public')
        self.writeln('    : ', self.base_class_name, '(')
        for i, param in enumerate(base_class_params):
            if i == len(base_class_params)-1: comma = ')'
            else: comma = ','
            self.writeln(param, comma, indent_offset=3)
        self.open_block()
        self.output_comment('Time units: ', self.free_vars[0].units, '\n')
        self.writeln('mpSystemInfo = OdeSystemInformation<',
                     self.class_name, '>::Instance();')
        self.writeln('Init();\n')
        #666 - initialise parameters
        for var in self.cell_parameters:
            if var.get_type() == VarTypes.Constant:
                self.writeln(self.code_name(var), self.EQ_ASSIGN, var.initial_value,
                             self.STMT_END)
        if self.use_lookup_tables and not self.separate_lut_class:
            self.output_lut_generation()
        self.close_block()
        return

    def output_lut_class(self):
        """Output a separate class for lookup tables.
        
        This will live entirely in the .cpp file."""
        # Lookup tables class
        self.writeln('class ', self.lt_class_name)
        self.writeln('{')
        self.writeln('public:')
        self.set_indent(1)
        # Method to get the table instance object
        self.writeln('static ', self.lt_class_name, '* Instance()')
        self.open_block()
        self.writeln('if (mpInstance == NULL)')
        self.writeln('{')
        self.writeln('mpInstance = new ', self.lt_class_name, ';',
                     indent_offset=1)
        self.writeln('}')
        self.writeln('return mpInstance;')
        self.close_block()
        # Table lookup methods
        self.output_lut_methods()
        # Make the class a singleton
        self.writeln('protected:', indent_level=0)
        self.writeln(self.lt_class_name, '(const ', self.lt_class_name,
                     '&);')
        self.writeln(self.lt_class_name, '& operator= (const ',
                     self.lt_class_name, '&);')
        # Constructor
        self.writeln(self.lt_class_name, '()')
        self.open_block()
        self.writeln('assert(mpInstance == NULL);')
        self.output_lut_generation()
        self.close_block()
        # Private data
        self.writeln('private:', indent_level=0)
        self.writeln('/** The single instance of the class */')
        self.writeln('static ', self.lt_class_name, ' *mpInstance;\n')
        if self.row_lookup_method:
            self.output_lut_row_lookup_memory()
        self.output_lut_declarations()
        # Close the class
        self.set_indent(0)
        self.writeln('};\n')
        # Initialise the instance pointer
        self.writeln(self.lt_class_name, '* ', self.lt_class_name,
                     '::mpInstance = NULL;')
        self.writeln()
        return

    def output_state_assignments(self, exclude_nonlinear=False,
                                 assign_rY=True,
                                 nodeset=None):
        """Output statements extracting state variables from their vector.

        If exclude_nonlinear is set to true, state variables appearing
        in the nonlinear system will not be included.

        If nodeset is given, only state variables appearing in nodeset
        will be included.
        """
        if assign_rY:
            self.writeln(self.TYPE_VECTOR_REF, ' rY = rGetStateVariables();')
        for i, var in enumerate(self.state_vars):
            if ( not exclude_nonlinear or 
                 self.code_name(var) not in self.nonlinear_system_vars) \
                 and (not nodeset or var in nodeset):
                if self.use_metadata and var.oxmeta_name:
                    value = self.modifier_call(var, self.vector_index('rY', i))
                else:
                    value = self.vector_index('rY', i)
                self.writeln(self.TYPE_DOUBLE, self.code_name(var),
                             self.EQ_ASSIGN, value, self.STMT_END)
                self.writeln(self.COMMENT_START, 'Units: ', var.units,
                             '; Initial value: ',
                             getattr(var, u'initial_value', 'Unknown'))
                #621 TODO: maybe convert if state var dimensions include time
        self.writeln()
        return
    
    def modifier_call(self, var, current_value):
        """Return code for a call to a modifier function for an oxmeta-annotated variable.
        
        The modifier function takes 2 parameters: the current value of the variable,
        and the current time.  It returns a modified value for the variable.
        """
        return ('this->' + var.oxmeta_name + '_modifier->calc(' +
                current_value + ', ' + self.code_name(self.free_vars[0]) + ')')
    
    def vector_index(self, vector, i):
        """Return code for accessing the i'th index of vector."""
        return vector + '[' + str(i) + ']'

    def output_nonlinear_state_assignments(self):
        """Output assignments for nonlinear state variables."""
        for i, varname in enumerate(self.nonlinear_system_vars):
            self.writeln('double ', varname, self.EQ_ASSIGN,
                         self.vector_index('rCurrentGuess', i), self.STMT_END)
            #621 TODO: maybe convert if state var dimensions include time
        self.writeln()
        return

    def output_equations(self, nodeset):
        """Output the mathematics described by nodeset.

        nodeset represents a subset of the assignments in the model.
        Output assignments in the order given by a topological sort,
        but only include those in nodeset.
        """
        for expr in (e for e in self.model.get_assignments() if e in nodeset):
            # Special-case the stimulus current
            if self.use_chaste_stimulus:
                if isinstance(expr, cellml_variable) and \
                        expr is self.doc._cml_config.i_stim_var:
                    self.writeln(self.TYPE_DOUBLE, nl=False)
                    self.write(self.code_name(expr))
                    self.write(self.EQ_ASSIGN)
                    #621: convert if free var is not in milliseconds
                    if self.conversion_factor:
                        conv = '(1.0/%.12g)*' % self.conversion_factor
                    else:
                        conv = ''
                    self.writeln('GetStimulus(', conv,
                                 self.code_name(self.free_vars[0]), ')',
                                 self.STMT_END, indent=False)
                elif not (isinstance(expr, mathml_apply) and
                          isinstance(expr.operator(), mathml_eq) and
                          isinstance(expr.eq.lhs, mathml_ci) and
                          expr.eq.lhs.variable is self.doc._cml_config.i_stim_var):
                    self.output_assignment(expr)
            else:
                self.output_assignment(expr)
        return

    def output_assignment(self, expr):
        """Output an assignment statement.

        Variables annotated with pe:keep should be stored as member
        variables (#666) so they can be read/set by users.  Hence if
        the variable assigned to is such a one, temporarily clear
        self.TYPE_DOUBLE and self.TYPE_CONST_DOUBLE before calling the
        base class method.
        """
        close_if = False #907
        if isinstance(expr, cellml_variable):
            assigned_var = expr
        else:
            if expr.eq.lhs.localName == 'ci':
                assigned_var = expr.eq.lhs.variable
            else:
                assigned_var = None # We don't store derivatives as members
                #907: Check if this is the derivative of the transmembrane
                # potential.
                if expr.eq.lhs.diff.dependent_variable == self.v_variable:
                    # Declare the variable for the derivative
                    self.writeln()
                    self.writeln(self.TYPE_DOUBLE, nl=False)
                    self.output_lhs(expr.eq.lhs)
                    self.writeln(self.STMT_END, indent=False)
                    # Fix to zero case
                    self.writeln('if (mSetVoltageDerivativeToZero)')
                    self.open_block()
                    self.writeln('', nl=False)
                    self.output_lhs(expr.eq.lhs)
                    self.writeln(self.EQ_ASSIGN, '0.0', self.STMT_END,
                                 indent=False)
                    self.close_block(blank_line=False)
                    self.writeln('else')
                    self.open_block()
                    close_if = True
        clear_type = (self.kept_vars_as_members and assigned_var and
                      assigned_var.pe_keep) or close_if
        if clear_type:
            self.TYPE_DOUBLE = self.TYPE_CONST_DOUBLE = ''
        if (self.use_metadata and assigned_var and assigned_var.oxmeta_name
            and assigned_var.get_type() == VarTypes.Constant):
            # "Constant" oxmeta-annotated parameters may be modified at run-time
            self.writeln(self.TYPE_DOUBLE, self.code_name(assigned_var), self.EQ_ASSIGN,
                         self.modifier_call(assigned_var, expr.initial_value), self.STMT_END)
        else:
            super(CellMLToChasteTranslator, self).output_assignment(expr)
        if clear_type:
            # Remove the instance attributes, thus reverting to the class members
            del self.TYPE_DOUBLE
            del self.TYPE_CONST_DOUBLE
        if close_if:
            self.close_block()
        return

    def output_mathematics(self):
        """Output the mathematics in this model.

        When backward Euler is used, we do so in 5 methods:
         * UpdateTransmembranePotential  does a forward Euler step for V
         * ComputeOneStepExceptVoltage  co-ordinates a backward Euler step
         * ComputeResidual and ComputeJacobian are used in the Newton iteration
         * GetIIonic returns the total ionic current

        For other solvers, only 2 methods are needed:
         * EvaluateYDerivatives computes the RHS of the ODE system
         * GetIIonic is as above
        """
        self.output_get_i_ionic()
        if self.use_backward_euler:
            self.output_backward_euler_mathematics()
        else:
            self.output_evaluate_y_derivatives()

    def output_get_i_ionic(self):
        """Output the GetIIonic method."""
        use_metadata = self.use_metadata
        self.use_metadata = False
        self.output_method_start('GetIIonic', [], self.TYPE_DOUBLE, access='public')
        self.open_block()
        # Output mathematics to calculate ionic current, using
        # solver_info.ionic_current.
        if hasattr(self.model, u'solver_info') and \
               hasattr(self.model.solver_info, u'ionic_current'):
            if not hasattr(self.model.solver_info.ionic_current, u'var'):
                raise ValueError('No ionic currents found; check your configuration file')
            nodes = map(lambda elt: self.varobj(unicode(elt)),
                        self.model.solver_info.ionic_current.var)
            nodeset = self.calculate_extended_dependencies(nodes)
            #print map(lambda v: v.fullname(), nodes)
            #print filter(lambda p: p[2]>0, map(debugexpr, nodeset))
            self.output_state_assignments(nodeset=nodeset)
            if self.use_lookup_tables:
                self.output_table_index_generation(
                    indexes_as_member=self.use_backward_euler)
            self.output_equations(nodeset)
            self.writeln()
            #621: check units of the ionic current
            conv = self.ionic_current_units_conversion_factor(
                self.get_var_units(nodes[0]))
            if conv == 1:
                conv = ''
            else:
                conv = str(conv) + ' * '
            self.writeln('return ', conv, '(', nl=False)
            plus = False
            for varelt in self.model.solver_info.ionic_current.var:
                if plus: self.write('+')
                else: plus = True
                self.output_variable(varelt)
            self.writeln(')', self.STMT_END, indent=False)
        else:
            self.writeln('return 0.0;')
        self.close_block()
        self.use_metadata = use_metadata

    def output_evaluate_y_derivatives(self, method_name='EvaluateYDerivatives'):
        """Output the EvaluateYDerivatives method."""
        # Work out what equations are needed to compute the derivatives
        nodes = map(lambda v: (v, self.free_vars[0]), self.state_vars)
        if self.use_chaste_stimulus:
            i_stim = self.doc._cml_config.i_stim_var
            nodeset = self.calculate_extended_dependencies(
                nodes, prune=[i_stim])
            nodeset.add(i_stim)
        else:
            nodeset = self.calculate_extended_dependencies(nodes)
        # Start code output
        self.output_method_start(method_name,
                                 [self.TYPE_DOUBLE + self.code_name(self.free_vars[0]),
                                  'const ' + self.TYPE_VECTOR_REF + ' rY',
                                  self.TYPE_VECTOR_REF + ' rDY'],
                                 'void', access='public')
        self.open_block()
        self.output_comment('Inputs:')
        self.output_comment('Time units: ', self.free_vars[0].units)
        #621: convert if free var is not in milliseconds
        if self.conversion_factor:
            self.writeln(self.code_name(self.free_vars[0]), ' *= ',
                         self.conversion_factor, self.STMT_END)
        self.output_state_assignments(assign_rY=False, nodeset=nodeset)
        self.writeln()
        if self.use_lookup_tables:
            self.output_table_index_generation(
                indexes_as_member=self.use_backward_euler)
        self.writeln(self.COMMENT_START, 'Mathematics')
        self.output_equations(nodeset)
        self.writeln()
        for i, var in enumerate(self.state_vars):
            deriv_assign = self.vector_index('rDY', i) + self.EQ_ASSIGN
            if self.conversion_factor:
                #621: convert if free var is not in milliseconds
                self.writeln(deriv_assign, self.conversion_factor,
                             '*', self.code_name(var, True), self.STMT_END)
            else:
                self.writeln(deriv_assign, self.code_name(var, True),
                             self.STMT_END)
        self.close_block()
        return

    def output_backward_euler_mathematics(self):
        """Output the mathematics methods used in a backward Euler cell.

        Outputs ComputeResidual, ComputeJacobian,
        UpdateTransmembranePotential and ComputeOneStepExceptVoltage.
        """
        # Residual
        ##########
        argsize = '[' + str(self.nonlinear_system_size) + ']'
        self.output_method_start('ComputeResidual',
                                 [self.TYPE_CONST_DOUBLE + 'rCurrentGuess' + argsize,
                                  self.TYPE_DOUBLE + 'rResidual' + argsize],
                                 'void', access='public')
        self.open_block()
        # Output mathematics for computing du/dt for each nonlinear state var u
        nodes = map(lambda u: (self.varobj(u), self.free_vars[0]),
                    self.nonlinear_system_vars)
        nodeset = self.calculate_extended_dependencies(nodes)
        self.output_state_assignments(exclude_nonlinear=True, nodeset=nodeset)
        self.output_nonlinear_state_assignments()
        self.output_equations(nodeset)
        self.writeln()
        # Fill in residual
        for i, var in enumerate(self.state_vars):
            try:
                j = self.nonlinear_system_vars.index(self.code_name(var))
            except ValueError:
                j = -1
            if j != -1:
                if self.conversion_factor:
                    #621: convert if free var is not in milliseconds
                    self.writeln('rResidual[', j, '] = rCurrentGuess[', j,
                                 '] - rY[', i, '] - mDt*',
                                 self.conversion_factor, '*',
                                 self.code_name(var, True), self.STMT_END)
                else:
                    self.writeln('rResidual[', j, '] = rCurrentGuess[', j,
                                 '] - rY[', i, '] - mDt*',
                                 self.code_name(var, ode=True), self.STMT_END)
        self.close_block()
        # Jacobian
        ##########
        self.output_method_start('ComputeJacobian',
                                 [self.TYPE_CONST_DOUBLE + 'rCurrentGuess' + argsize,
                                  self.TYPE_DOUBLE + 'rJacobian' + argsize + argsize],
                                 'void', access='public')
        self.open_block()
        # Mathematics that the Jacobian depends on
        used_vars = set()
        for entry in self.model.solver_info.jacobian.entry:
            used_vars.update(self._vars_in(entry.math))
        nodeset = self.calculate_extended_dependencies(used_vars)
        self.output_state_assignments(exclude_nonlinear=True, nodeset=nodeset)
        self.output_nonlinear_state_assignments()
        if self.conversion_factor:
            self.writeln('const double dt = ', self.conversion_factor,
                         ' * mDt;\n');
        else:
            self.writeln('const double dt = mDt;\n')
        self.output_equations(nodeset)
        self.writeln()
        # Jacobian entries
        for entry in self.model.solver_info.jacobian.entry:
            var_i, var_j = entry.var_i, entry.var_j
            i = self.nonlinear_system_vars.index(var_i)
            j = self.nonlinear_system_vars.index(var_j)
            self.writeln('rJacobian[', i, '][', j, '] = ', nl=False)
            if hasattr(entry.math, u'apply'):
                self.output_expr(entry.math.apply, False)
            else:
                self.output_expr(entry.math.cn, False)
            self.writeln(self.STMT_END, indent=False)
        self.close_block()
        # The other methods are protected
        self.writeln('protected:', indent_offset=-1)
        # UpdateTransmembranePotential
        ##############################
        self.output_method_start('UpdateTransmembranePotential',
                                 [self.TYPE_DOUBLE + self.code_name(self.free_vars[0])],
                                 'void', access='public')
        self.open_block()
        self.output_comment('Time units: ', self.free_vars[0].units)
        #621: convert if free var is not in milliseconds
        if self.conversion_factor:
            self.writeln(self.code_name(self.free_vars[0]), ' *= ',
                         self.conversion_factor, self.STMT_END)
        # Output mathematics to compute dV/dt
        nodes = [(self.state_vars[self.v_index], self.free_vars[0])]
        nodeset = self.calculate_extended_dependencies(nodes)
        self.output_state_assignments(nodeset=nodeset)
        if self.use_lookup_tables:
            self.output_table_index_generation(indexes_as_member=True)
        self.output_equations(nodeset)
        # Update V
        self.writeln()
        if self.conversion_factor:
            #621: convert if free var is not in milliseconds
            self.writeln('rY[', self.v_index, '] += mDt * ',
                         self.conversion_factor, '*',
                         self.code_name(self.state_vars[self.v_index],
                                        ode=True), self.STMT_END)
        else:
            self.writeln('rY[', self.v_index, '] += mDt * ',
                         self.code_name(self.state_vars[self.v_index],
                                        ode=True), self.STMT_END)
        self.close_block()
        # ComputeOneStepExceptVoltage
        ######################
        self.output_method_start('ComputeOneStepExceptVoltage',
                                 [self.TYPE_DOUBLE + self.code_name(self.free_vars[0])],
                                 'void', access='public')
        self.open_block()
        self.writeln(self.COMMENT_START, 'Time units: ',
                     self.free_vars[0].units)
        #621: convert if free var is not in milliseconds
        if self.conversion_factor:
            self.writeln(self.code_name(self.free_vars[0]), ' *= ',
                         self.conversion_factor, self.STMT_END)
        # Output mathematics to update linear state variables, using
        # solver_info.linear_odes.  Need to analyse maths to determine
        # which elements correspond to g and h, and what the linear vars
        # are.  Also need to use output_equations for variables used in
        # solver_info.linear_odes.
        linear_vars, ghs = [], []
        used_vars = set()
        for ode in self.model.solver_info.linear_odes.math.apply:
            varname = unicode(ode.apply.ci)
            g = ode.apply[1].operands().next()
            hu = list(ode.apply[1].operands())[1]
            h = hu.operands().next()
            linear_vars.append(self.varobj(varname))
            ghs.append((g, h))
            used_vars.update(self._vars_in(g))
            used_vars.update(self._vars_in(h))
        # Output required equations for used variables
        nodeset = self.calculate_extended_dependencies(used_vars)
        self.output_state_assignments(nodeset=nodeset)
        if self.use_lookup_tables:
            self.output_table_index_generation(indexes_as_member=True)
        self.output_equations(nodeset)
        # Output g and h calculations
        self.writeln()
        for i, gh in enumerate(ghs):
            g, h = gh
            self.writeln('const double _g_', i, ' = ', nl=False)
            self.output_expr(g, False)
            self.writeln(self.STMT_END, indent=False)
            self.writeln('const double _h_', i, ' = ', nl=False)
            self.output_expr(h, False)
            self.writeln(self.STMT_END, indent=False)
        # Update state variables:
        #   rY[i] = (rY[i] + _g_j*mDt) / (1 - _h_j*mDt)
        self.writeln()
        if self.conversion_factor:
            self.writeln(self.TYPE_CONST_DOUBLE, 'dt = mDt*',
                         self.conversion_factor, self.STMT_END)
        for i, var in enumerate(self.state_vars):
            try:
                j = linear_vars.index(var)
            except ValueError:
                j = -1
            if j != -1:
                if self.conversion_factor:
                    #621 TODO: convert if free var is not in milliseconds
                    self.writeln('rY[', i, '] = (rY[', i, '] + _g_', j,
                                 '*dt) / (1 - _h_', j, '*dt);')
                else:
                    self.writeln('rY[', i, '] = (rY[', i, '] + _g_', j,
                                 '*mDt) / (1 - _h_', j, '*mDt);')
        # Set up the Newton iteration
        self.writeln()
        self.writeln('double _guess[', self.nonlinear_system_size,
                     '] = {', nl=False)
        comma = False
        idx_map = [0] * self.nonlinear_system_size
        for i, var in enumerate(self.state_vars):
            try:
                j = self.nonlinear_system_vars.index(self.code_name(var))
                idx_map[j] = i
            except ValueError:
                pass
        for i in idx_map:
            if comma: self.write(',')
            else: comma = True
            self.write('rY[', i, ']')
        self.writeln('};', indent=False)
        # Solve
        self.writeln('CardiacNewtonSolver<', self.nonlinear_system_size,
                     '> *_solver = CardiacNewtonSolver<',
                     self.nonlinear_system_size, '>::Instance();')
        self.writeln('_solver->Solve(*this, _guess);')
        # Update state
        for j, i in enumerate(idx_map):
            self.writeln('rY[', i, '] = _guess[', j, '];')
        self.close_block()
        return

    def output_bottom_boilerplate(self):
        """Output bottom boilerplate.

        End class definition, output ODE system information (to .cpp) and
        serialization code (to .hpp), and end the file.
        """
        # End main class
        self.set_indent(offset=-1)
        self.writeln_hpp('};\n\n')
        # ODE system information
        self.writeln('template<>')
        self.writeln('void OdeSystemInformation<', self.class_name,
                     '>::Initialise(void)')
        self.open_block()
        self.output_comment('Time units: ', self.free_vars[0].units, '\n')
        for var in self.state_vars:
            if self.use_metadata and var.oxmeta_name:
                self.writeln('this->mVariableNames.push_back("', var.oxmeta_name, '");')    
            else:
                self.writeln('this->mVariableNames.push_back("', var.name, '");')
            self.writeln('this->mVariableUnits.push_back("', var.units, '");')
            init_val = getattr(var, u'initial_value', None)
            if init_val is None:
                init_comm = ' // Value not given in model'
                # Don't want compiler error, but shouldn't be a real number
                init_val = 'NAN'
            else:
                init_comm = ''
            self.writeln('this->mInitialConditions.push_back(', init_val, ');',
                       init_comm, '\n')
        self.writeln('this->mInitialised = true;')
        self.close_block()
        self.writeln()
        # Serialization
        if self.include_serialization:
            self.output_comment('Needs to be included last', subsidiary=True)
            self.writeln_hpp('#include "TemplatedExport.hpp"')
            self.writeln_hpp('CHASTE_CLASS_EXPORT(', self.class_name, ')')
            self.writeln_hpp()
            self.writeln_hpp('namespace boost')
            self.open_block(subsidiary=True)
            self.writeln_hpp('namespace serialization')
            self.open_block(subsidiary=True)
            # Save
            self.writeln_hpp('template<class Archive>')
            self.writeln_hpp('inline void save_construct_data(')
            self.writeln_hpp('Archive & ar, const ', self.class_name,
                             ' * t, const unsigned int fileVersion)',
                             indent_offset=1)
            self.open_block(subsidiary=True)
            self.writeln_hpp('const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();')
            self.writeln_hpp('const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();')
            self.writeln_hpp('ar << p_solver;')
            self.writeln_hpp('ar << p_stimulus;')
            self.close_block(subsidiary=True)
            # Load
            self.writeln_hpp('template<class Archive>')
            self.writeln_hpp('inline void load_construct_data(')
            self.writeln_hpp('Archive & ar, ', self.class_name,
                             ' * t, const unsigned int fileVersion)',
                             indent_offset=1)
            self.open_block(subsidiary=True)
            self.writeln_hpp('boost::shared_ptr<AbstractIvpOdeSolver> p_solver;')
            self.writeln_hpp('boost::shared_ptr<AbstractStimulusFunction> p_stimulus;')
            self.writeln_hpp('ar >> p_solver;')
            self.writeln_hpp('ar >> p_stimulus;')
            self.writeln_hpp('::new(t)', self.class_name, '(p_solver, p_stimulus);')
            self.close_block(subsidiary=True)
            self.close_block(subsidiary=True)
            self.close_block(subsidiary=True)
        if self.dynamically_loadable:
            # Write the C function to create instances of this cell model
            self.writeln('extern "C"')
            self.open_block()
            self.writeln('AbstractCardiacCell* MakeCardiacCell(')
            self.writeln('boost::shared_ptr<AbstractIvpOdeSolver> pSolver,', indent_offset=2)
            self.writeln('boost::shared_ptr<AbstractStimulusFunction> pStimulus)', indent_offset=2)
            self.open_block()
            self.writeln('return new ', self.class_name, '(pSolver, pStimulus);')
            self.close_block()
            self.close_block()
        # End file
        self.writeln_hpp('#endif // ', self.include_guard)
        return

    def output_lhs(self, expr):
        """Output the left hand side of an assignment expression."""
        if expr.localName == 'ci':
            self.output_variable(expr)
        elif expr.operator().localName == 'diff':
            ci_elt = expr.operands().next()
            self.output_variable(ci_elt, ode=True)
        return

    def output_variable(self, ci_elt, ode=False):
        """Output a ci element, i.e. a variable lookup."""
        if hasattr(ci_elt, '_cml_variable') and ci_elt._cml_variable:
            self.write(self.code_name(ci_elt.variable, ode=ode))
        else:
            # This ci element is in the solver_info section, thus
            # doesn't have all the extra annotations.  It is a fully
            # qualified name though.
            prefix = ['var_', 'd_dt_'][ode]
            varname = unicode(ci_elt)
            if varname[0] == '(':
                # (compname,varname)
                cname, vname = varname[1:-1].split(u',')
                if self.single_component:
                    varname = vname
                else:
                    varname = cname + '__' + vname
            elif varname == u'delta_t':
                # Special case for the timestep in ComputeJacobian
                prefix = ''
                varname = 'dt'
            else:
                # var_cname__vname
                varname = varname[4:]
            self.write(prefix + varname)
        return


class CellMLToCvodeTranslator(CellMLToChasteTranslator):
    """Translate a CellML model to C++ code for use with Chaste+CVODE."""
    
    TYPE_VECTOR_REF = 'N_Vector' # CVODE's vector is actually a pointer type
    
    def vector_index(self, vector, i):
        """Return code for accessing the i'th index of vector."""
        return 'NV_Ith_S(' + vector + ', ' + str(i) + ')'

    def output_top_boilerplate(self):
        """Output top boilerplate code.

        This method outputs #includes, and the start of the cell class
        with constructor, destructor, and LT methods.
        """
        # CVODE is optional in Chaste
        self.writeln("#ifdef CHASTE_CVODE")
        self.writeln_hpp("#ifdef CHASTE_CVODE")
        self.include_serialization = False
        self.use_backward_euler = False
        self.check_time_units()
        self.output_includes(base_class='AbstractCvodeCell')
        # Separate class for lookup tables?
        if self.use_lookup_tables and self.separate_lut_class:
            self.output_lut_class()
        # Start cell model class
        self.writeln_hpp('class ', self.class_name, ' : public ', self.base_class_name)
        self.open_block(subsidiary=True)
        # Parameter declarations, and set & get methods (#666)
        self.output_cell_parameters()
        # Constructor
        self.output_constructor(['boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver /* unused; should be empty */',
                                 'boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus'],
                                ['pOdeSolver', len(self.state_vars), self.v_index, 'pIntracellularStimulus'])
        # Destructor
        self.output_method_start('~'+self.class_name, [], '', access='public')
        self.open_block()
        self.close_block()
        # Lookup table declarations & methods
        if self.use_lookup_tables and not self.separate_lut_class:
            self.send_main_output_to_subsidiary()
            self.output_lut_declarations()
            self.output_lut_row_lookup_memory()
            self.output_lut_methods()
            self.send_main_output_to_subsidiary(False)
        # Verify state variables method; empty at present
        self.output_method_start('VerifyStateVariables', [], 'void', access='public')
        self.writeln('{}\n')
        return

    def output_mathematics(self):
        """Output the mathematics in this model.
        
        Two methods are needed:
         * EvaluateRhs computes the RHS of the ODE system
         * GetIIonic returns the total ionic current
        """
        self.output_get_i_ionic()
        self.output_evaluate_y_derivatives(method_name='EvaluateRhs')
    
    def output_bottom_boilerplate(self):
        """Call superclass method, then end the CHASTE_CVODE guard."""
        super(CellMLToCvodeTranslator, self).output_bottom_boilerplate()
        # CVODE is optional in Chaste
        self.writeln("#endif // CHASTE_CVODE")
        self.writeln_hpp("#endif // CHASTE_CVODE")

class CellMLToMapleTranslator(CellMLTranslator):
    """Translate a CellML model to Maple code."""

    # Language tokens that differ from the default
    EQ_ASSIGN = ' := '   # Assignment operator
    COMMENT_START = '# ' # Start of a 1 line comment
    # Types are determined automatically by Maple
    TYPE_DOUBLE = ''
    TYPE_CONST_DOUBLE = ''
    # Some constants are different
    PI = 'Pi'
    E = 'exp(1)'

    def __init__(self, omit_constants=False, compute_full_jacobian=False,
                 **kwargs):
        """Create a Maple translator.

        If omit_constants is set to true, assignments will not be
        generated for constant variables.  This should be used if
        these values will be altered at runtime, in order to prevent
        derivatives being calculated incorrectly.

        Set compute_full_jacobian to True to make Maple compute the
        Jacobian of the whole ODE system, rather than just the
        nonlinear portion.

        Other keyword arguments are all passed to the base class.
        """
        super(CellMLToMapleTranslator, self).__init__(**kwargs)
        # Maple translation doesn't support lookup tables
        self.use_lookup_tables = False
        # Translation parameters
        self.omit_constants = omit_constants
        self.compute_full_jacobian = compute_full_jacobian
        # Update some function names
        self.function_map = {}
        self.function_map.update(CellMLTranslator.function_map)
        del self.function_map['power']
        self.function_map.update(
            {'abs': 'abs', 'ln': 'ln', 'not': 'not',
             'sec': 'sec', 'csc': 'csc', 'cot': 'cot',
             'sech': 'sech', 'csch': 'csch', 'coth': 'coth',
             'arcsin': 'arcsin', 'arccos': 'arccos', 'arctan': 'arctan',
             'arcsec': 'arcsec', 'arccsc': 'arccsc', 'arccot': 'arccot',
             'arcsinh': 'arcsinh', 'arccosh': 'arccosh', 'arctanh': 'arctanh',
             'arcsech': 'arcsech', 'arccsch': 'arccsch', 'arccoth': 'arccoth'})
        self.recip_trig = {}
        self.nary_ops = {}
        self.nary_ops.update(CellMLTranslator.nary_ops)
        self.nary_ops.update({'and': 'and', 'or': 'or'})
        self.binary_ops = {}
        self.binary_ops.update(CellMLTranslator.binary_ops)
        self.binary_ops.update({'xor': 'xor', 'eq': '=', 'neq': '<>',
                                'power': '^'})
        self.special_roots = {}

    def output_file_name(self, model_filename):
        """Generate a name for our output file, based on the input file."""
        return os.path.splitext(model_filename)[0] + '.mpl'

    def output_top_boilerplate(self):
        """Output top boilerplate."""
        self.writeln('# Model: ', self.model.name)
        self.output_comment(version_comment(self.add_timestamp))
        self.writeln()
        self.writeln('interface(prettyprint=0);\n')
        return

    def output_bottom_boilerplate(self):
        """Output bottom boilerplate."""
        self.output_comment('\nJacobian calculation\n')
        if self.compute_full_jacobian:
            # Jacobian calculation for the whole ODE system
            # i.e. each df_i/du_j
            for i, var_i in enumerate(self.state_vars):
                for j, var_j in enumerate(self.state_vars):
                    self.writeln('print("--', self.code_name(var_i), '/',
                                 self.code_name(var_j), '--");')
                    self.writeln('diff(', self.code_name(var_i, ode=True), ', ',
                                 self.code_name(var_j), ');')
        elif hasattr(self.model, '_cml_nonlinear_system_variables'):
            # Jacobian calculation for Jon Whiteley's algorithm
            vars_text = self.model._cml_nonlinear_system_variables
            if type(vars_text) == type([]):
                var_objs = self.model._cml_nonlinear_system_variables
            else:
                # Get variable objects from names
                varnames = map(lambda s: s.split(','), vars_text.split(';'))
                var_objs = map(lambda (c, v):
                               self.model.get_variable_by_name(c, v),
                               varnames)
            # Output the Newton iteration expression for each variable
            for var in var_objs:
                self.writeln('g_', self.code_name(var), self.EQ_ASSIGN,
                             self.code_name(var), ' - ',
                             self.code_name(var), '_old - delta_t*',
                             self.code_name(var, ode=True), self.STMT_END)
            # Output the derivative calculations
            for var_i in var_objs:
                for var_j in var_objs:
                    self.writeln('print("--', self.code_name(var_i), '/',
                                 self.code_name(var_j), '--");')
                    self.writeln('diff(g_', self.code_name(var_i), ', ',
                                 self.code_name(var_j), ');')
        # Tell Maple to quit when done
        self.writeln()
        self.writeln('quit;')
        return

    def output_assignment(self, expr):
        """Output an assignment expression.

        Optionally, if this is an assignment of a constant, don't output,
        so that differentation doesn't optimise expressions away.
        """
        if isinstance(expr, cellml_variable):
            # This may be the assignment of a mapped variable, or a constant
            t = expr.get_type()
            if t == VarTypes.Mapped:
                self.writeln(self.TYPE_DOUBLE, self.code_name(expr),
                             self.EQ_ASSIGN,
                             self.code_name(expr.get_source_variable()),
                             self.STMT_END)
            elif t == VarTypes.Constant and not self.omit_constants:
                self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(expr),
                             self.EQ_ASSIGN, nl=False)
                self.output_number(expr.initial_value)
                self.writeln(self.STMT_END, indent=False)
        else:
            # This is a mathematical expression
            self.writeln(self.TYPE_DOUBLE, nl=False)
            opers = expr.operands()
            self.output_lhs(opers.next())
            self.write(self.EQ_ASSIGN)
            self.output_expr(opers.next(), False)
            self.writeln(self.STMT_END, indent=False)
        return

    def output_number(self, expr):
        """Output the plain number expr.
        
        With Maple, there is no need to make all constants parse as
        doubles to avoid problems with integer division or numbers too
        large for the int type.
        
        Negative numbers will be prefixed by a space to avoid unwanted
        decrement operations.
        """
        n = self.eval_number(expr)
        num = "%.12g" % n
        if num[0] == '-':
            num = ' ' + num
        self.write(num)

    def output_root(self, expr, paren):
        """Output a root taken to some degree.

        If a degree qualifier element is not provided, uses default 2.
        """
        if hasattr(expr, u'degree'):
            # A degree is given.  Compute x^(1/b)
            self.open_paren(paren)
            self.output_expr(expr.operands().next(), True)
            self.write('^(1/')
            self.output_expr(expr.degree, True)
            self.write(')')
            self.close_paren(paren)
        else:
            # Compute square root
            self.output_function('sqrt', expr.operands(), paren)

    def output_log(self, expr, paren):
        """Output a logarithm to the given base, which defaults to base 10."""
        if hasattr(expr, u'logbase'):
            # A base is provided.  Use the log[b](x) function.
            self.write('log[')
            self.output_expr(expr.logbase, False)
            self.write(']')
            self.output_function('', expr.operands(), paren)
        else:
            # Use base 10
            self.output_function('log10', expr.operands(), paren)

    def output_piecewise(self, expr, paren):
        """Output the piecewise expression expr.

        We use the if operator.
        """
        num_ifs = 0
        for piece in getattr(expr, u'piece', []):
            num_ifs += 1
            self.write('`if`(')
            self.output_expr(child_i(piece, 2), False) # Condition
            self.write(',')
            self.output_expr(child_i(piece, 1), False) # Result
            self.write(',')
        if hasattr(expr, u'otherwise'):
            self.output_expr(child_i(expr.otherwise, 1), paren) # Default case
        else:
            self.write('FAIL') # If this is hit, things get ugly
        for i in range(num_ifs):
            self.close_paren(True)

class CellMLToHaskellTranslator(CellMLTranslator):
    """Translate a CellML model to a Haskell version.

    This does not produce a 'runnable' version of the model, but
    rather an equivalent model in effectively an abstract syntax,
    which can then be interpreted by a suitable interpreter.
    
    This allows us to more easily specify an operational semantics for
    CellML, without having to worry about the XML issues in the
    interpreter itself.
    """

    STMT_END = ''
    COMMENT_START = '-- '
    TYPE_DOUBLE = ''
    TYPE_CONST_DOUBLE = ''
    E = '(exp 1)'
    PI = 'pi'
    TRUE = '(Bool True)'
    FALSE = '(Bool False)'

    def __init__(self, **kwargs):
        """Create a Haskell translator.

        Keyword arguments are all passed to the base class.
        """
        super(CellMLToHaskellTranslator, self).__init__(**kwargs)
        # We don't use lookup tables in Haskell code
        self.use_lookup_tables = False
        return

    def output_file_name(self, model_filename):
        """Generate a name for our output file, based on the input file."""
        return os.path.splitext(model_filename)[0] + '.hs'

    def stringify(self, s):
        """Quote a string."""
        return '"' + s + '"'

    def code_name(self, var, ode=False, full=False):
        """
        Return the full name of var in a form suitable for inclusion in a
        source file.
        
        The functionality of ode=True is implemented in output_apply
        rather than here, so this parameter must be False.

        If full is True, include the name of the owning component.

        If PE has been performed (there is only 1 component, and variable
        names have been munged) then transform the munged name to Haskell
        munged form.
        """
        if ode:
            raise NotImplemented # Never used; see output_apply.
        if self.single_component and '__' in var.name:
            name = var.name.replace('__', ',')
        elif full:
            name = var.xml_parent.name + ',' + var.name
        else:
            name = var.name
        return self.stringify(name)

    def output_top_boilerplate(self):
        """Output top boilerplate.

        Outputs the imports and model-level units definitions.
        Components and connections are output by output_mathematics.
        """
        self.class_name = self.class_name.lower()
        self.module_name = self.class_name[0].upper() + self.class_name[1:]
        self.writeln('module ', self.module_name, ' where')
        self.writeln('-- Model: ', self.model.name)
        self.output_comment(version_comment(self.add_timestamp))
        self.writeln()
        self.writeln('import CellML')
        self.writeln('import Units')
        self.writeln('import Environment')
        self.writeln()
        # Model definition
        self.writeln(self.class_name, ' = Model ',
                     self.stringify(self.class_name),
                     ' units components connections')
        self.writeln('  where')
        self.set_indent(offset=1)
        # Model-level units definitions
        self.writeln('units =')
        self.output_units(self.model)
        return

    def output_unit(self, udef):
        """Output a single units definition, recursively."""
        def prefix(uref):
            """Return a prefix of a units reference,
            as an integer to which 10 may be raised."""
            prefix = getattr(uref, u'prefix_', 0)
            if prefix in uref.SI_PREFIXES:
                prefix = uref.SI_PREFIXES[prefix]
            else:
                prefix = int(prefix)
            return prefix
        def num(n):
            """Wrap numbers for Haskell."""
            if n < 0:
                return "(" + str(n) + ")"
            else:
                return n

        self.write('(')
        if udef.is_base_unit():
            # Base unit
            self.write('BaseUnits ', self.stringify(udef.name))
        elif udef.is_simple():
            # Simple units
            self.write('SimpleUnits ', num(udef.get_multiplier()), ' ',
                       num(prefix(udef.unit)), ' ')
            self.output_unit(udef.unit.get_units_element())
            self.write(' ', num(udef.get_offset()))
        else:
            # Complex units
            self.write('ComplexUnits [')
            uref_comma = False
            for uref in udef.unit:
                if uref_comma: self.write(',')
                else: uref_comma = True
                self.write('Unit ', num(uref.get_multiplier()), ' ',
                           num(prefix(uref)), ' ')
                self.output_unit(uref.get_units_element())
                self.write(' ', num(uref.get_exponent()))
            self.write(']')
        self.write(')')
        return

    def output_units(self, units_parent):
        """Output the units definitions in this model or component."""
        self.open_list()
        comma, output = False, False
        for udef in getattr(units_parent, u'units', []):
            output = True
            if comma: self.writeln(',', indent=False)
            else: comma = True
            # Output a single definition
            self.writeln('UDef ', self.stringify(udef.name), nl=False)
            self.output_unit(udef)
        if output:
            self.writeln('', indent=False)
        self.close_list()
        return

    def output_mathematics(self):
        """Output the mathematics in this model.

        This method outputs the components and connections."""
        # Components
        self.writeln('components =')
        self.open_list()
        comma = False
        for comp in getattr(self.model, u'component', []):
            if comma: self.writeln(',')
            else: comma = True
            self.output_component(comp)
        self.writeln('')
        self.close_list()
        # Connections
        self.writeln('connections =')
        self.open_list()
        comma = False
        for var in (v for v in self.model.get_assignments()
                    if isinstance(v, cellml_variable)
                    if v.get_type() == VarTypes.Mapped):
            if comma: self.writeln(',', indent=False)
            else: comma = True
            self.output_connection(var)
        self.writeln('', indent=False)
        self.close_list()
        return

    def output_component(self, comp):
        """Output a single component."""
        self.writeln('MkComp ', self.stringify(comp.name))
        self.output_units(comp)
        # Variable declarations, associating units with var names
        self.open_list()
        comma = False
        for var in getattr(comp, u'variable', []):
            if comma: self.writeln(',')
            else: comma = True
            self.writeln('VarDecl ', self.code_name(var), ' ',
                         self.stringify(var.units))
        self.close_list()
        # And now the mathematics
        self.open_list()
        comma = False
        # Constants
        for var in (v for v in getattr(comp, u'variable', [])
                    if v.get_type() == VarTypes.Constant):
            if comma: self.writeln(',')
            else: comma = True
            self.output_assignment(var)
        # Expressions
        for math in getattr(comp, u'math', []):
            for expr in getattr(math, u'apply', []):
                if comma: self.writeln(',')
                else: comma = True
                self.output_assignment(expr)
        self.close_list()
        return

    def output_connection(self, conn):
        """Output a single connection."""
        to_var = conn
        from_var = conn.get_source_variable()
        self.writeln('VarMap', nl=False)
        self.write(' (', self.stringify(from_var.xml_parent.name), ',',
                   self.stringify(from_var.name), ')')
        self.write(' (', self.stringify(to_var.xml_parent.name), ',',
                   self.stringify(to_var.name), ')')
        return

    def output_bottom_boilerplate(self):
        """Output bottom boilerplate."""
        self.set_indent(offset=-1)
        self.writeln()
        self.output_comment('Evaluate derivatives at the start of time.')
        self.writeln()
        # Initial environment
        self.writeln('initial_environment :: Env')
        self.writeln('initial_environment = make_env')
        self.open_list()
        self.writeln('  (Var ', self.code_name(self.free_vars[0], full=True),
                     ', Val (Number 0))')
        for sv in self.state_vars:
            self.writeln(', (Var ', self.code_name(sv, full=True),
                         ', Val ', nl=False)
            self.output_number(sv.initial_value, as_value=True)
            self.writeln(')', indent=False)
        self.close_list()
        self.writeln('results = run_cellml ', self.class_name,
                     ' initial_environment')
        self.writeln()
        # Dynamic environment for PE
        self.writeln('dynamic_environment :: Env')
        self.writeln('dynamic_environment = foldr def initial_environment')
        self.open_list()
        # Include all variables marked as pe:keep
        comma = False
        for comp in getattr(self.model, u'component', []):
            for var in getattr(comp, u'variable', []):
                if var.pe_keep:
                    self.writeln([' ', ','][comma], ' (Var ',
                                 self.code_name(var, full=True),
                                 ', Val DynamicMarker)')
                    if not comma: comma = True
        self.close_list()
        self.writeln('where def (k,v) env = define env k v', indent_offset=1)
        self.writeln('pe_results = reduce_and_run_cellml ', self.class_name,
                     ' dynamic_environment')
        return

    def open_list(self):
        """Open a multi-line list."""
        self.set_indent(offset=1)
        self.writeln('[')
        self.set_indent(offset=1)
        return

    def close_list(self):
        """Close a multi-line list."""
        self.set_indent(offset=-1)
        self.writeln(']')
        self.set_indent(offset=-1)
        return

    def output_assignment(self, expr):
        """Output an assignment expression."""
        if isinstance(expr, cellml_variable):
            # Assignment of a constant
            self.writeln('Assign (Var ', self.code_name(expr), ') ', nl=False)
            self.output_number(expr.initial_value, units=expr.units)
            self.writeln(indent=False)
        else:
            # This is a mathematical expression
            opers = expr.operands()
            self.writeln('Assign ', nl=False)
            self.output_lhs(opers.next())
            self.write(' ')
            self.output_expr(opers.next(), True)
            self.writeln(indent=False)
        return

    def output_lhs(self, expr):
        """Output the left hand side of an assignment expression."""
        if expr.localName == 'ci':
            self.output_variable(expr, lhs=True)
        elif expr.operator().localName == 'diff':
            v1 = expr.operator().dependent_variable
            v2 = expr.operator().independent_variable
            self.write('(Ode ', self.code_name(v1), ' ', self.code_name(v2),
                       ')')
        return

    def output_variable(self, ci_elt, lhs=False):
        """Output a ci element, i.e. a variable lookup."""
        type_constructor = ['Variable', 'Var'][lhs]
        self.write('(', type_constructor, ' ',
                   self.code_name(ci_elt.variable), ')')
        return
    
    def output_number(self, expr, as_value=False, units=None):
        """Output the plain number expr.
        
        With Haskell there is no need to force numbers to parse as doubles.
        We do need to bracket negative numbers.
        """
        n = self.eval_number(expr)
        num = "%.12g" % n
        if num[0] == '-':
            num = "(" + num + ")"
        tc = ['Num', 'Number'][as_value]
        if not as_value:
            if units is None:
                units = getattr(expr, u'units', '')
            units = '(Left ' + self.stringify(units) + ')'
        else:
            units = ''
        self.write("(", tc, " ", num, " ", units, ")")
        return
    
    def output_apply(self, expr, paren):
        """Output an <apply> expression.
        
        paren is True if the context has requested parentheses.
        """
        op = expr.operator()
        op_name = op.localName.capitalize()
        self.open_paren(paren)
        # Some operators are special-cased, but most map directly
        if op_name == u'Root':
            self.output_special_apply(expr, u'Root', u'degree', u'Sqrt')
        elif op_name == u'Log':
            self.output_special_apply(expr, u'Log', u'logbase', u'Ln')
        elif op_name == u'Diff':
            if self.single_component:
                # A bit of a hack, due to the odd way this case is
                # handled by other translators - in effect the
                # translator has to do some PE...
                self.write('Apply Diff [Variable ',
                           self.code_name(
                    op.dependent_variable.get_source_variable(recurse=True)),
                           ', Variable ',
                           self.code_name(
                    op.independent_variable.get_source_variable(recurse=True)),
                           ']')
            else:
                self.really_output_apply(op_name, list(expr.operands()) +
                                         [expr.bvar.ci])
        else:
            self.really_output_apply(op_name, expr.operands())
        self.close_paren(paren)
        return

    def really_output_apply(self, operator, operands):
        """Actually output code for the application of an operator."""
        self.write('Apply ', operator, ' [')
        comma = False
        for operand in operands:
            if comma: self.write(',')
            else: comma = True
            self.output_expr(operand, False)
        self.write(']')
        return

    def output_special_apply(self, expr, op_name, qual_name, special_name):
        """Output a special-cased apply expression.

        op_name is the name of the general case operator.  If the
        expression has a qualifier called qual_name, this will be
        used, with the qualifier's value as second operand.  Otherwise,
        the operator special_name will be used, with a single operand.
        """
        if hasattr(expr, qual_name):
            self.really_output_apply(op_name, [expr.operands().next(),
                                               getattr(self, qual_name)])
        else:
            self.really_output_apply(special_name, expr.operands())
        return

    def output_piecewise(self, expr, paren):
        """Output the piecewise expression expr."""
        self.open_paren(paren)
        self.write('Piecewise [')
        comma = False
        for piece in getattr(expr, u'piece', []):
            if comma: self.write(',')
            else: comma = True
            self.write('Case ')
            self.output_expr(child_i(piece, 2), True) # Condition
            self.write(' ')
            self.output_expr(child_i(piece, 1), True) # Result
        self.write('] ')
        if hasattr(expr, u'otherwise'):
            self.write('(Just ')
            self.output_expr(child_i(expr.otherwise, 1), True) # Default case
            self.write(')')
        else:
            self.write('Nothing')
        self.close_paren(paren)
        return

class CellMLToMatlabTranslator(CellMLTranslator):
    """Translate a CellML model to Matlab code.

    The normal case generates a .m file such as could be used with ode45.
    When lookup tables are used (TODO), the file generated represents a
    function that returns a function handle, suitable for use with ODE
    solvers.
    """

    # Language tokens that differ from the default
    COMMENT_START = '% ' # Start of a 1 line comment
    # Types are determined automatically by Matlab
    TYPE_DOUBLE = ''
    TYPE_CONST_DOUBLE = ''
    # Some constants are different
    PI = 'pi'
    E = 'exp(1)'

    def __init__(self, **kwargs):
        super(CellMLToMatlabTranslator, self).__init__(**kwargs)
        # Update some function, etc. names
        self.function_map = {}
        self.function_map.update(CellMLTranslator.function_map)
        self.function_map.update(
            {'power': 'power', 'abs': 'abs',
             'xor': 'xor', 'not': '~',
             'sec': 'sec', 'csc': 'csc', 'cot': 'cot',
             'sech': 'sech', 'csch': 'csch', 'coth': 'coth',
             'arcsec': 'asec', 'arccsc': 'acsc', 'arccot': 'acot',
             'arcsech': 'asech', 'arccsch': 'acsch', 'arccoth': 'acoth'
             })
        self.recip_trig = {}
        self.binary_ops = {}
        self.binary_ops.update(CellMLTranslator.binary_ops)
        del self.binary_ops['xor']
        self.binary_ops['neq'] = '~='
        self.binary_ops['divide'] = './'
        self.nary_ops = {}
        self.nary_ops.update(CellMLTranslator.nary_ops)
        self.nary_ops['times'] = '.*'
        self.special_roots = {2: 'sqrt'}

    def output_file_name(self, model_filename):
        """Generate a name for our output file, based on the input file."""
        name = os.path.splitext(model_filename)[0] + '.m'
        # Matlab doesn't like long names :(
        if len(name) > 60:
            # Take end part so we get version/variant info if present
            name = name[-60:]
        return name

    def translate(self, doc, *args, **kwargs):
        """Generate code for the model or its Jacobian matrix."""
        self.variable_name_map = {}
        if hasattr(doc.model, u'solver_info') and \
               hasattr(doc.model.solver_info, u'jacobian'):
            kwargs['continuation'] = self.output_jacobian
        if 'output_filename' in kwargs and len(kwargs['output_filename'])>60:
            # Take end part so we get version/variant info if present
            kwargs['output_filename'] = kwargs['output_filename'][-60:]
        return super(CellMLToMatlabTranslator,
                     self).translate(doc, *args, **kwargs)

    def output_top_boilerplate(self):
        """Output top boilerplate."""
        self.output_comment(version_comment(self.add_timestamp))
        t = self.code_name(self.free_vars[0])
        # Matlab doesn't like long names :(
        if len(self.class_name) > 60:
            # Take end part so we get version/variant info if present
            self.class_name = self.class_name[-60:]
            # Strip leading underscores
            while self.class_name[0] == '_':
                self.class_name = self.class_name[1:]
        
        if self.use_lookup_tables:
            self.writeln('function dy_fun_ptr = ', self.class_name, '_lt(step)')
            self.output_comment('Generate a function to evaluate using '
                                'lookup tables the model ', self.model.name, '.')
            self.output_comment('The function returned is f, where dU/dt = f(t, U).')
            self.writeln()
            self.set_indent(offset=1)
            self.output_lut_generation()
            self.output_lut_lookups()
            self.writeln('tables = generate_tables(step);')
        else:
            self.writeln('function [dy_fun_ptr initial_values V_index t_units state_var_names] = ',
                         self.class_name, '()')
            self.output_comment('Get evaluation function and metadata for the model ',
                                self.model.name, '.')
            self.output_comment('\nReturns the function f (where dU/dt = f(t, U)),\n'
                                'suitable initial values for the system,\n'
                                'the index of the transmembrane potential within '
                                'the state variable vector,\n'
                                'the multiplicative factor of the time units,\n'
                                'and the names of the state variables.')
            self.set_indent(offset=1)
            self.writeln('V_index = ', self.v_index+1, ';')
            self.writeln('state_var_names = cell(1, ', len(self.state_vars), ');')
            self.writeln('initial_values = zeros(1, ', len(self.state_vars), ');')
            for i, var in enumerate(self.state_vars):
                self.writeln('state_var_names{', i+1, '}', self.EQ_ASSIGN,
                             "'", var.fullname(), "';")
                self.writeln('initial_values(', i+1, ')', self.EQ_ASSIGN,
                             getattr(var, u'initial_value', 'NaN'), ';')
            t_var = self.free_vars[0]
            t_units = t_var.component.get_units_by_name(t_var.units)
            self.writeln('t_units = ', t_units.get_multiplicative_factor(), ';')
        self.writeln('function dy = dy_fun(',t,', y)')
        self.set_indent(offset=1)
        self.output_comment('Time units: ', self.free_vars[0].units)
        self.writeln()
        for i, var in enumerate(self.state_vars):
            self.writeln(self.code_name(var), self.EQ_ASSIGN, 'y(', i+1,
                         ');')
            self.output_comment('Units: ', var.units, '; Initial value: ',
                                getattr(var, u'initial_value', 'Unknown'))
        self.writeln()
        if self.use_lookup_tables:
            for key, i in self.doc.lookup_table_indexes.iteritems():
                i = int(i) + 1
                min, max, step, var = key
                varname = self.code_name(var)
                self.writeln('table_values{', i, '} = lookup_', i,
                             '(tables, ', varname, ', step);')
        self.writeln()
        return

    def output_bottom_boilerplate(self):
        """Output bottom boilerplate."""
        self.writeln()
        self.writeln('dy = zeros(size(y));')
        for i, var in enumerate(self.state_vars):
            self.writeln('dy(', str(i+1), ') = ',
                         self.code_name(var, ode=True), ';')
        self.set_indent(offset=-1)
        self.writeln('end')
        self.writeln()
        self.writeln('dy_fun_ptr = @dy_fun;')
        self.set_indent(offset=-1)
        self.writeln('end')

    def code_name(self, var, ode=False, shorten=True):
        """Matlab has an upper limit on the length of variable names!"""
        full_name = super(CellMLToMatlabTranslator, self).code_name(var, ode)
        if shorten:
            full_name = self.shorten_name(full_name)
        return full_name

    def shorten_name(self, var_name):
        """If the name is too long for Matlab, shorten it."""
        if len(var_name) > 60:
            # Actual bound is 63, but let's be cautious
            try:
                return self.variable_name_map[var_name]
            except KeyError:
                new_name = 'shortened_var_' + str(len(self.variable_name_map))
                self.variable_name_map[var_name] = new_name
                return new_name
        else:
            return var_name

    def output_number(self, expr):
        """Output the plain number expr.
        
        With Matlab, there is no need to make all constants parse as
        doubles to avoid problems with integer division or numbers too
        large for the int type.
        
        Negative numbers will be prefixed by a space to avoid unwanted
        decrement operations.
        """
        n = self.eval_number(expr)
        num = "%.12g" % n
        if num[0] == '-':
            num = ' ' + num
        self.write(num)

    def output_root(self, expr, paren):
        """Output a root taken to some degree.

        If a degree qualifier element is not provided, uses default 2.
        """
        if hasattr(expr, u'degree'):
            # A degree is given.  Compute nthroot(x, b)
            x = expr.operands().next()
            b = expr.degree
            self.output_function('nthroot', [x, b], paren)
        else:
            # Compute square root
            self.output_function('sqrt', expr.operands(), paren)

    def output_piecewise(self, expr, paren):
        """Output the piecewise expression expr.

        Uses an ifexpr.m file to code if expressions.
        """
        num_ifs = 0
        for piece in getattr(expr, u'piece', []):
            num_ifs += 1
            self.write('ifexpr(')
            self.output_expr(child_i(piece, 2), False) # Condition
            self.write(',')
            self.output_expr(child_i(piece, 1), False) # Result
            self.write(',')
        if hasattr(expr, u'otherwise'):
            self.output_expr(child_i(expr.otherwise, 1), paren) # Default case
        else:
            self.write('NaN') # If this is hit, things get ugly
        for i in range(num_ifs):
            self.close_paren(True)

    def output_lut_generation(self):
        """Output code to generate lookup tables.

        There should be a list of suitable expressions available as
        self.doc.lookup_tables, to save having to search the whole
        model.
        """
        self.writeln('function tables = generate_tables(step)')
        self.set_indent(offset=1)
        self.output_comment('Generate all the lookup tables for this model.\n'
                            'Returns a cell array containing matrices, each column of '
                            'which contain one table.')
        self.use_lookup_tables = False
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            min, max, step, var = key
            i = int(idx) + 1
            table_extent = unicode(float(max) - float(min))
            num_tables = unicode(self.doc.lookup_tables_num_per_index[idx])
            self.writeln('tables{', i, '} = zeros(1+floor(', table_extent, '/step),',
                         num_tables, ');')
        for expr in self.doc.lookup_tables:
            j = int(expr.table_name) + 1
            i = int(expr.table_index) + 1
            var = expr.get_component().get_variable_by_name(expr.var)
            varname = self.code_name(var)
            self.writeln(varname, ' = [', expr.min, ':step:', expr.max, '];')
            self.writeln('tables{', i, '}(:,', j, ') = ', nl=False)
            self.output_expr(expr, False)
            self.writeln(';', indent=False)
        self.use_lookup_tables = True
        self.set_indent(offset=-1)
        self.writeln('end')
        self.writeln()

    def output_lut_lookups(self):
        """Output the functions that perform table lookups."""
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            i = int(idx) + 1
            min, max, step, var = key
            self.writeln('function val = lookup_', i, '(tables, var, step)')
            self.set_indent(offset=1)
            self.output_comment('Lookup all tables for variable var')
            self.writeln('if ~isreal(var)')
            self.writeln("error(['Index variable value ' num2str(var) ' is not real'])",
                         indent_offset=1)
            self.writeln('end')
            self.writeln('table_lower = ', min, ';')
            self.writeln('table_upper = ', max, ';')
            self.writeln('if var < table_lower || var >= table_upper')
            self.writeln("error(['Index variable value ' num2str(var) ' outside table bounds'])",
                         indent_offset=1)
            self.writeln('end')
            self.writeln('i = 1 + floor((var - table_lower)/step);')
            self.writeln('y1 = tables{', i, '}(i, :);')
            self.writeln('y2 = tables{', i, '}(i+1, :);')
            self.writeln('var_i = table_lower + step*(i-1);')
            self.writeln('val = y1 + (y2-y1) .* (var-var_i) ./ step;')
            self.set_indent(offset=-1)
            self.writeln('end')
            self.writeln()

    def output_table_lookup(self, expr, paren):
        """Output code to look up expr in the appropriate table."""
        i = int(expr.table_index) + 1
        j = int(expr.table_name) + 1
        self.write('table_values{', i, '}(', j, ')')

    def output_jacobian(self):
        """Generate code to compute the Jacobian matrix for this model."""
        t = self.code_name(self.free_vars[0])
        self.writeln('function J = jacobian(',t,', y)')
        self.set_indent(offset=1)
        self.output_comment('Jacobian matrix for the model ', self.model.name)
        self.output_comment('Evaluates the matrix J, where J(j,i) = d f_i / d u_j')
        # State variable assignments
        for i, var in enumerate(self.state_vars):
            self.writeln(self.code_name(var), self.EQ_ASSIGN, 'y(', str(i+1),
                         ');')
            self.output_comment('Units: ', var.units, '; Initial value: ',
                                getattr(var, u'initial_value', 'Unknown'))
        # Mathematics that the Jacobian depends on
        used_vars = set()
        for entry in self.model.solver_info.jacobian.entry:
            used_vars.update(self._vars_in(entry.math))
        nodeset = self.calculate_extended_dependencies(used_vars)
        self.output_equations(nodeset)
        self.writeln()
        # Jacobian entries
        state_var_names = map(lambda v: self.code_name(v, shorten=False),
                              self.state_vars)
        self.writeln('J = zeros(length(y));')
        for entry in self.model.solver_info.jacobian.entry:
            var_i, var_j = entry.var_i, entry.var_j
            i = state_var_names.index(var_i) + 1
            j = state_var_names.index(var_j) + 1
            self.writeln('J(', j, ',', i, ') = ', nl=False)
            if hasattr(entry.math, u'apply'):
                self.output_expr(entry.math.apply, False)
            else:
                self.output_expr(entry.math.cn, False)
            self.writeln(self.STMT_END, indent=False)
        self.set_indent(offset=-1)
        self.writeln('end')

    def output_variable(self, ci_elt, ode=False):
        """Output a ci element, i.e. a variable lookup."""
        if hasattr(ci_elt, '_cml_variable') and ci_elt._cml_variable:
            self.write(self.code_name(ci_elt.variable, ode=ode))
        else:
            # This ci element is in the solver_info section, thus
            # doesn't have all the extra annotations.  It is a fully
            # qualified name though.
            prefix = ['var_', 'd_dt_'][ode]
            varname = unicode(ci_elt)
            if varname[0] == '(':
                # (compname,varname)
                cname, vname = varname[1:-1].split(u',')
                if self.single_component:
                    varname = vname
                else:
                    varname = cname + '__' + vname
            elif varname == u'delta_t':
                # Special case for the timestep in ComputeJacobian
                prefix = ''
                varname = 'mDt'
            else:
                # var_cname__vname
                varname = varname[4:]
            self.write(self.shorten_name(prefix + varname))
        return


###############################################
# Register translation classes in this module #
###############################################

CellMLTranslator.register(CellMLTranslator, 'C++')
CellMLTranslator.register(CellMLToChasteTranslator, 'Chaste')
CellMLTranslator.register(CellMLToCvodeTranslator, 'CVODE')
CellMLTranslator.register(CellMLToMapleTranslator, 'Maple')
CellMLTranslator.register(CellMLToMatlabTranslator, 'Matlab')
CellMLTranslator.register(CellMLToHaskellTranslator, 'Haskell')


######################################################################
#                         Partial Evaluation                         #
######################################################################

def parteval(doc):
        # Logging
        logger = logging.getLogger('partial-evaluator')
        def debug(*args):
            logger.debug(' '.join(map(str, args)))
        # BTA
        doc.model.do_binding_time_analysis()
        def expr_lhs(expr):
            """Display the LHS of this expression."""
            lhs = expr.assigned_variable()
            if isinstance(lhs, cellml_variable):
                return lhs.fullname()
            else:
                return lhs[0].fullname() + u'/' + lhs[1].fullname()
        # Reduce/evaluate expressions
        while True:
            doc.model._pe_repeat = u'no'
            exprs = [e for e in doc.model.get_assignments()
                     if isinstance(e, mathml_apply)]
            for expr in exprs:
                if expr._get_binding_time() is BINDING_TIMES.static:
                    value = expr.evaluate()
                    debug("Evaluated", expr_lhs(expr),
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
            debug("----- looping -----")
        del doc.model._pe_repeat
        
        # Use canonical variable names on LHS of assignments
        def rename_vars(elt):
            if isinstance(elt, mathml_ci):
                var = elt.variable.get_source_variable(recurse=True)
                elt.xml_remove_child(unicode(elt))
                elt.xml_append(var.fullname(cellml=True))
                elt._cml_variable = var
                debug("Using canonical name", unicode(elt))
            else:
                for e in doc.model.xml_element_children(elt):
                    rename_vars(e)
        for expr in (e for e in doc.model.get_assignments()
                     if isinstance(e, mathml_apply)):
            rename_vars(expr.eq.lhs)

        # Tidy up kept variables, in case they aren't referenced in an eq'n.
        for comp in getattr(doc.model, u'component', []):
            for var in getattr(comp, u'variable', []):
                if var.pe_keep:
                    var._reduce()

        # Collapse into a single component
        new_comp = doc.xml_create_element(u'component', NSS[u'cml'],
                                          attributes={u'name': u'c'})
        old_comps = list(getattr(doc.model, u'component', []))
        doc.model.xml_append(new_comp)
        # We iterate over a copy of the component list so we can
        # delete components from the model in this loop, and so the
        # new component exists in the model so we can add content to
        # it.
        for comp in old_comps:
            # Move relevant contents into new_comp
            for units in getattr(comp, u'units', []):
                # Copy all <units> elements
                # TODO: Just generate the ones we need,
                # using _ensure_units_exist
                units.next_elem = None
                new_comp.xml_append(units)
            for var in getattr(comp, u'variable', []):
                # Only move used source variables
                if (var.get_usage_count() and
                    var.get_type() != VarTypes.Mapped) or var.pe_keep:
                    # Remove from where it was
                    comp._del_variable(var)
                    # Set name to canonical version
                    var.name = var.fullname(cellml=True)
                    # Place in new component
                    new_comp._add_variable(var)
            # Don't copy reactions
            for math in getattr(comp, u'math', []):
                # Copy all <math> elements with content
                if math.xml_children:
                    math.next_elem = None
                    new_comp.xml_append(math)
                    # Invalidate cached component links
                    math._unset_component_links()
            doc.model.xml_remove_child(comp)
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
            rhs = expr.eq.rhs
            expr._cml_depends_on = list(expr.vars_in(rhs))
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



class ConfigurationError(ValueError):
    """Error thrown if configuration file is invalid."""
    pass

class ConfigurationStore(object):
    """
    A container for configuration information, read in from XML
    configuration files.  The file structure is described in the
    read_configuration_file method.
    """
    def __init__(self, doc, options=None):
        """Create a new store.

        doc specifies a CellML document, the processing of which this
        configuration store will affect.

        If given, options should be an optparse.Values instance
        containing command-line options.
        """
        self.doc = doc
        doc._cml_config = self
        self.options = options
        # Transmembrane potential
        self.V_definitions = [u'membrane,V']
        self.V_variable = None
        # Lookup table configuration
        self.lut_config = {}
        self.lut_config_keys = []
        # Ionic currents configuration
        self.i_stim_definitions = [u'membrane,i_Stim']
        self.i_stim_var = None
        self.i_ionic_definitions = [u'membrane,i_.*']
        self.i_ionic_vars = []
        return

    def read_configuration_file(self, config_file):
        """Read configuration stored in config_file.

        The configuration file is expected to be XML, conforming to
        the following structure.  Currently little checking is done on
        the file format; incorrectly formatted files are unlikely to
        give particularly helpful error messages.

        The root element may contain a 'global' element, giving global
        configuration options.  These include:

         * 'lookup_tables'
           Contains one or more 'lookup_table' elements, one for each
           type of lookup table available.  These contain (a selection of)
           the elements:
           * 'var' - the variable to key on.  The component name given
             should be that from which the variable is exported.  Must be
             present.
           * 'min', 'max', 'step' - table bounds parameters.  Optional.
           Default values are used for parameters that are not present.
         * 'currents'
           Defines which variables hold the ionic and stimulus currents,
           if any.  It contains 2 elements:
           * 'stimulus' - the full name of the stimulus current variable
           * 'ionic_match' - a regular expression matching full names of
             ionic current variables.  It may also match the stimulus
             current, but the stimulus will never be considered an ionic
             current.  The value is split on ','; the first part is then
             matched against component names, and the second against
             variables in matching components.
         * 'transmembrane_potential'
           Defines which variable holds the transmembrane potential.
           Defaults to 'membrane,V' if not present.
           
        The root element also contains 0 or more 'for_model' elements,
        which contain settings for individual models.  These must have
        at least one of an 'id' or 'name' attribute, which specify the
        model in question.  They can also contain:
        
         * a 'gating_vars' element, which has a sequence of 'var' elements
           listing gating variables;
         * a 'newton_vars' element, which has a sequence of 'var' elements
           listing variables to be solved for using Newton's method,
           or a similar non-linear solver;
         * and anything allowable as global configuration options.  Options
           given here override those specified globally.

        Configuration which is identical for groups of models may be given
        using the 'for_models' element.  This has an 'ids' element as its
        first child, which contains 'id' elements holding either the name
        or id of a model.  The remaining contents of the 'for_models'
        element are as for 'for_model'.

        There are 3 ways of specifying variables:
        1. By name (var type='name')
           Variable names are given in full form, i.e. component_name,variable_name
        2. By standardised name (var type='oxmeta')
           Use the name from the oxmeta annotations
        3. By reference to a section of this config file (when defining lookup table keys)
           e.g. <var type='config-name'>transmembrane_potential</var>

        Within any element that specifies a variable, a list of <var> elements can be
        provided.  Each will be tried in turn to see if a match can be found in the model,
        and the first match wins.
        Alternatively these elements may have text content which must be of the form
        'component_name,variable_name'.

        Some items are overridden if oxmeta annotations are present and the
        use_metadata option is set.  Currently this just applies to the
        stimulus current and transmembrane potential.
        """
        rules = [bt.ws_strip_element_rule(u'*')]
        config_doc = amara_parse(config_file, rules=rules)
        # Parse global configuration
        glo = config_doc.xml_xpath(u'/*/global')
        if glo:
            glo = glo[0]
            if hasattr(glo, u'lookup_tables'):
                # Lookup tables configuration
                self._parse_lookup_tables(glo.lookup_tables)
            if hasattr(glo, u'currents'):
                # Configure which vars are ionic and stimulus currents
                self._parse_currents(glo.currents)
            if hasattr(glo, u'transmembrane_potential'):
                # Configure the transmembrane potential variable
                self._parse_Vm(glo.transmembrane_potential)
        # Get the config section(s) for our model.  Sections
        # specifically for this model come after sections covering
        # multiple models, so they take precedence.
        model_id = getattr(self.doc.model, u'id', self.doc.model.name)
        sections = config_doc.xml_xpath(
            u'/*/for_models[ids/id="%s" or ids/id="%s"]'
            % (self.doc.model.name, model_id))
        sections.extend(config_doc.xml_xpath(
            u'/*/for_model[@name="%s" or @id="%s"]'
            % (self.doc.model.name, model_id)))
        # Newton variables
        newton_vars = []
        for section in sections:
            if hasattr(section, 'newton_vars'):
                for varname in section.newton_vars.var:
                    newton_vars.append(unicode(varname).strip())
        if newton_vars:
            #print "Non-linear vars:", newton_vars
            self.doc.model._cml_nonlinear_system_variables = u';'.join(newton_vars)
        # Overrides for global configuration
        for section in sections:
            if hasattr(section, u'lookup_tables'):
                self._parse_lookup_tables(section.lookup_tables)
            if hasattr(section, u'currents'):
                self._parse_currents(section.currents)
            if hasattr(section, u'transmembrane_potential'):
                self._parse_Vm(section.transmembrane_potential)
        # Now identify the variables in the model
        self.find_current_vars()
        self.find_transmembrane_potential()
        return

    def _create_var_def(self, content, defn_type):
        """Create a variable definition object."""
        xml_fragment = '<var type="%s">%s</var>' % (defn_type, content)
        return amara.parse(str(xml_fragment)).var

    def _check_var_def(self, var_elt, var_desc):
        """Check a variable definition is syntactically valid.
        
        If type == 'name', it must have text content of the form "component_name,variable_name".
        If type == 'oxmeta', it must have text content drawn from METADATA_NAMES.
        If type == 'config-name', it must have text content either 'stimulus' or 'transmembrane_potential'.
        """
        defn_type = getattr(var_elt, u'type', u'name')
        if defn_type == u'name':
            name_parts = unicode(var_elt).strip().split(',')
            if len(name_parts) != 2:
                raise ConfigurationError('Invalid definition of ' + var_desc + ': '
                                         + unicode(var_elt))
        elif defn_type == u'oxmeta':
            if unicode(var_elt) not in METADATA_NAMES:
                raise ConfigurationError('"' + unicode(var_elt) + '" is not a valid oxmeta name')
        elif defn_type == u'config-name':
            if unicode(var_elt) not in [u'stimulus', u'transmembrane_potential']:
                raise ConfigurationError('"' + unicode(var_elt) + '" is not a name known to the config file')
        else:
            raise ConfigurationError('"' + defn_type + '" is not a valid variable definition type')
        return

    def _parse_Vm(self, vm_elt):
        """Parse definition of variable holding the transmembrane potential."""
        if hasattr(vm_elt, 'var'):
            # List of possibilities
            self.V_definitions = []
            for vardef in vm_elt.var:
                self._check_var_def(vardef, 'transmembrane potential')
                self.V_definitions.append(vardef)
        else:
            # Old style - single variable given by text content
            self._check_var_def(vm_elt, 'transmembrane potential')
            self.V_definitions = [vm_elt]
        return

    def _parse_currents(self, currents):
        """Parse definitions of ionic and stimulus currents."""
        if hasattr(currents, u'stimulus'):
            if hasattr(currents.stimulus, u'var'):
                # List of possibilities
                self.i_stim_definitions = []
                for var in currents.stimulus.var:
                    self._check_var_def(var, 'stimulus current')
                    self.i_stim_definitions.append(var)
            else:
                self._check_var_def(currents.stimulus, 'stimulus current')
                self.i_stim_definitions = [currents.stimulus]
        if hasattr(currents, u'ionic_match'):
            if hasattr(currents.ionic_match, u'var'):
                # List of possibilities
                self.i_ionic_definitions = []
                for var in currents.ionic_match.var:
                    self._check_var_def(var, 'ionic currents')
                    self.i_ionic_definitions.append(var)
            else:
                self._check_var_def(currents.ionic_match, 'ionic currents')
                self.i_ionic_definitions = [currents.ionic_match]
        return
    
    def _find_variable(self, defn, pe_done=False):
        """Find a variable matching the given definition.

        If pe_done is True, then partial evaluation has been performed
        on the model, so looking for variables by name needs to look for
        variables called compname__varname in the single component.
        """
        defn_type = getattr(defn, u'type', u'name')
        if defn_type == u'name':
            name_parts = unicode(defn).strip().split(',')
            if pe_done:
                try:
                    var = self.doc.model.component.get_variable_by_name(u'__'.join(name_parts))
                except KeyError:
                    var = None
            else:
                var = self.doc.model.xml_xpath(u'cml:component[@name="%s"]/cml:variable[@name="%s"]'
                                               % tuple(name_parts))
                if var:
                    var = var[0]
        elif defn_type == u'oxmeta':
            var = self.doc.model.xml_xpath(u'cml:component/cml:variable[@oxmeta:name="%s"]'
                                           % unicode(defn))
            if var:
                var = var[0]
        elif defn_type == u'config-name':
            if unicode(defn) == u'stimulus':
                var = self.i_stim_var
            elif unicode(defn) == u'transmembrane_potential':
                var = self.V_variable
        else:
            raise ConfigurationError('"' + defn_type + '" is not a valid variable definition type')
        return var

    def find_current_vars(self):
        """Find the variables representing currents."""
        # Try metadata first if specified on command line
        if self.options.use_metadata:
            i_stim = self.doc.model.xml_xpath(u'cml:component/cml:variable[@oxmeta:name="membrane_stimulus_current"]')
            if i_stim:
                self.i_stim_var = i_stim[0]
            else:
                DEBUG('metadata', 'membrane_stimulus_current not found in metadata enabled cellml file.')
            # TODO: add ionic currents here? (when using RDF metadata)
        if not self.i_stim_var:
            for defn in self.i_stim_definitions:
                var = self._find_variable(defn)
                if var:
                    self.i_stim_var = var
                    break
            else:
                # No match :(
                print "No stimulus current found; you'll have trouble generating Chaste code"
        # Other ionic currents just set from config file
        self.i_ionic_vars = []
        for defn in self.i_ionic_definitions:
            if getattr(defn, u'type', u'name') != u'name':
                raise ConfigurationError('Ionic current definitions have to have type "name"')
            regexps = unicode(defn).strip().split(',')
            comp_re = re.compile(regexps[0] + '$')
            var_re = re.compile(regexps[1] + '$')
            for component in getattr(self.doc.model, u'component', []):
                if comp_re.match(unicode(component.name).strip()):
                    for var in getattr(component, u'variable', []):
                        if (var is not self.i_stim_var and
                            var_re.match(unicode(var.name).strip())):
                            self.i_ionic_vars.append(var)
        if not self.i_ionic_vars:
            print "No ionic currents found; you'll have trouble generating Chaste code"
        return

    def _parse_lookup_tables(self, lookup_tables):
        """Parse a lookup_tables element."""
        for lt in lookup_tables.lookup_table:
            var_type = getattr(lt.var, u'type', u'name')
            var_name = unicode(lt.var).strip()
            config_key = (var_type, var_name)
            if not self.lut_config.has_key(config_key):
                self.lut_config[config_key] = {}
                self._set_lut_defaults(self.lut_config[config_key])
                self.lut_config_keys.append(config_key)
            for elt in [e for e in lt.xml_children
                        if isinstance(e, amara.bindery.element_base)
                        and e.localName != u'var']:
                self.lut_config[config_key]['table_' + elt.localName] = unicode(elt).strip()
        return

    def _set_lut_defaults(self, lut_dict):
        """Set default configuration for a lookup table."""
        def_dict = LookupTableAnalyser._LT_DEFAULTS
        for k, v in def_dict.iteritems():
            if k != 'table_var':
                lut_dict[k] = v
        return

    def annotate_currents_for_pe(self):
        "Annotate ionic & stimulus current vars so PE doesn't remove them."
        if self.i_stim_var:
            self.i_stim_var.set_pe_keep(True)
        for var in self.i_ionic_vars:
            var.set_pe_keep(True)
        return
    
    def annotate_metadata_for_pe(self):
        "Annotate all vars tagged with metadata so PE doesn't remove them."
        for var in self.metadata_vars:
            var.set_pe_keep(True)        
        return

    def find_transmembrane_potential(self):
        """Find and store the variable object representing V.

        Uses the name given in the command line options, if present.
        Otherwise uses first metadata, if present, then the configuration file.
        """
        if not self.options:
            raise ValueError('No command line options given')
        # Check command line option
        if self.options.transmembrane_potential:
            self.V_definitions = [self.options.transmembrane_potential.strip().split(',')]
            if len(self.V_definitions[0]) != 2:
                raise ConfigurationError('The name of V must contain both '
                                         'component and variable name')
        # Check for metadata annotation
        elif self.options.use_metadata:
            var = self.doc.model.xml_xpath(u'cml:component/cml:variable[@oxmeta:name="membrane_voltage"]')
            if var:
                self.V_variable = var[0]
        if not self.V_variable:
            for defn in self.V_definitions:
                var = self._find_variable(defn)
                if var:
                    self.V_variable = var
                    break
        if not self.V_variable:
            raise ConfigurationError('No transmembrane potential found; check your configuration')
        return self.V_variable

    def find_lookup_variables(self, pe_done=False):
        """Find the variable objects used as lookup table keys.

        This method translates the variable names given in the
        configuration file into objects in the document, and then uses
        those objects as keys in our lut_config dictionary.

        If pe_done is True, then partial evaluation has been performed
        on the model, so look for variables called compname__varname
        in the single component.  Otherwise, we look for variables
        called varname in component compname.
        """
        new_config = {}
        for key in self.lut_config_keys:
            defn_type, content = key
            defn = self._create_var_def(content, defn_type)
            var = self._find_variable(defn, pe_done)
            if not var:
                # Variable doesn't exist, so we can't index on it
                LOG('lookup-tables', logging.WARNING, 'Variable', content,
                    'not found, so not using as table index.')
            else:
                if not var in new_config:
                    new_config[var] = {}
                new_config[var].update(self.lut_config[key])
        self.lut_config = new_config
        return

    # TODO - move into seperate metadata class?
    def validate_metadata(self):
        """Check that the metadata annotations are 'sensible'.
        
        Ensures that only names we know are used, and that the same name isn't used for multiple
        variables.
        """
        self.metadata_vars = self.doc.xml_xpath(u'/*/*/cml:variable[@oxmeta:name]')
        names_used = [var.oxmeta_name for var in self.metadata_vars]
        # Check all metadata is allowed
        if frozenset(names_used) <= METADATA_NAMES:
            DEBUG('metadata', 'Metadata values are valid')
        else:
            DEBUG('metadata', 'Metadata values are NOT valid')
            raise ConfigurationError('Metadata values are NOT valid, try running with --assume-valid')
        # Check for duplicates
        d = {}
        for name in names_used:
            if name in d:
                raise ConfigurationError(name + ' metadata attribute is duplicated in the cellml file.')
            else:
                d[name] = name

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


######################################################################
#                    For running as an executable                    #
######################################################################

def get_options(args):
    """get_options(args):
    Process our command-line options.

    args is a list of options & positional arguments.
    """
    usage = 'usage: %prog [options] <cellml file or URI>'
    parser = optparse.OptionParser(version="%%prog %s" % __version__,
                                   usage=usage)
    parser.add_option('-o', dest='outfilename', metavar='OUTFILE',
                      help="write program code to OUTFILE [default action is"
                      " to use the input filename with a different extension]")
    parser.add_option('-l', '--lookup-tables',
                      dest='lut', action='store_true', default=False,
                      help="perform a lookup table analysis")
    parser.add_option('-p', '--pe', '--partial-evaluation',
                      dest='pe', action='store_true', default=False,
                      help="partially evaluate the model")
    parser.add_option('-u', '--units-conversions',
                      action='store_true', default=False,
                      help="add explicit units conversion mathematics")
    parser.add_option('-c', '--class-name',
                      help="explicitly set the name of the generated class")
    parser.add_option('-a', '--augment-class-name',
                      dest='augment_class_name', action='store_true',
                      default=False,
                      help="alter the class name to show what transformations"
                      " are used")
    parser.add_option('-C', '--output-cellml',
                      dest='translate', action='store_false',
                      help="output an annotated CellML file, on stdout unless"
                      " -o specified")
    parser.add_option('-T', '--translate',
                      dest='translate', action='store_true',
                      default=True,
                      help="output computer code [default]")
    translators = sorted(CellMLTranslator.translators)
    parser.add_option('-t', '--translate-type',
                      type='choice', choices=translators,
                      default='Chaste', metavar='TYPE',
                      help="the type of code to output [default: %default].  "
                      "Choices: " + str(translators))
    parser.add_option('-d', '--debug', action='store_true', default=False,
                      help="output debug info to stderr")
    parser.add_option('--Wu', '--warn-on-units-errors',
                      action='store_true', default=False,
                      dest='warn_on_units_errors',
                      help="give a warning instead of an error for"
                      " dimensional inconsistencies")
    parser.add_option('--no-timestamp',
                      action='store_true', default=False,
                      help="don't add a timestamp comment to generated files")
    parser.add_option('--config-file',
                      action='store',
                      help="pathname of configuration file")
    parser.add_option('--omit-constants',
                      action='store_true', default=False,
                      help="when generating Maple code, don't include "
                      "assignments of constants")
    parser.add_option('--compute-full-jacobian',
                      action='store_true', default=False,
                      help="make generated Maple code compute the full Jacobian"
                      " matrix, rather than just that for the nonlinear portion"
                      " of the ODE system")
    parser.add_option('-J', '--do-jacobian-analysis',
                      action='store_true', default=False,
                      help="experimental Jacobian analysis; implies -t Maple")
    parser.add_option('-V', '--transmembrane-potential',
                      default=None, metavar='POT_VAR',
                      help=
                      "POT_VAR is the full name of the variable representing"
                      " the transmembrane potential.  If not specified here,"
                      " the configuration file will be used.  Defaults to "
                      "'membrane,V'.")
    parser.add_option('-j', '--maple-output',
                      metavar='FILENAME', default=None,
                      help="file containing output from a Maple script "
                      "generated using -J.  The generated code/CellML will "
                      "then contain a symbolic Jacobian as computed by Maple.")
    parser.add_option('--use-chaste-stimulus',
                      action='store_true', default=False,
                      help="when generating Chaste code, use Chaste's stimulus"
                      " rather than that defined in the model")
    parser.add_option('--no-separate-lut-class', dest='separate_lut_class',
                      action='store_false', default=True,
                      help="don't put lookup tables in a separate class")
    parser.add_option('--row-lookup-method',
                      action='store_true', default=False,
                      help="add a method to look up a whole row of a table")
    parser.add_option('--lt-index-uses-floor',
                      action='store_true', default=False,
                      help="use floor() to calculate LT indices, instead of "
                      "just casting")
    parser.add_option('--constrain-table-indices',
                      action='store_true', default=False,
                      help="constraint lookup table index variables to remain"
                      " within the bounds specified, rather than throwing an"
                      " exception if they go outside the bounds")
    parser.add_option('-i', '--convert-interfaces',
                      action='store_true', default=False,
                      help="perform units conversions at interfaces to Chaste. "
                      "(only works if -t Chaste is used)")
    parser.add_option('--assume-valid',
                      action='store_true', default=False,
                      help="skip some of the model validation checks")
    parser.add_option('--no-member-vars', dest='kept_vars_as_members',
                      action='store_false', default=True,
                      help="[debug] don't store kept variables as members")
    parser.add_option('--bad-lt-layout-for-cache', dest='bad_tables_for_cache',
                      action='store_true', default=False,
                      help="[debug] use the old LT layout, with poorer cache"
                      " performance")
    parser.add_option('-m', '--use-metadata',
                      dest='use_metadata', action='store_true', default=False,
                      help="obtain variable names from metadata embedded within "
                      "the CellML file, rather than config.xml")
    parser.add_option('-y', '--dll', '--dynamically-loadable',
                      dest='dynamically_loadable',
                      action='store_true', default=False,
                      help="add code to allow the model to be compiled to a "
                      "shared library and dynamically loaded")

    options, args = parser.parse_args(args)
    if len(args) != 1:
        parser.error("exactly one input CellML file must be specified")
    return options, args[0]


class NotifyHandler(logging.Handler):
    """
    A logging handler that just notes if any messages are logged.
    """
    def __init__(self, level=logging.NOTSET):
        logging.Handler.__init__(self, level=level)
        self.reset()

    def emit(self, record):
        self.messages = True

    def reset(self):
        """Reset the handler, as if no messages have occurred."""
        self.messages = False

def run():
    # Translate the file given on the command line
    options, model_file = get_options(sys.argv[1:])

    if options.debug:
        formatter = logging.Formatter(fmt="%(name)s: %(message)s")
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(formatter)
        handler.addFilter(OnlyDebugFilter())
        logging.getLogger().addHandler(handler)
        logging.getLogger().setLevel(logging.DEBUG)        

    # We can't translate if some warnings occur, as well as if the
    # model is invalid
    notifier = NotifyHandler(level=logging.WARNING_TRANSLATE_ERROR)
    logging.getLogger('validator').addHandler(notifier)
    v = validator.CellMLValidator()
    valid, doc = v.validate(model_file, True,
                            warn_on_units_errors=options.warn_on_units_errors,
                            assume_valid=options.assume_valid)
    v.quit()
    del v

    if not valid or notifier.messages:
        print >>sys.stderr, model_file,
        if not valid:
            print >>sys.stderr, "is not a valid CellML file"
        else:
            print >>sys.stderr, "contains untranslatable constructs"
        sys.exit(1)

    config = ConfigurationStore(doc, options=options)
    if options.config_file:
        config.read_configuration_file(options.config_file)
    else:
        # Use defaults
        config.find_current_vars()
        config.find_transmembrane_potential()

    if options.use_metadata:
        DEBUG('metadata', "metadata enabled, using for ConfigStore")
        if not options.assume_valid:
            config.validate_metadata()

    class_name = getattr(options, 'class_name', None)
    if options.augment_class_name and not class_name:
        class_name = u'CML_' + doc.model.name.replace('-', '_')
        if options.pe:
            class_name += '_pe'
        if options.lut:
            class_name += '_lut'
        if options.maple_output:
            class_name += '_be'
        if options.use_metadata:
            class_name += '_sens'

    output_filename = getattr(options, 'outfilename', None)
    if not options.translate and not output_filename:
        output_filename = 'stdout'

    if options.units_conversions:
        doc.model.add_units_conversions()

    if options.do_jacobian_analysis:
        lin = LinearityAnalyser()
        lin.analyse_for_jacobian(doc, V=config.V_variable)
        options.translate_type = 'Maple'

    if options.pe:
        # Do partial evaluation
        if options.config_file:
            config.annotate_currents_for_pe()
        # Need to ensure pe doesn't remove metadata-annotated variables
        if options.use_metadata:
            config.annotate_metadata_for_pe()
        parteval(doc)

    if options.lut:
        # Do the lookup table analysis
        lut = LookupTableAnalyser()
        if options.config_file:
            config.find_lookup_variables(options.pe)
        elif options.pe:
            lut.set_params(table_var=u'membrane__V')
        lut.analyse_model(doc)

    if options.maple_output:
        # Parse Jacobian matrix
        from maple_parser import MapleParser
        mp = MapleParser()
        jacobian_file = file(options.maple_output) # TODO: Error checking
        doc.model._cml_jacobian = mp.parse(jacobian_file)
        jacobian_file.close()
        # Rearrange linear ODEs
        lin = LinearityAnalyser()
        lin.analyse_for_jacobian(doc, V=config.V_variable)
        lin.rearrange_linear_odes(doc)
        # Add info as XML
        doc.model._add_solver_info()
        # TODO: Analyse the XML, adding cellml_variable references, etc.

    if options.translate:
        # Translate to code
        klasses = CellMLTranslator.translators
        klass = klasses[options.translate_type]
        initargs = {'add_timestamp': not options.no_timestamp}
        transargs = {'v_variable': config.V_variable}
        transargs['row_lookup_method'] = options.row_lookup_method
        transargs['lt_index_uses_floor'] = options.lt_index_uses_floor
        transargs['constrain_table_indices'] = options.constrain_table_indices
        transargs['bad_tables_for_cache'] = options.bad_tables_for_cache
        if issubclass(klass, CellMLToMapleTranslator):
            initargs['omit_constants'] = options.omit_constants
            initargs['compute_full_jacobian'] = options.compute_full_jacobian
        elif issubclass(klass, CellMLToChasteTranslator):
            doc.model._add_solver_info_ionic_current()
            transargs['use_chaste_stimulus'] = options.use_chaste_stimulus
            transargs['separate_lut_class'] = options.separate_lut_class
            transargs['convert_interfaces'] = options.convert_interfaces
            transargs['kept_vars_as_members'] = options.kept_vars_as_members
            transargs['use_metadata'] = options.use_metadata
            transargs['dynamically_loadable'] = options.dynamically_loadable
        t = klass(**initargs)
        t.translate(doc, model_file, output_filename, class_name=class_name, **transargs)
    else:
        # Add a comment element
        comment = pycml.comment_base(
            body=u'\n' + version_comment(not options.no_timestamp) + u'\n')
        doc.xml_insert_before(doc.model, comment)
        # Output annotated model
        stream = open_output_stream(output_filename)
        doc.xml(indent=u'yes', stream=stream)
        close_output_stream(stream)


if __name__ == '__main__':
    run()

    # For use in testing
    def euler(nsteps=1000, dt=0.01):
        global tvar, state_vars, exprs
        tvar = t.free_vars[0]
        state_vars = t.state_vars
        for var in state_vars:
            var.set_value(float(var.initial_value))
        tvar.set_value(0.0)
        exprs = [e for e in doc.model.get_assignments()
                 if isinstance(e, mathml_apply)]
        for i in range(nsteps):
            for expr in exprs:
                expr.evaluate()
            tvar.set_value(tvar.get_value() + dt)
            for var in state_vars:
                var.set_value(var.get_value() +
                              dt * var.get_value(ode=tvar))
        return

    def writefile(outfn='test.cml'):
        # Write out CellML file
        st = open_output_stream(outfn)
        doc.xml(indent=1, stream=st)
        st.close()
        return

    def show_usage():
        for comp in doc.model.component:
            for var in comp.variable:
                print var.fullname(), var._cml_usage_count


    def fix_divide_by_zero():
        """
        Several models have equations of a form that may give rise to
        a divide by zero error on simulation, especially when lookup
        tables are used.  The general form is:
        
        (a * (V - v0)) / (exp(b * (V - v0)) - 1)
        
        When V = v0 this is undefined, however the limit of the
        function as V approaches v0 from either side is well-defined,
        and each limit is the same.  We approximate the limit by
        linear interpolation between values of the expression for
        (V-v0) = +/- 1e-10.
        """
        divides = [d.xml_parent
                   for d in doc.xml_xpath(u'//m:apply/m:divide')]
        for divide in divides:
            pass
        return
