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
This part of PyCml deals with converting CellML models into
programming language code, primarily C++ compatible with Chaste, but
supporting a few other languages also (and easily extensible).
"""

import optparse
import os
import re
import time
import sys

# Make sure PyCml is on sys.path
pycml_path = os.path.dirname(os.path.realpath(__file__))
sys.path[0:0] = [pycml_path]

# Common CellML processing stuff
import pycml
from pycml import *  # Put contents in the local namespace as well
import validator
import optimize


__version__ = "$Revision$"[11:-2]

def version_comment(note_time=True):
    """Generate a version comment, with optional time info."""
    if note_time:
        t = '\non ' + time.asctime()
    else:
        t = ''
    text = """Processed by pycml - CellML Tools in Python
    (translate: %s, pycml: %s, optimize: %s)%s""" % (
        __version__, pycml.__version__, optimize.__version__, t)
    return text


def debugexpr(e):
    "For debugging."
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
    """Error thrown if CellML translation fails."""
    pass

class ConfigurationError(ValueError):
    """Error thrown if configuration file is invalid."""
    pass

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
    TYPE_VOID = 'void '
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

    def __init__(self, add_timestamp=True, options=None):
        """Create a translator."""
        self.options = options
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
                        "Free vars:" + str(self.free_vars)])
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
        start = kwargs.get('start', self.COMMENT_START)
        if kwargs.get('pad', False):
            start = ' ' + start
        comment = ''.join(map(str, args))
        lines = comment.split('\n')
        for line in lines:
            self.writeln(start, line, **kwargs)

    def output_doxygen(self, *args, **kwargs):
        """Output a (multi-line) string as a Doxygen comment."""
        kwargs['start'] = self.DOXYGEN_COMMENT_START
        self.output_comment(*args, **kwargs)

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
        if var.component.name == u'':
            # Special case variables, that don't really exist in the model
            name = var.name
        elif self.single_component:
            name = prefix + var.name
        else:
            name = prefix + var.component.name + '__' + var.name
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
                self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(expr),
                             self.EQ_ASSIGN,
                             self.code_name(expr.get_source_variable()),
                             self.STMT_END, nl=False)
                self.output_comment(expr.units, indent=False, pad=True)
            elif t == VarTypes.Constant:
                self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(expr),
                             self.EQ_ASSIGN, nl=False)
                self.output_number(expr.initial_value)
                self.writeln(self.STMT_END, indent=False, nl=False)
                self.output_comment(expr.units, indent=False, pad=True)
        else:
            # This is a mathematical expression
            self.writeln(self.TYPE_CONST_DOUBLE, nl=False)
            opers = expr.operands()
            self.output_lhs(opers.next())
            self.write(self.EQ_ASSIGN)
            self.output_expr(opers.next(), False)
            self.writeln(self.STMT_END, indent=False, nl=False)
            #1365: add a comment with the LHS units
            self.output_comment(expr._get_element_units(expr.eq.lhs, return_set=False).description(),
                                indent=False, pad=True)

    def output_lhs(self, expr):
        """Output the left hand side of an assignment expression."""
        if expr.localName == 'ci':
            self.output_variable(expr)
        elif expr.operator().localName == 'diff':
            self.write(self.code_name(expr.operator().dependent_variable,
                                      ode=True))

    def output_variable(self, ci_elt, ode=False):
        """Output a ci element, i.e. a variable lookup."""
        self.write(self.code_name(ci_elt.variable, ode=ode))

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
    def calculate_extended_dependencies(self, nodes, prune=[], prune_deps=[]):
        """Method moved to cellml_model."""
        return self.model.calculate_extended_dependencies(nodes, prune, prune_deps)

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
            if self.single_component and cname != self.model.component.name:
                name = cname + u'__' + vname
            else:
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

        Will include state variables.  Also if an expression is being
        replaced by a lookup table, will only return the table key variable.
        """
        res = set()
        if self.use_lookup_tables and isinstance(expr, mathml) and self.is_lookup_table(expr):
            key_var = self.varobj(expr.getAttributeNS(NSS['lut'], u'var'))
            key_var = key_var.get_source_variable(recurse=True)
            res.add(key_var)
        elif isinstance(expr, mathml_ci):
            varname = unicode(expr)
            varobj = self.varobj(varname.strip())
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
    # same .cpp file though).  This can then be made a singleton class
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
                self.writeln('unsigned _table_index_', idx, self.STMT_END)
                self.writeln('double _factor_', idx, self.STMT_END)
                if self.row_lookup_method:
                    self.writeln('double* _lt_', idx, '_row', self.STMT_END)
        self.writeln()
        return

    def output_lut_indices(self):
        """Output declarations for the lookup table indices."""
        self.output_comment('Lookup table indices')
        for key, idx in self.doc.lookup_table_indexes.iteritems():
            self.writeln('unsigned _table_index_', idx, self.STMT_END)
            self.writeln('double _factor_', idx, self.STMT_END)
            if self.row_lookup_method:
                self.writeln('double* _lt_', idx, '_row', self.STMT_END)
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
    
    def contained_table_indices(self, node):
        """Return all lookup tables used directly in computing this node.

        If this is an expression node, checks all its children for table
        lookups, and returns the set of table indices used.
        """
        result = set()
        if isinstance(node, amara.bindery.element_base):
            if self.is_lookup_table(node):
                result.add(node.table_index)
            else:
                for child in node.xml_children:
                    result.update(self.contained_table_indices(child))
        return result

    def output_table_lookup(self, expr, paren):
        """Output code to look up expr in the appropriate table."""
        i = expr.table_index
        if self.row_lookup_method:
            self.write('_lt_', i, '_row[', expr.table_name, ']')
        else:
            self.write(self.lookup_method_prefix, '_lookup_', expr.table_name,
                       '(_table_index_', i, ', _factor_', i, ')')
        return

    def output_table_index_generation(self, indexes_as_member=False, nodeset=set()):
        """Output code to calculate indexes into any lookup tables.
        
        If indexes_as_member is True then the index variables should be
        considered to be member variables, rather than locally declared.
        
        If nodeset is given, then filter the table indices calculated so
        that only those needed to compute the nodes in nodeset are defined.
        """
        self.output_comment('Lookup table indexing')
        if indexes_as_member:
            index_type = ''
            factor_type = ''
            row_type = ''
        else:
            index_type = 'unsigned '
            factor_type = 'double '
            row_type = 'double* '
        tables_to_index = set()
        for node in nodeset:
            tables_to_index.update(self.contained_table_indices(node))
        for key, i in self.doc.lookup_table_indexes.iteritems():
            if nodeset and i not in tables_to_index:
                continue
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
                self.writeln(row_type, '_lt_', i, '_row = ',
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

    # Type of (a reference to) the state variable vector
    TYPE_VECTOR = 'std::vector<double> '
    TYPE_VECTOR_REF = 'std::vector<double>& '
    
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
                      'use_modifiers': False,
                      'dynamically_loadable': False,
                      'use_protocol': False
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
    
    def determine_units_conversion(self, varobj, to_chaste=True):
        """Determine if we needs to units-convert the given quantity.
        
        Will check if varobj has units dimensionally equivalent to those of time
        or the stimulus current.  If it does, then a conversion may be required.
        Returns an appropriate conversion expression (just the code_name of the
        variable if no conversion is required) paired with a list of nodes used
        in computing this conversion.
        
        By default converts from model units to Chaste units.  If to_chaste is
        False, performs the reverse conversion.
        """
        var_units = self.get_var_units(varobj)
        conversion = self.code_name(varobj)
        nodes_used = []
        # Units to check against
        ms = cellml_units.create_new(
            self.model, 'milliseconds',
            [{'units': 'second', 'prefix': 'milli'}])
        current_units = self.get_var_units(self.doc._cml_config.i_stim_var)
        if var_units.dimensionally_equivalent(ms): # Check against time
            if not var_units.equals(ms):
                if not to_chaste:
                    conversion_factor = (ms.get_multiplicative_factor() /
                                         var_units.get_multiplicative_factor())
                else:
                    conversion_factor = (var_units.get_multiplicative_factor() /
                                         ms.get_multiplicative_factor())
                conversion += ' * ' + str(conversion_factor)
        elif var_units.dimensionally_equivalent(current_units): # Check against current
            conversion, nodes_used = self.ionic_current_units_conversion(conversion, var_units, to_chaste)
        return conversion, nodes_used
    
    def ionic_current_units_conversion(self, varname, all_units, to_chaste=True):
        """
        Check whether the units of the transmembrane currents are as expected by
        Chaste, and return an appropriate conversion expression if not, along with
        a list of nodes used in computing this conversion.
        
        By default converts from model units to Chaste units.  If to_chaste is
        False, performs the reverse conversion.
        
        There are three main cases:
        
        1. Current in amps/area, capacitance in farads/area.
           In this case the dimensions are as expected by Chaste, so we just
           compute the units conversion factor.
        2. Current in amps/farads.
           In this case we convert to uA/uF then multiply by Chaste's value
           for the membrane capacitance (in uF/cm^2).
        3. Current in amps, capacitance in farads.
           We assume the cell model conceptually represents a cell, and hence
           that its membrane capacitance is supposed to represent the same
           thing as Chaste's.  Thus convert current to uA, capacitance to uF,
           and return current/capacitance * Chaste's capacitance.
        
        If the model doesn't match one of these cases, print a warning message,
        or die if running fully automatically.
        """
        # We just use the first units element; if the model passes units checking,
        # this is OK.  If it doesn't, then it's pot luck as to whether this
        # method does anything sensible at all...
        model_units = all_units[0]
        conversion = varname # Default to no conversion
        nodes_used = []
        chaste_units = cellml_units.create_new(
            self.model, 'uA_per_cm2',
            [{'units': 'ampere', 'prefix': 'micro'},
             {'units': 'metre', 'prefix': 'centi', 'exponent': '-2'}])
        microamps = cellml_units.create_new(self.model, u'microamps',
                                            [{'units':'ampere', 'prefix':'micro'}])
        microfarads = cellml_units.create_new(self.model, u'microfarads',
                                              [{'units':'farad', 'prefix':'micro'}])
        A_per_F = cellml_units.create_new(self.model, 'A_per_F',
                                          [{'units': 'ampere'},
                                           {'units': 'farad', 'exponent': '-1'}])
        
        if to_chaste:
            times, divide = ' * ', ' / '
        else:
            times, divide = ' / ', ' * '
        
        def problem(msg):
            if self.options.fully_automatic:
                raise TranslationError(msg)
            else:
                print >>sys.stderr, msg
        
        if chaste_units.dimensionally_equivalent(model_units):
            # Case 1
            conv = (model_units.get_multiplicative_factor() /
                    chaste_units.get_multiplicative_factor())
            if conv != 1:
                conversion += times + str(conv)
        elif A_per_F.dimensionally_equivalent(model_units):
            # Case 2
            conv = (model_units.get_multiplicative_factor() /
                    A_per_F.get_multiplicative_factor())
            if conv != 1:
                conversion += times + str(conv)
            conversion += times + 'HeartConfig::Instance()->GetCapacitance()'
        elif microamps.dimensionally_equivalent(model_units):
            # Case 3, probably
            conv = (model_units.get_multiplicative_factor() /
                    microamps.get_multiplicative_factor())
            if conv != 1:
                conversion += times + str(conv)
            model_Cm = self.doc._cml_config.Cm_variable
            if not model_Cm:
                problem("Cannot convert ionic current from amps to uA/cm^2 "
                        "without knowing which variable in the cell model "
                        "represents the membrane capacitance.")
            else:
                nodes_used.append(model_Cm)
                Cm_units = self.get_var_units(model_Cm)
                Cm_value = self.code_name(model_Cm)
                conv = (Cm_units.get_multiplicative_factor() /
                        microfarads.get_multiplicative_factor())
                if conv != 1:
                    conversion += divide + '(' + str(conv) + ' * ' + Cm_value + ')'
                else:
                    conversion += divide + Cm_value
            conversion += times + 'HeartConfig::Instance()->GetCapacitance()'
        else:
            problem("Units of the ionic current are not in the "
                    "dimensions expected by Chaste (uA/cm^2) and cannot "
                    "be converted automatically.")
        return conversion, nodes_used
    
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
            elif not ms.equals(time_units):
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
        
        If self.dynamically_loadable is set, includes extra headers needed
        for that case.
        
        Reads self.include_serialization and self.use_backward_euler.
        Sets self.base_class_name and self.class_inheritance.
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
            self.writeln_hpp('#include "ChasteSerialization.hpp"')
            self.writeln_hpp('#include <boost/serialization/base_object.hpp>')
        self.writeln('#include <cmath>')
        self.writeln('#include <cassert>')
        self.writeln('#include <memory>')
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
        if self.use_modifiers:
            self.writeln_hpp('#include "AbstractCardiacCellWithModifiers.hpp"')
            self.writeln_hpp('#include "AbstractModifier.hpp"')
            # Modify the base class name
            self.base_class_name = 'AbstractCardiacCellWithModifiers<' + self.base_class_name + ' >'
        self.class_inheritance = ' : public ' + self.base_class_name
        if self.dynamically_loadable:
            self.writeln_hpp('#include "AbstractDynamicallyLoadableEntity.hpp"')
            self.class_inheritance += ', public AbstractDynamicallyLoadableEntity'
        if self.use_protocol:
            self.writeln_hpp('#include "AbstractSystemWithOutputs.hpp"')
            self.class_inheritance += ', public AbstractSystemWithOutputs<' + self.TYPE_VECTOR + '>'
        self.writeln('#include "Exception.hpp"')
        self.writeln('#include "OdeSystemInformation.hpp"')
        self.writeln('#include "RegularStimulus.hpp"')
        self.writeln_hpp('#include "AbstractStimulusFunction.hpp"')
        self.writeln('#include "HeartConfig.hpp"')
        self.writeln('#include "IsNan.hpp"')
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
        args_string = ', '.join(filter(None, map(str, args)))
        self.writeln_hpp(ret_type, method_name, '(', args_string, ')', self.STMT_END)
        self.writeln(ret_type, self.class_name, '::', method_name, '(', args_string, ')')

    def output_derived_quantities(self):
        """Output a ComputeDerivedQuantities method if any such quantities exist.
        
        Looks for variables annotated with pycml:derived-quantity=yes, and generates
        a method to compute all these variables from a given state.
        """
        dqs = self.derived_quantities
        if dqs:
            self.output_method_start('ComputeDerivedQuantities',
                                     [self.TYPE_DOUBLE + self.code_name(self.free_vars[0]),
                                      'const ' + self.TYPE_VECTOR + '& rY'], # We need it to really be a reference
                                     self.TYPE_VECTOR, access='public')
            self.open_block()
            self.output_comment('Inputs:')
            self.output_comment('Time units: ', self.free_vars[0].units)
            #621: convert if free var is not in milliseconds
            if self.conversion_factor:
                self.writeln(self.code_name(self.free_vars[0]), ' *= ',
                             self.conversion_factor, self.STMT_END)
            # Work out what equations are needed
            nodeset = self.calculate_extended_dependencies(dqs)
            # State variable inputs
            self.output_state_assignments(assign_rY=False, nodeset=nodeset)
            self.writeln()
            if self.use_lookup_tables:
                self.output_table_index_generation(
                    indexes_as_member=self.use_backward_euler,
                    nodeset=nodeset)
            # Output equations
            self.output_comment('Mathematics')
            self.output_equations(nodeset)
            self.writeln()
            # Assign to results vector
            self.writeln(self.vector_create('dqs', len(dqs)))
            for i, var in enumerate(dqs):
                # TODO: #621 conversions
                self.writeln(self.vector_index('dqs', i), self.EQ_ASSIGN,
                             self.code_name(var), self.STMT_END)
            self.writeln('return dqs', self.STMT_END)
            self.close_block(blank_line=True)
        
       
    def output_cell_parameters(self):
        """Output declarations, set & get methods for cell parameters.
        
        Sets self.cell_parameters to be those constant variables annotated with
        pycml:modifiable-parameter.  These use the mParameters functionality in
        Chaste.
        
        Also handles variables annotated with pe:keep (see #666).
        They can be real parameters, which have both set & get methods,
        or computed variables, which just have get methods.
        
        Also collects any variables annotated with an RDF oxmeta name into
        self.metadata_vars. Only constants and state variables are included.
        """
        # Find annotated parameters
        self.cell_parameters = filter(
            lambda v: v.is_modifiable_parameter,
            cellml_metadata.find_variables(self.model,
                                           ('pycml:modifiable-parameter', NSS['pycml']),
                                           'yes'))
        if self.kept_vars_as_members:
            kept_vars = cellml_metadata.find_variables(self.model, ('pe:keep', NSS[u'pe']), 'yes')
            kept_vars = list(set(kept_vars) - set(self.cell_parameters))
        else:
            kept_vars = []
        # Reduce intra-run variation
        kept_vars.sort(key=lambda v: v.fullname())
        self.cell_parameters.sort(key=lambda v: v.fullname())
        
        for i, var in enumerate(self.cell_parameters):
            # Remember the var's index
            var._cml_param_index = i

        # Create set of all oxmeta-annotated variables
        vars = cellml_metadata.find_variables(self.model, ('bqbiol:is', NSS[u'bqbiol']))
        # Keep only the variables with an oxmeta name
        vars = filter(lambda v: v.oxmeta_name, vars)
        # We're actually only interested in constants or state variables
        self.metadata_vars = set([v for v in vars if v.get_type() == VarTypes.Constant
                                  or v in self.state_vars])
        
        # #1464 Create a set of metadata variables that will have modifiers
        # We want to avoid writing out metadata for stimulus current as it is used once and then discarded.
        # \todo - use protocol information to put only the required modifiers into this list.
        stimulus_names = set('membrane_stimulus_current_'+ v for v in ['duration', 'amplitude', 'period', 'offset'])
        self.modifier_vars = set([v for v in self.metadata_vars if v.oxmeta_name not in stimulus_names])

        #1544: temporary check for possible Bad Things
        if (((self.use_modifiers and self.modifier_vars) or self.cell_parameters) and
            self.use_backward_euler):
            msg = "Backward Euler cell models cannot safely use modifiable parameters at present (see #1544)."
            print >>sys.stderr, msg # Really want raise ConfigurationError(msg) but this breaks things!

        # Generate member variable declarations
        self.set_access('private')
        if kept_vars or self.metadata_vars:
            self.output_comment('\nSettable parameters and readable variables\n',
                                subsidiary=True)
        for var in kept_vars:
            self.writeln_hpp(self.TYPE_DOUBLE, self.code_name(var), self.STMT_END)
            
        # Write out the modifier member variables. 
        if self.use_modifiers:
            for var in self.modifier_vars:
                self.writeln_hpp('boost::shared_ptr<AbstractModifier> mp_' + var.oxmeta_name + '_modifier', self.STMT_END)    
        # Generate Set & Get methods
        self.set_access('public')
        for var in kept_vars:
            # Generate Get method
            self.output_method_start('Get_' + self.code_name(var, prefix=''), [], self.TYPE_DOUBLE)
            self.open_block()
            self.writeln('return ', self.code_name(var), ';')
            self.close_block()
        
        # Methods associated with oxmeta annotated variables
        for var in self.metadata_vars:
            if var.get_type() == VarTypes.Constant:
                self.output_method_start('Get_' + var.oxmeta_name + '_constant', [], self.TYPE_DOUBLE)
                self.open_block()
                self.output_comment('Constant value given in CellML')
                self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(var), self.EQ_ASSIGN,
                             var.initial_value, self.STMT_END)
                conversion, nodes_used = self.determine_units_conversion(var)
                if nodes_used:
                    nodeset = self.calculate_extended_dependencies(nodes_used)
                    self.output_equations(nodeset)
                self.writeln('return ', conversion, self.STMT_END)
                self.close_block()
                self.writeln()
        self.output_default_stimulus()
        self.output_intracellular_calcium()
        
        # Find & store derived quantities, for use elsewhere
        self.derived_quantities = cellml_metadata.find_variables(self.model,
                                                                 ('pycml:derived-quantity', NSS['pycml']),
                                                                 'yes')
        # Reduce intra-run variation
        self.derived_quantities.sort(key=lambda v: v.fullname())
                
    def output_default_stimulus(self):
        """
        Output a default cell stimulus from the metadata specification
        as long as the following metadata exists:
         * membrane_stimulus_current_amplitude
         * membrane_stimulus_current_duration
         * membrane_stimulus_current_period
        and optionally:
         * membrane_stimulus_current_offset
        
        Ensures that the amplitude of the generated RegularStimulus is negative.
        """
        mandatory_args = set('membrane_stimulus_current_'+v for v in ['duration', 'amplitude', 'period'])
        optional_arg = 'membrane_stimulus_current_offset'
        var_names = set(v.oxmeta_name for v in self.metadata_vars)
        if len(var_names & mandatory_args) != 3:
            return

        self.output_method_start('UseCellMLDefaultStimulus', [], self.TYPE_VOID, 'public')
        self.open_block()
        self.output_comment('Use the default stimulus specified by CellML metadata')
        self.writeln('mpIntracellularStimulus.reset(new RegularStimulus(')
        self.writeln('        -fabs(Get_membrane_stimulus_current_amplitude_constant()),')
        self.writeln('        Get_membrane_stimulus_current_duration_constant(),')
        self.writeln('        Get_membrane_stimulus_current_period_constant(),')
        if optional_arg in var_names:
            self.writeln('        Get_membrane_stimulus_current_offset_constant()));')
        else :
            self.writeln('        0.0));')
        self.close_block(blank_line=True)
    
    def output_intracellular_calcium(self):
        """
        If a (state) variable has been annotated as cytosolic_calcium_concentration,
        generate a GetIntracellularCalciumConcentration method.
        """
        # Find cytosolic_calcium_concentration
        cai = self.doc.model.get_variable_by_oxmeta_name('cytosolic_calcium_concentration', throw=False)
        if cai and cai in self.state_vars:
            i = self.state_vars.index(cai[0])
            self.output_method_start('GetIntracellularCalciumConcentration', [], self.TYPE_DOUBLE, 'public')
            self.open_block()
            self.writeln('return ', self.vector_index('mStateVariables', i), self.STMT_END)
            self.close_block(blank_line=True)
        
    def code_name(self, var, *args, **kwargs):
        """
        Return the full name of var in a form suitable for inclusion in a
        source file.
        
        Overrides the base class version to access mParameters for parameters.
        """
        if hasattr(var, '_cml_param_index'):
            return self.vector_index('mParameters', var._cml_param_index)
        else:
            return super(CellMLToChasteTranslator, self).code_name(var, *args, **kwargs)

    def output_top_boilerplate(self):
        """Output top boilerplate.
        
        This method outputs the constructor and destructor of the cell
        class, and also lookup table declarations and lookup methods.
        It also calls output_verify_state_variables.
        """
        self.include_serialization = not self.use_modifiers # TODO: Implement
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
        self.writeln_hpp('class ', self.class_name, self.class_inheritance)
        self.open_block(subsidiary=True)
        # Serialization
        if self.include_serialization:
            self.writeln_hpp('friend class boost::serialization::access;')
            self.writeln_hpp('template<class Archive>')
            self.writeln_hpp('void serialize(Archive & archive, const unsigned int version)')
            self.open_block(subsidiary=True)
            self.writeln_hpp('archive & boost::serialization::base_object<', self.base_class_name,
                             ' >(*this);')
            if self.dynamically_loadable:
                self.writeln_hpp('archive & boost::serialization::base_object<AbstractDynamicallyLoadableEntity>(*this);')
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
                    self.set_access('private')
                    self.output_lut_indices()
                    self.set_access('public')
            else:
                self.output_lut_declarations()
                self.output_lut_row_lookup_memory()
                self.output_lut_methods()
            self.send_main_output_to_subsidiary(False)
        self.output_verify_state_variables()
        return
    
    def output_verify_state_variables(self):
        """Output the VerifyStateVariables method.
        
        This will look for state variables annotated with pycml:range-low and/or pycml:range-high,
        which specify allowable ranges for these variables.  The generated method will check that
        they are within the range.  Both limits are included, i.e. they specify a closed interval.
        """
        # Verify state variables method; empty at present
        self.output_method_start('VerifyStateVariables', [], 'void')
        self.open_block()
        low_prop = ('pycml:range-low', NSS['pycml'])
        high_prop = ('pycml:range-high', NSS['pycml'])
        low_range_vars = filter(
            lambda v: v.get_type() == VarTypes.State,
            cellml_metadata.find_variables(self.model, low_prop))
        high_range_vars = filter(
            lambda v: v.get_type() == VarTypes.State,
            cellml_metadata.find_variables(self.model, high_prop))
        nodeset = set(low_range_vars + high_range_vars)
        self.output_state_assignments(nodeset=nodeset)
        error_template = 'EXCEPTION(DumpState("State variable %s has gone out of range. Check model parameters, for example spatial stepsize"));'
        for var in low_range_vars:
            self.writeln('if (', self.code_name(var), ' < ', var.get_rdf_annotation(low_prop), ')')
            self.writeln(error_template % self.code_name(var), indent_offset=1)
        for var in high_range_vars:
            self.writeln('if (', self.code_name(var), ' > ', var.get_rdf_annotation(high_prop), ')')
            self.writeln(error_template % self.code_name(var), indent_offset=1)
        self.close_block(True)
  
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
        # Filter out empty params, to make backward Euler happy
        base_class_params = filter(None, map(str, base_class_params))
        for i, param in enumerate(base_class_params):
            if i == len(base_class_params)-1: comma = ')'
            else: comma = ','
            self.writeln(param, comma, indent_offset=3)
        self.open_block()
        self.output_comment('Time units: ', self.free_vars[0].units, '\n')
        self.writeln('this->mpSystemInfo = OdeSystemInformation<',
                     self.class_name, '>::Instance();')
        self.writeln('Init();\n')
        #1464 - cleverer modifiers...
        if self.use_modifiers and self.modifier_vars:
            self.output_comment('These will get initialised to DummyModifiers in the base class method.')
            for var in self.modifier_vars:
                self.writeln('this->AddModifier("' + var.oxmeta_name + '",')
                self.writeln('                  mp_' + var.oxmeta_name + '_modifier)', self.STMT_END)        
        #666 - initialise parameters
        for var in self.cell_parameters:
            if var.get_type() == VarTypes.Constant:
                self.writeln(self.vector_index('this->mParameters', var._cml_param_index),
                             self.EQ_ASSIGN, var.initial_value, self.STMT_END, ' ',
                             self.COMMENT_START, var.fullname(), ' [', var.units, ']')
        #1354 - specify protocol outputs
        if self.use_protocol:
            outputs = cellml_metadata.find_variables(self.model,
                                                     ('pycml:output-variable', NSS['pycml']),
                                                     'yes')
            if outputs:
                self.output_comment('Protocol outputs')
                self.writeln(self.vector_initialise('this->mOutputsInfo', len(outputs)))
                for i, output in enumerate(outputs):
                    self.writeln(self.vector_index('this->mOutputsInfo', i), self.EQ_ASSIGN,
                                 'std::make_pair(', nl=False)
                    if output.get_type() == VarTypes.Free:
                        self.writeln('UNSIGNED_UNSET, FREE', indent=False, nl=False)
                    elif output.get_type() == VarTypes.State:
                        self.writeln(self.state_vars.index(output), ', STATE', indent=False, nl=False)
                    elif output.is_derived_quantity:
                        self.writeln(self.derived_quantities.index(output), ', DERIVED', indent=False, nl=False)
                    elif output.is_modifiable_parameter:
                        self.writeln(self.cell_parameters.index(output), ', PARAMETER', indent=False, nl=False)
                    else:
                        raise ValueError('Unexpected protocol output: ' + str(output))
                    self.writeln(')', self.STMT_END, indent=False)
                self.writeln()
        # Lookup table generation, if not in a singleton
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
        self.writeln('if (mpInstance.get() == NULL)')
        self.writeln('{')
        self.writeln('mpInstance.reset(new ', self.lt_class_name, ');',
                     indent_offset=1)
        self.writeln('}')
        self.writeln('return mpInstance.get();')
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
        self.writeln('assert(mpInstance.get() == NULL);')
        self.output_lut_generation()
        self.close_block()
        # Private data
        self.writeln('private:', indent_level=0)
        self.writeln('/** The single instance of the class */')
        self.writeln('static std::auto_ptr<', self.lt_class_name, '> mpInstance;\n')
        if self.row_lookup_method:
            self.output_lut_row_lookup_memory()
        self.output_lut_declarations()
        # Close the class
        self.set_indent(0)
        self.writeln('};\n')
        # Define the instance pointer
        self.writeln('std::auto_ptr<', self.lt_class_name, '> ',
                     self.lt_class_name, '::mpInstance;')
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
        used_vars = set()
        for var in self.state_vars:
            if ((not exclude_nonlinear or 
                 self.code_name(var) not in self.nonlinear_system_vars)
                 and (nodeset is None or var in nodeset)):
                used_vars.add(var)
        if assign_rY and used_vars:
            self.writeln(self.TYPE_VECTOR_REF, 'rY = rGetStateVariables();')
        for i, var in enumerate(self.state_vars):
            if var in used_vars:
                if self.use_modifiers and var in self.modifier_vars:
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
        return ('mp_' + var.oxmeta_name + '_modifier->Calc(' +
                current_value + ', ' + self.code_name(self.free_vars[0]) + ')')
    
    def vector_index(self, vector, i):
        """Return code for accessing the i'th index of vector."""
        return vector + '[' + str(i) + ']'
    
    def vector_create(self, vector, size):
        """Return code for creating a new vector with the given size."""
        return ''.join(map(str, [self.TYPE_VECTOR, vector, '(', size, ')', self.STMT_END]))
    
    def vector_initialise(self, vector, size):
        """Return code for creating an already-declared vector with the given size."""
        return ''.join(map(str, [vector, '.resize(', size, ')', self.STMT_END]))

    def output_nonlinear_state_assignments(self, nodeset=None):
        """Output assignments for nonlinear state variables."""
        if nodeset:
            var_names = set([self.code_name(node) for node in nodeset if isinstance(node, cellml_variable)])
        for i, varname in enumerate(self.nonlinear_system_vars):
            if not nodeset or varname in var_names:
                self.writeln('double ', varname, self.EQ_ASSIGN,
                             self.vector_index('rCurrentGuess', i), self.STMT_END)
                #621 TODO: maybe convert if state var dimensions include time
        self.writeln()
        return
    
    def get_stimulus_assignment(self):
        """Return code for getting Chaste's stimulus current, with units conversions.
        
        Returns a tuple, the first item of which is the code, and the second is a
        set of nodes required for computing the conversion.
        """
        expr = self.doc._cml_config.i_stim_var
        output = self.code_name(expr) + self.EQ_ASSIGN
        #621: convert if free var is not in milliseconds
        if self.conversion_factor:
            conv_time = '(1.0/%.12g)*' % self.conversion_factor
        else:
            conv_time = ''
        get_stim = 'GetIntracellularAreaStimulus(' + conv_time + self.code_name(self.free_vars[0]) + ')'
        # Convert from Chaste stimulus units (uA/cm2) to model units
        stim_units = self.get_var_units(expr)
        conversion, conv_nodes = self.ionic_current_units_conversion(get_stim, stim_units, False)
        return output + conversion + self.STMT_END, conv_nodes

    def output_equations(self, nodeset, zero_stimulus=False):
        """Output the mathematics described by nodeset.

        nodeset represents a subset of the assignments in the model.
        Output assignments in the order given by a topological sort,
        but only include those in nodeset.
        """
        # Special case for the stimulus current
        if self.doc._cml_config.i_stim_var in nodeset:
            if zero_stimulus:
                i_stim = self.doc._cml_config.i_stim_var
                stim_assignment = self.code_name(i_stim) + self.EQ_ASSIGN + '0.0' + self.STMT_END
            else:
                stim_assignment, conv_nodes = self.get_stimulus_assignment()
                conv_nodes = self.calculate_extended_dependencies(conv_nodes)
                conv_nodes_new = conv_nodes - nodeset
                if conv_nodes_new:
                    if conv_nodes != conv_nodes_new:
                        # Some of nodeset is used to compute the conversion; ordering becomes tricky
                        stim_index = self.model.get_assignments().index(self.doc._cml_config.i_stim_var)
                        for node in conv_nodes - conv_nodes_new:
                            if self.model.get_assignments().index(node) > stim_index:
                                raise NotImplemented
                    self.output_equations(conv_nodes)
        for expr in (e for e in self.model.get_assignments() if e in nodeset):
            # Special-case the stimulus current
            if self.use_chaste_stimulus:
                if isinstance(expr, cellml_variable) and \
                        expr is self.doc._cml_config.i_stim_var:
                    clear_type = (self.kept_vars_as_members and expr.pe_keep)
                    if clear_type:
                        self.TYPE_CONST_DOUBLE = ''
                    self.writeln(self.TYPE_CONST_DOUBLE, stim_assignment)
                    if clear_type:
                        # Remove the instance attribute, thus reverting to the class member
                        del self.TYPE_CONST_DOUBLE
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
        
        Also has overrides for modifiable parameters and modifier calls.
        """
        clear_type = False
        # Figure out what is being assigned to
        if isinstance(expr, cellml_variable):
            assigned_var = expr
        else:
            if expr.eq.lhs.localName == 'ci':
                assigned_var = expr.eq.lhs.variable
            else:
                assigned_var = None # We don't store derivatives as members
                #907: Check if this is the derivative of the transmembrane potential
                if not self.use_backward_euler and expr.eq.lhs.diff.dependent_variable == self.v_variable:
                    clear_type = True
        # Parameters don't need assigning
        if assigned_var in self.cell_parameters:
            return
        # Is the variable assigned to stored as a class member?
        clear_type = (clear_type or
                      (self.kept_vars_as_members and assigned_var and assigned_var.pe_keep
                       and not assigned_var.is_derived_quantity))
        if clear_type:
            self.TYPE_DOUBLE = self.TYPE_CONST_DOUBLE = ''
        if (assigned_var and assigned_var.get_type() == VarTypes.Constant and
            self.use_modifiers and assigned_var in self.modifier_vars):
            # "Constant" oxmeta-annotated parameters may be modified at run-time
            self.writeln(self.TYPE_CONST_DOUBLE, self.code_name(assigned_var), self.EQ_ASSIGN,
                         self.modifier_call(assigned_var, expr.initial_value), self.STMT_END, nl=False)
            self.output_comment(assigned_var.units, indent=False, pad=True)
        else:
            super(CellMLToChasteTranslator, self).output_assignment(expr)
        if clear_type:
            # Remove the instance attributes, thus reverting to the class members
            del self.TYPE_DOUBLE
            del self.TYPE_CONST_DOUBLE
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
        
        Where derived-quantity annotations are present, we also generate a
        ComputeDerivedQuantities method.
        """
        self.output_get_i_ionic()
        if self.use_backward_euler:
            self.output_backward_euler_mathematics()
        else:
            self.output_evaluate_y_derivatives()
        self.output_derived_quantities()

    def output_get_i_ionic(self):
        """Output the GetIIonic method."""
        use_modifiers = self.use_modifiers
        self.use_modifiers = False
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
            # GetIIonic must not include the stimulus current
            i_stim = self.doc._cml_config.i_stim_var
            nodeset = self.calculate_extended_dependencies(nodes, prune_deps=[i_stim])
            #print map(lambda v: v.fullname(), nodes)
            #print filter(lambda p: p[2]>0, map(debugexpr, nodeset))
            #621: check units of the ionic current
            units_objs = map(self.get_var_units, nodes)
            conversion, conv_nodes = self.ionic_current_units_conversion('i_ionic', units_objs)
            conv_nodes = self.calculate_extended_dependencies(conv_nodes, prune=nodeset, prune_deps=[i_stim])
            all_nodes = conv_nodes|nodeset
            # Output main part of maths
            self.output_state_assignments(nodeset=all_nodes)
            if self.use_lookup_tables:
                self.output_table_index_generation(
                    indexes_as_member=self.use_backward_euler,
                    nodeset=all_nodes)
            self.output_equations(nodeset, zero_stimulus=True)
            self.writeln()
            # Assign the total current to a temporary so we can check for NaN and
            # do units conversion if needed.
            self.writeln(self.TYPE_DOUBLE, 'i_ionic', self.EQ_ASSIGN, nl=False)
            if self.doc._cml_config.i_ionic_negated:
                self.writeln('-(', nl=False, indent=False)
            plus = False
            for varelt in self.model.solver_info.ionic_current.var:
                if plus: self.write('+')
                else: plus = True
                self.output_variable(varelt)
            if self.doc._cml_config.i_ionic_negated:
                self.writeln(')', nl=False, indent=False)
            self.writeln(self.STMT_END, indent=False)
            self.writeln('EXCEPT_IF_NOT(!std::isnan(i_ionic));')
            self.output_equations(conv_nodes)
            self.writeln('return ', conversion, self.STMT_END)
        else:
            self.writeln('return 0.0;')
        self.close_block()
        self.use_modifiers = use_modifiers

    def output_evaluate_y_derivatives(self, method_name='EvaluateYDerivatives'):
        """Output the EvaluateYDerivatives method."""
        # Start code output
        self.output_method_start(method_name,
                                 [self.TYPE_DOUBLE + self.code_name(self.free_vars[0]),
                                  'const ' + self.TYPE_VECTOR_REF + 'rY',
                                  self.TYPE_VECTOR_REF + 'rDY'],
                                 'void', access='public')
        self.open_block()
        if not self.state_vars:
            # This isn't an ODE model!
            self.close_block()
            return
        self.output_comment('Inputs:')
        self.output_comment('Time units: ', self.free_vars[0].units)
        #621: convert if free var is not in milliseconds
        if self.conversion_factor:
            self.writeln(self.code_name(self.free_vars[0]), ' *= ',
                         self.conversion_factor, self.STMT_END)
        # Work out what equations are needed to compute the derivatives
        derivs = set(map(lambda v: (v, self.free_vars[0]), self.state_vars))
        if self.v_variable in self.state_vars:
            dvdt = (self.v_variable, self.free_vars[0])
            derivs.remove(dvdt) #907: Consider dV/dt separately
        else:
            dvdt = None
        if self.use_chaste_stimulus:
            i_stim = [self.doc._cml_config.i_stim_var]
        else:
            i_stim = []
        nonv_nodeset = self.calculate_extended_dependencies(
            derivs, prune_deps=i_stim)
        if dvdt:
            v_nodeset = self.calculate_extended_dependencies(
                [dvdt], prune=nonv_nodeset, prune_deps=i_stim)
        else:
            v_nodeset = set()
        # State variable inputs
        all_nodes = nonv_nodeset|v_nodeset
        self.output_state_assignments(assign_rY=False, nodeset=all_nodes)
        self.writeln()
        if self.use_lookup_tables:
            # TODO: Filter tables by available state/protocol variables?
            self.output_table_index_generation(
                indexes_as_member=self.use_backward_euler,
                nodeset=all_nodes)
        self.output_comment('Mathematics')
        #907: Declare dV/dt
        if dvdt:
            self.writeln(self.TYPE_DOUBLE, self.code_name(self.v_variable, ode=True), self.STMT_END)
        # Output mathematics required for non-dV/dt derivatives (which may include dV/dt)
        self.output_equations(nonv_nodeset)
        self.writeln()
        #907: Calculation of dV/dt
        if dvdt:
            self.writeln('if (mSetVoltageDerivativeToZero)')
            self.open_block()
            self.writeln(self.code_name(self.v_variable, ode=True), self.EQ_ASSIGN, '0.0', self.STMT_END)
            self.close_block(blank_line=False)
            self.writeln('else')
            self.open_block()
            self.output_equations(v_nodeset)
            self.close_block()
        # Assign to derivatives vector
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
                                 [self.TYPE_DOUBLE + self.code_name(self.free_vars[0]),
                                  self.TYPE_CONST_DOUBLE + 'rCurrentGuess' + argsize,
                                  self.TYPE_DOUBLE + 'rResidual' + argsize],
                                 'void', access='public')
        self.open_block()
        # Output mathematics for computing du/dt for each nonlinear state var u
        nodes = map(lambda u: (self.varobj(u), self.free_vars[0]),
                    self.nonlinear_system_vars)
        nodeset = self.calculate_extended_dependencies(nodes, prune_deps=[self.doc._cml_config.i_stim_var])
        self.output_state_assignments(exclude_nonlinear=True, nodeset=nodeset)
        self.output_nonlinear_state_assignments(nodeset=nodeset)
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
                                 self.code_name(var, ode=True), self.STMT_END)
                else:
                    self.writeln('rResidual[', j, '] = rCurrentGuess[', j,
                                 '] - rY[', i, '] - mDt*',
                                 self.code_name(var, ode=True), self.STMT_END)
        self.close_block()
        # Jacobian
        ##########
        self.output_method_start('ComputeJacobian',
                                 [self.TYPE_DOUBLE + self.code_name(self.free_vars[0]),
                                  self.TYPE_CONST_DOUBLE + 'rCurrentGuess' + argsize,
                                  self.TYPE_DOUBLE + 'rJacobian' + argsize + argsize],
                                 'void', access='public')
        self.open_block()
        # Mathematics that the Jacobian depends on
        used_vars = set()
        for entry in self.model.solver_info.jacobian.entry:
            used_vars.update(self._vars_in(entry.math))
        nodeset = self.calculate_extended_dependencies(used_vars, prune_deps=[self.doc._cml_config.i_stim_var])
        self.output_state_assignments(exclude_nonlinear=True, nodeset=nodeset)
        self.output_nonlinear_state_assignments(nodeset=nodeset)
        if self.conversion_factor:
            self.writeln('const double dt = ', self.conversion_factor, ' * mDt;\n');
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
            elif hasattr(entry.math, u'cn'):
                self.output_expr(entry.math.cn, False)
            elif hasattr(entry.math, u'ci'):
                self.output_expr(entry.math.ci, False)
            else:
                raise ValueError('Unexpected entry: ' + entry.xml())
            self.writeln(self.STMT_END, indent=False)
        self.close_block()
        # The other methods are protected
        self.writeln_hpp('protected:', indent_offset=-1)
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
        nodeset = self.calculate_extended_dependencies(nodes, prune_deps=[self.doc._cml_config.i_stim_var])
        self.output_state_assignments(nodeset=nodeset)
        if self.use_lookup_tables:
            self.output_table_index_generation(indexes_as_member=True,
                                               nodeset=nodeset)
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
        used_vars = set() # NB: Also contains g&h if they are mathml_applys so table index generation works
        for ode in self.model.solver_info.linear_odes.math.apply:
            varname = unicode(ode.apply.ci)
            g = ode.apply[1].operands().next()
            hu = list(ode.apply[1].operands())[1]
            h = hu.operands().next()
            linear_vars.append(self.varobj(varname))
            ghs.append((g, h))
            if not isinstance(g, mathml_cn): used_vars.add(g)
            if not isinstance(h, mathml_cn): used_vars.add(h)
            used_vars.update(self._vars_in(g))
            used_vars.update(self._vars_in(h))
        # Output required equations for used variables
        nodeset = self.calculate_extended_dependencies(used_vars, prune_deps=[self.doc._cml_config.i_stim_var])
        self.output_state_assignments(nodeset=nodeset)
        if self.use_lookup_tables:
            self.output_table_index_generation(indexes_as_member=True,
                                               nodeset=nodeset)
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
        self.writeln('_solver->Solve(*this, ', self.code_name(self.free_vars[0]), ', _guess);')
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
        self.writeln('this->mSystemName', self.EQ_ASSIGN, '"', self.model.name, '"', self.STMT_END)
        self.writeln()
        def output_var(vector, var):
            if var.oxmeta_name:
                self.writeln('this->m', vector, 'Names.push_back("', var.oxmeta_name, '");')   
            elif hasattr(var, u'id') and var.id:
                self.writeln('this->m', vector, 'Names.push_back("', var.id, '");') 
            else:
                self.writeln('this->m', vector, 'Names.push_back("', var.name, '");')
            self.writeln('this->m', vector, 'Units.push_back("', var.units, '");')
        for var in self.state_vars:
            output_var('Variable', var)
            init_val = getattr(var, u'initial_value', None)
            if init_val is None:
                init_comm = ' // Value not given in model'
                # Don't want compiler error, but shouldn't be a real number
                init_val = 'NAN'
            else:
                init_comm = ''
            self.writeln('this->mInitialConditions.push_back(', init_val, ');',
                       init_comm, '\n')
        # Model parameters
        for var in self.cell_parameters:
            if var.get_type() == VarTypes.Constant:
                output_var('Parameter', var)
                self.writeln()
        # Derived quantities
        for var in self.derived_quantities:
            output_var('DerivedQuantity', var)
        self.writeln('this->mInitialised = true;')
        self.close_block()
        self.writeln()
        # Serialization
        if self.include_serialization:
            self.output_comment('Needs to be included last', subsidiary=True)
            self.writeln_hpp('#include "SerializationExportWrapper.hpp"')
            self.writeln_hpp('CHASTE_CLASS_EXPORT(', self.class_name, ')')
            self.output_comment('Serialization for Boost >= 1.36')
            self.writeln('#include "SerializationExportWrapperForCpp.hpp"')
            self.writeln('CHASTE_CLASS_EXPORT(', self.class_name, ')')
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
            self.writeln('AbstractCardiacCellInterface* MakeCardiacCell(')
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
            # This ci element doesn't have all the extra annotations.  It is a fully
            # qualified name though.  This is typically because PE has been done.
            prefix = ['var_', 'd_dt_'][ode]
            varname = unicode(ci_elt)
            try:
                var = self.varobj(varname)
            except KeyError:
                var = None
            if var:
                self.write(self.code_name(var, ode=ode))
            else:
                if varname == u'delta_t':
                    # Special case for the timestep in ComputeJacobian
                    prefix = ''
                    varname = 'dt'
                else:
                    # Assume it's a suitable name
                    pass
                self.write(prefix + varname)
        return


class CellMLToCvodeTranslator(CellMLToChasteTranslator):
    """Translate a CellML model to C++ code for use with Chaste+CVODE."""
    
    # Type of (a reference to) the state variable vector
    TYPE_VECTOR = 'N_Vector '
    TYPE_VECTOR_REF = 'N_Vector ' # CVODE's vector is actually a pointer type
    
    def vector_index(self, vector, i):
        """Return code for accessing the i'th index of vector."""
        return 'NV_Ith_S(' + vector + ', ' + str(i) + ')'
    
    def vector_create(self, vector, size):
        """Return code for creating a new vector with the given size."""
        return ''.join(map(str, [self.TYPE_VECTOR, vector, self.EQ_ASSIGN,
                                 'N_VNew_Serial(', size, ')', self.STMT_END]))

    def vector_initialise(self, vector, size):
        """Return code for creating an already-declared vector with the given size."""
        return ''.join(map(str, [vector, self.EQ_ASSIGN, 'N_VNew_Serial(', size, ')', self.STMT_END]))

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
        self.writeln_hpp('class ', self.class_name, self.class_inheritance)
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
        self.output_derived_quantities()
    
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
            for var_i in self.state_vars:
                for var_j in self.state_vars:
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




class SolverInfo(object):
    """Add information for specialised translator classes into a model."""
    def __init__(self, model, force=False):
        """Add information for the solvers as XML.

        The Jacobian and linearity analyses store their results in
        Python data structures as attributes of this object.
        Transcribe these into XML in a child <solver_info> element.

        If any of these elements exist in the model they will be left
        unaltered, unless force is set to True.
        
        This constructor just sets up the container element; call one
        of the add_* methods to actually add the information to it.
        """
        self._model = model
        if force and hasattr(model, u'solver_info'):
            model.xml_remove_child(model.solver_info)
        if hasattr(model, u'solver_info'):
            solver_info = model.solver_info
        else:
            solver_info = model.xml_create_element(u'solver_info', NSS[u'solver'])
            model.xml_append(solver_info)
        self._solver_info = solver_info
        self._component = None
    
    def add_all_info(self):
        """Actually add the info."""
        self.add_transmembrane_potential_name()
        self.add_membrane_ionic_current()
        self.add_linearised_odes()
        self.add_jacobian_matrix()
    
    def add_transmembrane_potential_name(self):
        """The name of the transmembrane potential."""
        solver_info = self._solver_info
        model = self._model
        if not hasattr(solver_info, u'transmembrane_potential'):
            v_elt = model.xml_create_element(
                u'transmembrane_potential', NSS[u'solver'],
                content=model._cml_transmembrane_potential.fullname())
            solver_info.xml_append(v_elt)
    
    def add_linearised_odes(self):
        """Linearised ODEs - where du/dt = g + hu (and g, h are not functions of u).
        
        Structure looks like:
        <linear_odes>
            <math>
                <apply><eq/>
                    <apply><diff/>
                        <bvar><ci>t</ci></bvar>
                        <ci>u</ci>
                    </apply>
                    <apply><plus/>
                        g
                        <apply><times/>
                            h
                            <ci>u</ci>
                        </apply>
                    </apply>
                </apply>
                .
                .
                .
            </math>
        </linear_odes>
        """
        solver_info = self._solver_info
        model = self._model
        if not hasattr(solver_info, u'linear_odes'):
            odes_elt = model.xml_create_element(u'linear_odes', NSS[u'solver'])
            solver_info.xml_append(odes_elt)
            odes_math = model.xml_create_element(u'math', NSS[u'm'])
            odes_elt.xml_append(odes_math)
            linear_vars = model._cml_linear_update_exprs.keys()
            linear_vars.sort(key=lambda v: v.fullname())
            free_var = model._cml_free_var
            for var in linear_vars:
                g, h = model._cml_linear_update_exprs[var]
                hu = mathml_apply.create_new(model, u'times', [h, var.fullname()])
                rhs = mathml_apply.create_new(model, u'plus', [g, hu])
                odes_math.xml_append(mathml_diff.create_new(
                    model, free_var.fullname(), var.fullname(), rhs))
    
    def add_jacobian_matrix(self):
        """Jacobian matrix elements.
        
        Structure looks like:
        <jacobian>
            <entry var_i='varname' var_j='varname'>
                <math> apply|cn|ci ...</math>
            </entry>
        </jacobian>
        """
        solver_info = self._solver_info
        model = self._model
        if not hasattr(solver_info, u'jacobian'):
            jac_elt = model.xml_create_element(u'jacobian', NSS[u'solver'])
            solver_info.xml_append(jac_elt)
            jac_vars = model._cml_jacobian.keys()
            jac_vars.sort() # Will sort by variable name
            rules = [bt.ws_strip_element_rule(u'*')]
            for v_i, v_j in jac_vars:
                # Add (i,j)-th entry
                binder = make_xml_binder()
                attrs = {u'var_i': unicode(v_i),
                         u'var_j': unicode(v_j)}
                entry = model.xml_create_element(u'entry', NSS[u'solver'], attributes=attrs)
                jac_elt.xml_append(entry)
                entry_doc = amara_parse(model._cml_jacobian[(v_i, v_j)].xml(),
                                        rules=rules, binderobj=binder)
                entry.xml_append(entry_doc.math)
        return

    def add_membrane_ionic_current(self):
        """Add ionic current information as XML for solvers to use."""
        solver_info = self._solver_info
        model = self._model
        # The total ionic current.  This relies on having a
        # configuration store.
        if hasattr(model.xml_parent, '_cml_config') and \
               not hasattr(solver_info, u'ionic_current'):
            conf = model.xml_parent._cml_config
            ionic_elt = model.xml_create_element(u'ionic_current', NSS[u'solver'])
            # Adds each ionic var to the xml doc from the config store
            for var in conf.i_ionic_vars:
                DEBUG("translate", var.name, var.xml_parent.name, var.fullname())
                varelt = model.xml_create_element(u'var', NSS[u'solver'],
                                                  content=var.fullname())
                ionic_elt.xml_append(varelt)
            solver_info.xml_append(ionic_elt)
        return
    
    def add_variable_links(self):
        """Link ci elements in the added XML to cellml_variable objects.
        
        This analyses the names in the ci elements to determine which variable in
        the model they refer to.
        """
        self._process_mathematics(self._add_variable_links)
    
    def do_binding_time_analysis(self):
        """Do a binding time analysis on the additional mathematics.
        
        This requires self.add_variable_links to have been called already.
        """
        self._process_mathematics(lambda elt: elt._get_binding_time())
        
    def _process_mathematics(self, func):
        """Apply func to each top-level mathematical construct in the solver info blocks.
        
        func must be able to accept mathml_apply, mathml_ci and mathml_cn elements.
        """
        solver_info = self._solver_info
        # Jacobian
        if hasattr(solver_info, u'jacobian'):
            for entry in solver_info.jacobian.entry.math:
                for elt in entry.xml_children:
                    if getattr(elt, 'nodeType', None) == Node.ELEMENT_NODE:
                        func(elt)
        # Linearised ODEs
        if hasattr(solver_info, u'linear_odes'):
            for elt in solver_info.linear_odes.math.xml_children:
                if getattr(elt, 'nodeType', None) == Node.ELEMENT_NODE:
                    func(elt)
    
    def get_modifiable_mathematics(self):
        """Get an iterable over mathematical constructs in the solver info blocks that can be changed.
        
        Returned elements will be mathml_apply, mathml_ci or mathml_cn instances.
        """
        solver_info = self._solver_info
        # Jacobian - entry definitions can be changed
        if hasattr(solver_info, u'jacobian'):
            for entry in solver_info.jacobian.entry.math:
                for elt in entry.xml_children:
                    if getattr(elt, 'nodeType', None) == Node.ELEMENT_NODE:
                        yield elt
        # Linearised ODEs - only g & h can be changed
        if hasattr(solver_info, u'linear_odes'):
            for ode in solver_info.linear_odes.math.apply:
                rhs = list(ode.operands())[1]
                opers = rhs.operands()
                g = opers.next()
                h = opers.next().operands().next()
                yield g
                yield h
    
    def _add_variable_links(self, elt):
        """Recursively link ci elements in the given XML tree to cellml_variable objects.
        
        Also sets component links: for ci elements, to the component containing the linked
        variable, and for cn elements, to the first component in the model.
        """
        if isinstance(elt, mathml_ci):
            var = self._get_variable(unicode(elt))
            elt._cml_variable = var
            elt._cml_component = var.component
        elif isinstance(elt, mathml_cn):
            # Fake a component, since it doesn't really have one
            elt._cml_component = elt.model.component
        elif hasattr(elt, 'xml_children'):
            for child in elt.xml_children:
                self._add_variable_links(child)

    def _get_variable(self, varname):
        """Return the variable in the model with name varname."""
        if varname[0] == '(':
            # (compname,varname)
            cname, vname = varname[1:-1].split(u',')
        elif '__' in varname:
            # [var_]cname__vname
            if varname.startswith('var_'):
                varname = varname[4:]
            cname, vname = varname.split(u'__')
        elif varname == u'delta_t':
            # Special case for the timestep in ComputeJacobian
            return self._get_special_variable(u'dt', VarTypes.Free)
        else:
            raise ValueError("Unrecognised variable name in SolverInfo: " + varname)
        # Determine the variable object from cname,vname
        try:
            comp = self._model.get_component_by_name(cname)
            var = comp.get_variable_by_name(vname)
        except KeyError:
            if len(list(self._model.component)) > 1:
                raise ValueError("Cannot find component '%s' needed by ci element '%s'"
                                 % (cname, varname))
            elif self._model.component.ignore_component_name:
                # Done PE already
                try:
                    var = self._model.component.get_variable_by_name(cname+'__'+vname)
                except KeyError:
                    raise ValueError("Cannot find variable '%s' in SolverInfo" % varname)
            else:
                raise ValueError("Cannot find variable '%s' in SolverInfo" % varname)
        return var

    def _get_special_variable(self, varname, ptype=VarTypes.Unknown):
        """Get or create a special variable object that doesn't really exist in the model."""
        comp = self._get_special_component()
        try:
            var = comp.get_variable_by_name(varname)
        except KeyError:
            var = cellml_variable.create_new(self._model, varname, u'dimensionless')
            comp._add_variable(var)
            var._set_type(ptype)
        return var

    def _get_special_component(self):
        """Get or create a special component for containing special variables."""
        if not self._component:
            self._component = cellml_component.create_new(self._model, u'')
            self._model._add_component(self._component, special=True)
        return self._component



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
        # Membrane capacitance
        self.Cm_definitions = []
        self.Cm_variable = None
        # Lookup table configuration
        self.lut_config = {}
        self.lut_config_keys = []
        # Ionic currents configuration
        self.i_stim_definitions = [u'membrane,i_Stim']
        self.i_stim_var = None
        self.i_ionic_definitions = [u'membrane,i_.*']
        self.i_ionic_vars = []
        # Whether GetIIonic will need to negate the sum of i_ionic_vars
        self.i_ionic_negated = False
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
             
             This is mostly redundant now, because the equation for dV/dt
             is used first to determine the ionic currents (see documentation
             for _find_transmembrane_currents_from_voltage_ode), and only
             if this fails to find suitable currents will the ionic_match
             definition be used.
         * 'transmembrane_potential'
           Defines which variable holds the transmembrane potential.
           Defaults to 'membrane,V' if not present.
         * 'membrane_capacitance'
           Defines which variable holds the cell membrane capacitance.
           
        The root element also contains 0 or more 'for_model' elements,
        which contain settings for individual models.  These must have
        at least one of an 'id' or 'name' attribute, which specify the
        model in question.  They can also contain anything allowable as
        global configuration options.  Options given here override those
        specified globally.

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
        3. By reference to a section of the config file (when defining lookup table keys)
           e.g. <var type='config-name'>transmembrane_potential</var>

        Within any element that specifies a variable, a list of <var> elements can be
        provided.  Each will be tried in turn to see if a match can be found in the model,
        and the first match wins.
        Alternatively these elements may have text content which must be of the form
        'component_name,variable_name'.

        Some items are overridden if oxmeta annotations are present in the model, with
        the annotated variable taking precedence over the config file specification.
        """
        DEBUG('config', "Reading configuration from ", config_file)
        rules = [bt.ws_strip_element_rule(u'*')]
        config_doc = amara_parse(config_file, rules=rules)
        # Overrides for command-line options
        if self.options and hasattr(config_doc.pycml_config, 'command_line_args'):
            args = map(str, config_doc.pycml_config.command_line_args.arg)
            args.append('dummy-file')
            get_options(args, self.options)
        # Sections to use in configuration; later sections take precedence
        sections = []
        # Use global configuration?
        glo = config_doc.xml_xpath(u'/*/global')
        if glo:
            sections.append(glo[0])
        # Get the config section(s) for our model.  Sections
        # specifically for this model come after sections covering
        # multiple models, so they take precedence.
        model_id = getattr(self.doc.model, u'id', self.doc.model.name)
        sections.extend(config_doc.xml_xpath(
            u'/*/for_models[ids/id="%s" or ids/id="%s"]'
            % (self.doc.model.name, model_id)))
        sections.extend(config_doc.xml_xpath(
            u'/*/for_model[@name="%s" or @id="%s"]'
            % (self.doc.model.name, model_id)))
        # Main items of configuration
        for section in sections:
            if hasattr(section, u'lookup_tables'):
                self._parse_lookup_tables(section.lookup_tables)
            if hasattr(section, u'currents'):
                self._parse_currents(section.currents)
            if hasattr(section, u'transmembrane_potential'):
                self._parse_Vm(section.transmembrane_potential)
            if hasattr(section, u'membrane_capacitance'):
                self._parse_Cm(section.membrane_capacitance)
    
    def finalize_config(self):
        """Having read all the configuration files, apply to the model."""
        # Identify the variables in the model
        self.find_transmembrane_potential()
        self.find_membrane_capacitance()
        self.find_current_vars()

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
            if unicode(var_elt) not in cellml_metadata.METADATA_NAMES:
                raise ConfigurationError('"' + unicode(var_elt) + '" is not a valid oxmeta name')
        elif defn_type == u'config-name':
            if unicode(var_elt) not in [u'stimulus', u'transmembrane_potential', u'membrane_capacitance']:
                raise ConfigurationError('"' + unicode(var_elt) + '" is not a name known to the config file')
        else:
            raise ConfigurationError('"' + defn_type + '" is not a valid variable definition type')
        return

    def _parse_var(self, elt, name):
        """Parse definition of a special variable."""
        if hasattr(elt, 'var'):
            # List of possibilities
            defs = []
            for vardef in elt.var:
                self._check_var_def(vardef, name)
                defs.append(vardef)
        else:
            # Old style - single variable given by text content
            self._check_var_def(elt, name)
            defs = [elt]
        return defs

    def _parse_Vm(self, vm_elt):
        """Parse definition of variable holding the transmembrane potential."""
        self.V_definitions = self._parse_var(vm_elt, 'transmembrane_potential')
    
    def _parse_Cm(self, cm_elt):
        """Parse definition of variable holding the cell membrane capacitance."""
        self.Cm_definitions = self._parse_var(cm_elt, 'membrane_capacitance')

    def _parse_currents(self, currents):
        """Parse definitions of ionic and stimulus currents."""
        if hasattr(currents, u'stimulus'):
            self.i_stim_definitions = self._parse_var(currents.stimulus, 'stimulus current')
        if hasattr(currents, u'ionic_match'):
            self.i_ionic_definitions = self._parse_var(currents.ionic_match, 'ionic currents')
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
            var = self.doc.model.get_variable_by_oxmeta_name(str(defn), throw=False)
        elif defn_type == u'config-name':
            if unicode(defn) == u'stimulus':
                var = self.i_stim_var
            elif unicode(defn) == u'transmembrane_potential':
                var = self.V_variable
            elif unicode(defn) == u'membrane_capacitance':
                var = self.Cm_variable
            else:
                raise ConfigurationError('"' + str(defn) + '" is not a valid configuration file variable name')
        else:
            raise ConfigurationError('"' + defn_type + '" is not a valid variable definition type')
        return var
    
    def _find_transmembrane_currents_from_voltage_ode(self):
        """Analyse the expression for dV/dt to determine the transmembrane currents.
        
        Looks for an equation defining dV/d(something) and assumes the something is
        time; this will be checked during code generation for Chaste.  It then finds
        all variables on the RHS of this equation which have the same units as the
        stimulus current (self.i_stim_var) and identifies these as transmembrane
        currents.  Will automatically exclude the stimulus current.
        
        If self.V_variable or self.i_stim_var are not set, returns the empty list.
        """
        if not self.V_variable:
            DEBUG('config', "Transmembrane potential not configured, so can't "
                  "determine currents from its ODE")
            return []
        if not self.i_stim_var:
            DEBUG('config', "Stimulus current not configured, so can't determine "
                  "transmembrane currents from dV/dt")
            return []
        stim_units = self.i_stim_var.component.get_units_by_name(self.i_stim_var.units)
        def search_expr(expr):
            """Recursively search for ci elements."""
            if isinstance(expr, mathml_ci):
                v = expr.variable.get_source_variable(recurse=True)
                if v is not self.i_stim_var:
                    # Check units
                    u = v.component.get_units_by_name(v.units)
                    if u.dimensionally_equivalent(stim_units):
                        ionic_vars.append(v)
                # Fake this variable being 1 so we can check the sign of GetIIonic
                expr.variable.set_value(1.0)
            elif isinstance(expr, mathml_apply):
                for o in expr.operands():
                    search_expr(o)
        def clear_values(expr):
            """Recursively clear saved values for variables in this expression."""
            if isinstance(expr, mathml_ci):
                expr.variable.unset_values()
            else:
                for elt in getattr(expr, 'xml_children', []):
                    clear_values(elt)
        ionic_vars = []
        # Iterate over all expressions in the model, to find the one for dV/d(something)
        for expr in (e for e in self.doc.model.get_assignments() if isinstance(e, mathml_apply) and e.is_ode()):
            # Assume the independent variable is time; if it isn't, we'll catch this later
            (dep_var, time_var) = expr.assigned_variable()
            if dep_var.get_source_variable(recurse=True) is self.V_variable:
                # Recursively search for ci elements
                search_expr(expr.eq.rhs)
                # Check the sign of the RHS
                self.i_ionic_negated = expr.eq.rhs.evaluate() > 0.0
                clear_values(expr.eq.rhs)
                # Found dV/d(something); don't check any more expressions
                break
        DEBUG('config', "Found ionic currents from dV/dt: ", ionic_vars)
        return ionic_vars
    
    def _find_var(self, oxmeta_name, definitions):
        """Find the variable object in the model for a particular concept.
        
        Will look for a variable annotated with the given oxmeta_name first, then
        try the list of definitions from the configuration file in turn.
        """
        var = None
        # Prepend an oxmeta definition
        oxmeta_defn = self._create_var_def(oxmeta_name, 'oxmeta')
        for defn in [oxmeta_defn] + definitions:
            var = self._find_variable(defn)
            if var:
                break
        return var

    def find_current_vars(self):
        """Find the variables representing currents."""
        # Find the stimulus current, if it exists for this kind of model (some are self-excitatory)
        if not self.doc.model.is_self_excitatory():
            self.i_stim_var = self._find_var('membrane_stimulus_current', self.i_stim_definitions)
            if not self.i_stim_var:
                # No match :(
                msg = "No stimulus current found; you'll have trouble generating Chaste code"
                if self.options.fully_automatic:
                    raise ConfigurationError(msg)
                else:
                    print >>sys.stderr, msg
        # For other ionic currents, try using the equation for dV/dt first
        self.i_ionic_vars = self._find_transmembrane_currents_from_voltage_ode()
        # Otherwise use the config file
        if not self.i_ionic_vars:
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
            msg = "No ionic currents found; you'll have trouble generating Chaste code"
            if self.options.fully_automatic:
                raise ConfigurationError(msg)
            else:
                print >>sys.stderr, msg
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
        def_dict = optimize.LookupTableAnalyser._LT_DEFAULTS
        for k, v in def_dict.iteritems():
            if k != 'table_var':
                lut_dict[k] = v
        return

    def annotate_currents_for_pe(self):
        """Annotate ionic & stimulus current vars so PE doesn't remove them.
        Also annotate the membrane capacitance, if defined."""
        if self.i_stim_var:
            self.i_stim_var.set_pe_keep(True)
        for var in self.i_ionic_vars:
            var.set_pe_keep(True)
        if self.Cm_variable:
            self.Cm_variable.set_pe_keep(True)
        return
    
    def annotate_metadata_for_pe(self):
        "Annotate all vars tagged with metadata so PE doesn't remove them."
        for var in self.metadata_vars:
            var.set_pe_keep(True)
        return

    def find_transmembrane_potential(self):
        """Find and store the variable object representing V.

        Tries metadata annotation first.  If that fails, uses the name given in
        the command line options, if present.  If that fails, uses the config file.
        """
        if not self.options:
            raise ValueError('No command line options given')
        # Check command line option before config file
        if self.options.transmembrane_potential:
            self.V_definitions[0:0] = [self.options.transmembrane_potential.strip().split(',')]
            if len(self.V_definitions[0]) != 2:
                raise ConfigurationError('The name of V must contain both '
                                         'component and variable name')
        self.V_variable = self._find_var('membrane_voltage', self.V_definitions)
        if not self.V_variable:
            raise ConfigurationError('No transmembrane potential found; check your configuration')
        return self.V_variable
    
    def find_membrane_capacitance(self):
        """Find and store the variable object representing the cell membrane capacitance.
        
        Uses first metadata, if present, then the configuration file."""
        self.Cm_variable = self._find_var('membrane_capacitance', self.Cm_definitions)

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
    def validate_metadata(self, assume_valid=False):
        """Check that the metadata annotations are 'sensible'.
        
        Ensures that only names we know are used, and that the same name isn't used for multiple
        variables.
        """
        vars = cellml_metadata.find_variables(self.doc.model, ('bqbiol:is', NSS['bqbiol']))
        self.metadata_vars = filter(lambda v: v.oxmeta_name, vars)
        if assume_valid:
            return
        names_used = [var.oxmeta_name for var in self.metadata_vars]
        DEBUG('metadata', 'Names found: ', names_used)
        # Check all metadata is allowed
        unknown_names = frozenset(names_used) - cellml_metadata.METADATA_NAMES
        if unknown_names:
            msg = ['Unrecognised oxmeta variable names found (run with --assume-valid to ignore):']
            msg.extend(sorted(unknown_names))
            raise ConfigurationError('\n  '.join(msg))
        # Check for duplicates
        d = {}
        for name in names_used:
            if name in d:
                raise ConfigurationError(name + ' metadata attribute is duplicated in the cellml file.')
            else:
                d[name] = name


######################################################################
#                    For running as an executable                    #
######################################################################

def get_options(args, default_options=None):
    """get_options(args):
    Process our command-line options.

    args is a list of options & positional arguments.
    
    default_options, if given, is an instance of optparse.Values created by a
    previous call to this function.
    """
    usage = 'usage: %prog [options] <cellml file or URI>'
    parser = optparse.OptionParser(version="%%prog %s" % __version__,
                                   usage=usage)
    # What type of translation is being performed
    parser.add_option('-T', '--translate',
                      dest='translate', action='store_true',
                      default=True,
                      help="output computer code [default]")
    parser.add_option('-C', '--output-cellml',
                      dest='translate', action='store_false',
                      help="output an annotated CellML file instead of translating, "
                      " on stdout unless -o specified")
    translators = sorted(CellMLTranslator.translators)
    parser.add_option('-t', '--translate-type',
                      type='choice', choices=translators,
                      default='Chaste', metavar='TYPE',
                      help="the type of code to output [default: %default].  "
                      "Choices: " + str(translators))
    parser.add_option('-o', dest='outfilename', metavar='OUTFILE',
                      help="write program code to OUTFILE [default action is"
                      " to use the input filename with a different extension]")
    # Global adjustment settings
    parser.add_option('--config-file',
                      action='append',
                      help="pathname of configuration file")
    parser.add_option('-A', '--fully-automatic',
                      action='store_true', default=False,
                      help="if human intervention is required, fail noisily")
    parser.add_option('--assume-valid',
                      action='store_true', default=False,
                      help="skip some of the model validation checks")
    parser.add_option('-V', '--transmembrane-potential',
                      default=None, metavar='POT_VAR',
                      help=
                      "POT_VAR is the full name of the variable representing"
                      " the transmembrane potential.  If not specified here,"
                      " the configuration file will be used.  Defaults to "
                      "'membrane,V'.")
    parser.add_option('-d', '--debug', action='store_true', default=False,
                      help="output debug info to stderr")
    parser.add_option('-D', '--debug-source', action='append',
                      help="only show debug info from the specified part of the"
                      " code.  This option may appear more than once to select"
                      " multiple sources.  Implies -d.")
    # What optimisations/transformations to do
    parser.add_option('-l', '--lookup-tables',
                      dest='lut', action='store_true', default=False,
                      help="perform a lookup table analysis")
    parser.add_option('-p', '--pe', '--partial-evaluation',
                      dest='pe', action='store_true', default=False,
                      help="partially evaluate the model")
    parser.add_option('-u', '--units-conversions',
                      action='store_true', default=False,
                      help="add explicit units conversion mathematics")
    parser.add_option('--Wu', '--warn-on-units-errors',
                      action='store_true', default=False,
                      dest='warn_on_units_errors',
                      help="give a warning instead of an error for"
                      " dimensional inconsistencies")
    parser.add_option('-j', '--maple-output',
                      metavar='FILENAME', default=None,
                      help="file containing output from a Maple script "
                      "generated using -J.  The generated code/CellML will "
                      "then contain a symbolic Jacobian as computed by Maple.")
    # Settings tweaking the generated code
    parser.add_option('-c', '--class-name', default=None,
                      help="explicitly set the name of the generated class")
    parser.add_option('-a', '--augment-class-name',
                      dest='augment_class_name', action='store_true',
                      default=False,
                      help="alter the class name to show what transformations"
                      " are used")
    parser.add_option('--no-timestamp',
                      action='store_true', default=False,
                      help="don't add a timestamp comment to generated files")
    parser.add_option('-J', '--do-jacobian-analysis',
                      action='store_true', default=False,
                      help="experimental Jacobian analysis; implies -t Maple")
    # Options specific to Maple output
    parser.add_option('--omit-constants',
                      action='store_true', default=False,
                      help="when generating Maple code, don't include "
                      "assignments of constants")
    parser.add_option('--compute-full-jacobian',
                      action='store_true', default=False,
                      help="make generated Maple code compute the full Jacobian"
                      " matrix, rather than just that for the nonlinear portion"
                      " of the ODE system")
    # Options specific to Chaste output
    parser.add_option('-y', '--dll', '--dynamically-loadable',
                      dest='dynamically_loadable',
                      action='store_true', default=False,
                      help="add code to allow the model to be compiled to a "
                      "shared library and dynamically loaded"
                      " (only works if -t Chaste is used)")
    parser.add_option('--use-chaste-stimulus',
                      action='store_true', default=False,
                      help="when generating Chaste code, use Chaste's stimulus"
                      " rather than that defined in the model")
    parser.add_option('-i', '--convert-interfaces',
                      action='store_true', default=False,
                      help="perform units conversions at interfaces to Chaste."
                      " (only works if -t Chaste is used)")
    parser.add_option('-m', '--use-modifiers',
                      action='store_true', default=False,
                      help="[experimental] add modifier functions for certain"
                      " metadata-annotated variables for use in sensitivity analysis"
                      " (only works if -t Chaste is used)")
    parser.add_option('--protocol',
                      help="[experimental] specify a simulation protocol to apply to"
                      " the model prior to translation")
    # Settings for lookup tables
    parser.add_option('--no-separate-lut-class', dest='separate_lut_class',
                      action='store_false', default=True,
                      help="don't put lookup tables in a separate class")
    parser.add_option('--row-lookup-method',
                      action='store_true', default=True,
                      help="add and use a method to look up a whole row of a table")
    parser.add_option('--no-row-lookup-method',
                      action='store_false', dest='row_lookup_method',
                      help="don't add and use a method to look up a whole row of a table")
    parser.add_option('--lt-index-uses-floor',
                      action='store_true', default=False,
                      help="use floor() to calculate LT indices, instead of "
                      "just casting")
    parser.add_option('--constrain-table-indices',
                      action='store_true', default=False,
                      help="constraint lookup table index variables to remain"
                      " within the bounds specified, rather than throwing an"
                      " exception if they go outside the bounds")
    parser.add_option('--bad-lt-layout-for-cache', dest='bad_tables_for_cache',
                      action='store_true', default=False,
                      help="[debug] use the old LT layout, with poorer cache"
                      " performance")
    # Settings for partial evaluation
    parser.add_option('--no-member-vars', dest='kept_vars_as_members',
                      action='store_false', default=True,
                      help="[debug] don't store kept variables as members")

    options, args = parser.parse_args(args, values=default_options)
    if len(args) != 1:
        parser.error("exactly one input CellML file must be specified")
    return options, args[0]


def load_model(model_file, options):
    """Load and validate a CellML model."""
    # Setup logging
    if options.debug_source:
        options.debug = True
    if options.debug:
        formatter = logging.Formatter(fmt="%(name)s: %(message)s")
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(formatter)
        handler.addFilter(OnlyDebugFilter())
        if options.debug_source:
            handler.addFilter(OnlyTheseSourcesFilter(options.debug_source))
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
    
    return doc

def run():
    """Translate the file given on the command line."""
    options, model_file = get_options(sys.argv[1:])
    doc = load_model(model_file, options)
    
    # Apply protocol, if given
    if options.protocol:
        import protocol
        protocol.apply_protocol_file(doc, options.protocol)

    config = ConfigurationStore(doc, options=options)
    if options.config_file:
        for config_file in options.config_file:
            config.read_configuration_file(config_file)
        config.finalize_config()
    else:
        # Use defaults
        config.find_transmembrane_potential()
        config.find_current_vars()

    config.validate_metadata(options.assume_valid)

    # These bits could do with improving, as they annotate more than is really needed!
    # We need to ensure PE doesn't remove ionic currents needed for GetIIonic
    config.annotate_currents_for_pe()
    # "Need" to ensure pe doesn't remove metadata-annotated variables (when using modifiers or default stimulus?)
    config.annotate_metadata_for_pe()

    class_name = options.class_name
    if options.augment_class_name and not class_name:
        class_name = u'CML_' + doc.model.name.replace('-', '_')
        if options.pe:
            class_name += '_pe'
        if options.lut:
            class_name += '_lut'
        if options.maple_output:
            class_name += '_be'
        if options.use_modifiers:
            class_name += '_sens'

    solver_info = SolverInfo(doc.model)

    output_filename = getattr(options, 'outfilename', None)
    if not options.translate and not output_filename:
        output_filename = 'stdout'

    if options.units_conversions:
        doc.model.add_units_conversions()

    if options.do_jacobian_analysis:
        lin = optimize.LinearityAnalyser()
        lin.analyse_for_jacobian(doc, V=config.V_variable)
        options.translate_type = 'Maple'
        options.maple_output = False # Just in case...!

    if options.maple_output:
        # Parse Jacobian matrix
        from maple_parser import MapleParser
        mp = MapleParser()
        jacobian_file = file(options.maple_output) # TODO: Error checking
        doc.model._cml_jacobian = mp.parse(jacobian_file)
        jacobian_file.close()
        # Rearrange linear ODEs
        lin = optimize.LinearityAnalyser()
        lin.analyse_for_jacobian(doc, V=config.V_variable)
        lin.rearrange_linear_odes(doc)
        # Add info as XML
        solver_info.add_all_info()
        # Analyse the XML, adding cellml_variable references, etc.
        solver_info.add_variable_links()

    if options.pe:
        # Do partial evaluation
        pe = optimize.PartialEvaluator()
        pe.parteval(doc, solver_info)

    if options.lut:
        # Do the lookup table analysis
        lut = optimize.LookupTableAnalyser()
        if options.config_file:
            config.find_lookup_variables(options.pe)
        elif options.pe:
            lut.set_params(table_var=u'membrane__V')
        lut.analyse_model(doc, solver_info)

    if options.translate:
        # Translate to code
        klasses = CellMLTranslator.translators
        klass = klasses[options.translate_type]
        initargs = {'add_timestamp': not options.no_timestamp,
                    'options': options}
        transargs = {'v_variable': config.V_variable}
        transargs['row_lookup_method'] = options.row_lookup_method
        transargs['lt_index_uses_floor'] = options.lt_index_uses_floor
        transargs['constrain_table_indices'] = options.constrain_table_indices
        transargs['bad_tables_for_cache'] = options.bad_tables_for_cache
        if issubclass(klass, CellMLToMapleTranslator):
            initargs['omit_constants'] = options.omit_constants
            initargs['compute_full_jacobian'] = options.compute_full_jacobian
        elif issubclass(klass, CellMLToChasteTranslator):
            solver_info.add_membrane_ionic_current()
            transargs['use_chaste_stimulus'] = options.use_chaste_stimulus
            transargs['separate_lut_class'] = options.separate_lut_class
            transargs['convert_interfaces'] = options.convert_interfaces
            transargs['kept_vars_as_members'] = options.kept_vars_as_members
            transargs['use_modifiers'] = options.use_modifiers
            transargs['dynamically_loadable'] = options.dynamically_loadable
            transargs['use_protocol'] = bool(options.protocol)
        t = klass(**initargs)
        t.translate(doc, model_file, output_filename, class_name=class_name, **transargs)
        cellml_metadata.remove_model(doc.model)
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
    def euler(doc, t, nsteps=1000, dt=0.01):
        global tvar, state_vars, exprs
        tvar = t.free_vars[0]
        state_vars = t.state_vars
        for var in state_vars:
            var.set_value(float(var.initial_value))
        tvar.set_value(0.0)
        exprs = [e for e in doc.model.get_assignments()
                 if isinstance(e, mathml_apply)]
        for _ in range(nsteps):
            for expr in exprs:
                expr.evaluate()
            tvar.set_value(tvar.get_value() + dt)
            for var in state_vars:
                var.set_value(var.get_value() +
                              dt * var.get_value(ode=tvar))
        return

    def writefile(doc, outfn='test.cml'):
        # Write out CellML file
        st = open_output_stream(outfn)
        doc.xml(indent=1, stream=st)
        st.close()
        return

    def show_usage(doc):
        for comp in doc.model.component:
            for var in comp.variable:
                print var.fullname(), var._cml_usage_count


    def fix_divide_by_zero(doc):
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
