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

import validator
    
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
        
    def _reanalyse_model(self, error_handler):
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


class InterfaceGenerator(ModelModifier):
    """Class for generating an interface between a CellML model and external code.
    
    This contains functionality for users to describe the interface desired by the external code, i.e.
    which variables are inputs and/or outputs, and expected units.  It will then create a new component
    within the CellML model containing these variables, and add units conversions where required.  The
    external code then only needs to interact with this new component.
    """
    pass
