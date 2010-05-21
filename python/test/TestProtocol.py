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

import os
import sys
import unittest

# Get PyCml modules
sys.path[0:0] = ['python/pycml']
import protocol
import translate
import pycml

class TestProtocol(unittest.TestCase):
    """Tests of the Protocol system."""
    def LoadModel(self, model_filename, options=[]):
        args = ['-C', '-A', '--assume-valid', model_filename] + options
        options, model_file = translate.get_options(args)
        doc = translate.load_model(model_file, options, pycml_path='python/pycml')
        return doc
    
    def TestChangeInitialValue(self):
        doc = self.LoadModel('heart/src/odes/cellml/luo_rudy_1991.cellml')
        p = protocol.Protocol(doc.model, multi_stage=True)
        # Change initial value for V
        v_init = u'-80'
        new_V = pycml.cellml_variable.create_new(doc.model, u'membrane,V', u'millivolt',
                                                 initial_value=v_init,
                                                 interfaces={u'public': u'out'})
        p.inputs.append(new_V)
        # Check we are actually changing things
        V = doc.model.get_variable_by_name('membrane', 'V')
        self.failIfEqual(V.initial_value, v_init)
        time = doc.model.get_variable_by_name(u'environment', u'time')
        dv_dt = V._get_ode_dependency(time)
        # Apply protocol to model
        p.modify_model()
        # Check the changes
        V = doc.model.get_variable_by_name('membrane', 'V')
        self.assert_(V is new_V)
        self.assertEqual(V.initial_value, v_init)
        self.assertEqual(V.get_type(), pycml.VarTypes.State)
        self.assertEqual(V._get_dependencies(), [])
        self.assert_(V._get_ode_dependency(time) is dv_dt)
