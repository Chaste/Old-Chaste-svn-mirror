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
        
    def TestChangeMathsErrors(self):
        """Test that changing mathematics reports errors for wrong inputs.
        
        In particular, it should throw if we try to define a Mapped variable
        using mathematics. 
        """
        doc = self.LoadModel('heart/src/odes/cellml/luo_rudy_1991.cellml')
        p = protocol.Protocol(doc.model, multi_stage=True)
        # T is a Mapped variable in time_dependent_potassium_current
        T = doc.model.get_variable_by_name(u'time_dependent_potassium_current', u'T')
        T_const = pycml.mathml_apply.create_new(T, u'eq', [u'time_dependent_potassium_current,T',
                                                           (unicode(T.get_value()), T.units)])
        p.inputs.append(T_const)
        self.assertRaises(protocol.ProtocolError, p.modify_model)
        # Same should apply if we replace by an ODE
        T_ode = pycml.mathml_diff.create_new(T, u'time_dependent_potassium_current,time',
                                             u'time_dependent_potassium_current,T',
                                             (u'1', u'dimensionless'))
        p.inputs = [T_ode]
        self.assertRaises(protocol.ProtocolError, p.modify_model)

    def TestChangeMathsToComputed(self):
        """Test changing the mathematical definition of a variable to type Computed.
        
        Covers the cases of variables originally of type State, Computed, and Constant.
        """
        doc = self.LoadModel('heart/src/odes/cellml/luo_rudy_1991.cellml')
        p = protocol.Protocol(doc.model, multi_stage=True)
        # Add maths setting Cai to a constant, thus replacing the ODE (state->computed)
        Cai = doc.model.get_variable_by_name(u'intracellular_calcium_concentration', u'Cai')
        time = doc.model.get_variable_by_name(u'environment', u'time')
        ode = Cai._get_ode_dependency(time)
        Cai_const = pycml.mathml_apply.create_new(ode, u'eq', [u'intracellular_calcium_concentration,Cai',
                                                               (u'0.0002', u'millimolar')])
        p.inputs.append(Cai_const)
        self.failUnlessEqual(Cai.get_type(), pycml.VarTypes.State)
        # Change computed defintion for g_K (computed->computed)
        g_K = doc.model.get_variable_by_name(u'time_dependent_potassium_current', u'g_K')
        g_K_const = pycml.mathml_apply.create_new(g_K, u'eq', [u'time_dependent_potassium_current,g_K',
                                                               (u'0.282', u'milliS_per_cm2')])
        p.inputs.append(g_K_const)
        old_g_K = g_K._get_dependencies()[0]
        self.failUnlessEqual(g_K.get_type(), pycml.VarTypes.Computed)
        # Change PR_NaK from constant->computed (with same value)
        PR_NaK = doc.model.get_variable_by_name(u'time_dependent_potassium_current', u'PR_NaK')
        PR_NaK_const = pycml.mathml_apply.create_new(PR_NaK, u'eq',
                                                     [u'time_dependent_potassium_current,PR_NaK',
                                                      (PR_NaK.initial_value, PR_NaK.units)])
        p.inputs.append(PR_NaK_const)
        self.failUnlessEqual(PR_NaK.get_type(), pycml.VarTypes.Constant)
        # Apply protocol to model
        p.modify_model()
        # Check the changes
        self.assertEqual(Cai, doc.model.get_variable_by_name(u'intracellular_calcium_concentration', u'Cai'))
        self.assertEqual(ode.xml_parent, None)
        self.assertRaises(translate.MathsError, Cai._get_ode_dependency, time)
        self.assertEqual(Cai._get_dependencies()[0], Cai_const)
        self.assertEqual(Cai.get_type(), pycml.VarTypes.Computed)
        self.assertEqual(old_g_K.xml_parent, None)
        self.assertEqual(g_K._get_dependencies()[0], g_K_const)
        self.assertEqual(g_K.get_type(), pycml.VarTypes.Computed)
        self.assertEqual(PR_NaK, doc.model.get_variable_by_name(u'time_dependent_potassium_current', u'PR_NaK'))
        self.failIf(hasattr(PR_NaK, u'initial_value'))
        self.assertEqual(PR_NaK._get_dependencies()[0], PR_NaK_const)
        self.assertEqual(PR_NaK.get_type(), pycml.VarTypes.Computed)

    def TestChangeMathsToState(self):
        """Test changing the mathematical definition of a variable to type State.
        
        Covers the cases of variables originally of type State, Computed, and Constant.
        """
        doc = self.LoadModel('heart/src/odes/cellml/luo_rudy_1991.cellml')
        p = protocol.Protocol(doc.model, multi_stage=True)
        # state->state: membrane,V
        V = doc.model.get_variable_by_name(u'membrane', u'V')
        self.failUnlessEqual(V.get_type(), pycml.VarTypes.State)
        time = doc.model.get_variable_by_name(u'membrane', u'time')
        V_old = V._get_ode_dependency(time)
        currents = map(lambda i: u'membrane,i_' + i, ['Na', 'si', 'K', 'b'])
        rhs = pycml.mathml_apply.create_new(V, u'plus', currents)
        V_new = pycml.mathml_diff.create_new(V, u'membrane,time', u'membrane,V', rhs)
        V_new_str = str(V_new)
        p.inputs.append(V_new)
        # constant->state: membrane,T
        T = doc.model.get_variable_by_name(u'membrane', u'T')
        self.assertEqual(T.get_type(), pycml.VarTypes.Constant)
        self.assertEqual(T._get_dependencies(), [])
        T_new = pycml.mathml_diff.create_new(T, u'membrane,time', u'membrane,T',
                                             (u'1', u'dimensionless')) # TODO: proper units
        p.inputs.append(T_new)
        # computed->state: membrane,FonRT
        FonRT = doc.model.get_variable_by_name(u'membrane', u'FonRT')
        self.assertEqual(FonRT.get_type(), pycml.VarTypes.Computed)
        FonRT_old = FonRT._get_dependencies()[0]
        FonRT_new = pycml.mathml_diff.create_new(FonRT, u'membrane,time', u'membrane,FonRT',
                                                 (u'0', u'per_millivolt_millisecond'))
        # TODO: FonRT now has no initial condition; this is a validation warning so OK...
        p.inputs.append(FonRT_new)
        # Apply protocol to model
        p.modify_model()
        # Check the changes
        self.assertEqual(V, doc.model.get_variable_by_name(u'membrane', u'V'))
        self.assertEqual(V_old.xml_parent, None)
        self.assertNotEqual(V_old, V._get_ode_dependency(time))
        self.assertEqual(V_new, V._get_ode_dependency(time))
        self.assertEqual(V.get_type(), pycml.VarTypes.State)
        self.assertNotEqual(str(V_new), V_new_str) # Variable refs get renamed
        self.assertEqual(str(V_new), 'timeVi_Nai_sii_Ki_b')
        self.assertEqual(T.get_type(), pycml.VarTypes.State)
        self.assertEqual(T_new, T._get_ode_dependency(time))
        self.assertEqual(FonRT.id, u'FonRT')
        self.assertEqual(FonRT.get_type(), pycml.VarTypes.State)
        self.assertEqual(FonRT_new, FonRT._get_ode_dependency(time))
        self.assertEqual(FonRT._get_dependencies(), [])
        self.assertEqual(FonRT_old.xml_parent, None)

    def TestAddNewVariables(self):
        """Test we can add new variables along with defining mathematics.
        """
        pass