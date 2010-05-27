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
        self._doc = doc
        return doc
    
    def CreateLr91Test(self):
        doc = self.LoadModel('heart/src/odes/cellml/luo_rudy_1991.cellml')
        p = protocol.Protocol(doc.model, multi_stage=True)
        return p
    
    def NewVariable(self, name, units, id=None, initial_value=None, interfaces={}):
        return pycml.cellml_variable.create_new(self._doc, name, units,
                                                id=id, initial_value=initial_value,
                                                interfaces=interfaces)
    
    def NewApply(self, operator, operands):
        return pycml.mathml_apply.create_new(self._doc, operator, operands)
    
    def NewAssign(self, lhs, rhs):
        return self.NewApply(u'eq', [lhs, rhs])
    
    def TestChangeInitialValue(self):
        doc = self.LoadModel('heart/src/odes/cellml/luo_rudy_1991.cellml')
        p = protocol.Protocol(doc.model, multi_stage=True)
        # Change initial value for V
        v_init = u'-80'
        new_V = self.NewVariable(u'membrane,V', u'millivolt',
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
        # An error case: changing interfaces is not allowed
        new_V = self.NewVariable(u'membrane,V', u'millivolt',
                                 initial_value=v_init,
                                 interfaces={u'public': u'none'})
        p.inputs = [new_V]
        self.assertRaises(AssertionError, p.modify_model)
        new_V = self.NewVariable(u'membrane,V', u'millivolt',
                                 initial_value=v_init)
        p.inputs = [new_V]
        self.assertRaises(AssertionError, p.modify_model)
        
    def TestChangeMathsErrors(self):
        """Test that changing mathematics reports errors for wrong inputs.
        
        In particular, it should throw if we try to define a Mapped variable
        using mathematics, or if we reference non-existent variables.
        """
        doc = self.LoadModel('heart/src/odes/cellml/luo_rudy_1991.cellml')
        p = protocol.Protocol(doc.model, multi_stage=True)
        # T is a Mapped variable in time_dependent_potassium_current
        T = doc.model.get_variable_by_name(u'time_dependent_potassium_current', u'T')
        T_const = self.NewAssign(u'time_dependent_potassium_current,T',
                                 (unicode(T.get_value()), T.units))
        p.inputs.append(T_const)
        self.assertRaises(protocol.ProtocolError, p.modify_model)
        # Same should apply if we replace by an ODE
        T_ode = pycml.mathml_diff.create_new(T, u'time_dependent_potassium_current,time',
                                             u'time_dependent_potassium_current,T',
                                             (u'1', u'dimensionless'))
        p.inputs = [T_ode]
        self.assertRaises(protocol.ProtocolError, p.modify_model)
        # Referencing a missing variable is also an error
        T = doc.model.get_variable_by_name(u'membrane', u'T')
        T_const = self.NewAssign(u'membrane,T', u'membrane,missing')
        p.inputs = [T_const]
        self.assertRaises(KeyError, p.modify_model)
        # As is assigning to a missing variable
        missing = self.NewAssign(u'membrane,missing', u'membrane,T')
        p.inputs = [missing]
        self.assertRaises(KeyError, p.modify_model)

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
        Cai_const = self.NewAssign(u'intracellular_calcium_concentration,Cai',
                                   (u'0.0002', u'millimolar'))
        p.inputs.append(Cai_const)
        self.failUnlessEqual(Cai.get_type(), pycml.VarTypes.State)
        # Change computed defintion for g_K (computed->computed)
        g_K = doc.model.get_variable_by_name(u'time_dependent_potassium_current', u'g_K')
        g_K_const = self.NewAssign(u'time_dependent_potassium_current,g_K',
                                   (u'0.282', u'milliS_per_cm2'))
        p.inputs.append(g_K_const)
        old_g_K = g_K._get_dependencies()[0]
        self.failUnlessEqual(g_K.get_type(), pycml.VarTypes.Computed)
        # Change PR_NaK from constant->computed (with same value)
        PR_NaK = doc.model.get_variable_by_name(u'time_dependent_potassium_current', u'PR_NaK')
        PR_NaK_const = self.NewAssign(u'time_dependent_potassium_current,PR_NaK',
                                      (PR_NaK.initial_value, PR_NaK.units))
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
        rhs = self.NewApply(u'plus', currents)
        V_new = pycml.mathml_diff.create_new(V, u'environment,time', u'membrane,V', rhs)
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
        doc = self.LoadModel('heart/src/odes/cellml/luo_rudy_1991.cellml')
        p = protocol.Protocol(doc.model, multi_stage=True)
        # Add a new 'I_total' variable and computation
        I_total = self.NewVariable(u'membrane,I_total', u'microA_per_cm2')
        p.inputs.append(I_total)
        currents = map(lambda i: u'membrane,i_' + i, ['Na', 'si', 'b'])
        currents.append(u'membrane,I_stim')
        currents.append(u'membrane,potassium_currents')
        rhs = self.NewApply(u'plus', currents)
        defn = self.NewAssign(u'membrane,I_total', rhs)
        p.inputs.append(defn)
        # Add a variable to the protocol component
        one = self.NewVariable(u'one', u'dimensionless', initial_value=u'1')
        p.inputs.append(one)
        self.assertRaises(KeyError, doc.model.get_component_by_name, u'protocol')
        # Apply protocol to model
        p.modify_model()
        # Check the changes
        self.assertEqual(I_total, doc.model.get_variable_by_name(u'membrane', u'I_total'))
        self.assertEqual(I_total.component, doc.model.get_component_by_name(u'membrane'))
        self.assertEqual(I_total.component, I_total.xml_parent)
        self.assertEqual(I_total.get_type(), pycml.VarTypes.Computed)
        self.assertEqual(I_total._get_dependencies()[0], defn)
        self.assertEqual(I_total.component, defn.component)
        self.assertEqual(str(defn), 'I_totali_Nai_sii_bI_stimpotassium_currents')
        self.assertEqual(one, doc.model.get_variable_by_name(u'protocol', u'one'))
        self.assertEqual(one.component, doc.model.get_component_by_name(u'protocol'))
        self.assertEqual(one.component.name, u'protocol')
        self.assertEqual(one.model, doc.model)
        self.assertEqual(one.get_type(), pycml.VarTypes.Constant)
        # Now apply a new protocol, to create a protocol_ component
        # Also checks connection creation
        p2 = protocol.Protocol(doc.model, multi_stage=True)
        two = self.NewVariable(u'two', u'dimensionless')
        two_defn = self.NewAssign(u'two', u'membrane,time')
        p2.inputs.append(two)
        p2.inputs.append(two_defn)
        self.assertRaises(KeyError, doc.model.get_component_by_name, u'protocol_')
        p2.modify_model()
        self.assertEqual(two, doc.model.get_variable_by_name(u'protocol_', u'two'))
        self.assertEqual(two.component, doc.model.get_component_by_name(u'protocol_'))
        self.assertEqual(two.component.name, u'protocol_')
        self.assertEqual(two.model, doc.model)
        self.assertNotEqual(one.component, two.component)

    def TestConnectionMaking(self):
        """Test that creating new connections between variables works.
        
        This is done a little in other tests; here we try to exercise all cases.
        """
        p = self.CreateLr91Test()
        # Add a new variable in an encapsulated component
        src = self.NewVariable(u'time_dependent_potassium_current_Xi_gate,new_src', u'dimensionless')
        p.inputs.append(src)
        # Add mathematics using it in another encapsulated component
        targ1 = self.NewVariable(u'time_independent_potassium_current_K1_gate,target1', u'dimensionless')
        use1 = self.NewAssign(targ1.name, u'time_dependent_potassium_current_Xi_gate,new_src')
        p.inputs.extend([targ1, use1])
        # Add mathematics using it in a top-level component
        targ2 = self.NewVariable(u'target2', u'dimensionless')
        use2 = self.NewAssign(targ2.name, u'time_dependent_potassium_current_Xi_gate,new_src')
        p.inputs.extend([use2, targ2])
        # Add 2 new mapped variables with chained dependency
        stim1 = self.NewVariable(u'ionic_concentrations,stim1', u'microA_per_cm2')
        stim2 = self.NewVariable(u'intracellular_calcium_concentration,stim2', u'microA_per_cm2')
        p.inputs.extend([stim1, stim2])
        p.inputs.append(self.NewAssign(stim1.name, u'membrane,I_stim'))
        p.inputs.append(self.NewAssign(stim2.name, stim1.name))
        # Apply protocol and test
        p.modify_model()
        model = self._doc.model
        self.assertEqual(src.model, model)
        src1 = model.get_variable_by_name(u'time_independent_potassium_current_K1_gate', u'new_src')
        src2 = model.get_variable_by_name(u'protocol', u'new_src')
        self.assertNotEqual(src, src1)
        self.assertEqual(src1.get_source_variable(recurse=True), src)
        self.assertNotEqual(src1.get_source_variable(), src)
        self.assertEqual(src2.get_source_variable(recurse=True), src)
        self.assertNotEqual(src2.get_source_variable(), src)
        self.assertEqual(targ1._get_dependencies()[0], use1)
        self.assertEqual(use1._get_dependencies()[0], src1)
        self.assertEqual(targ2._get_dependencies()[0], use2)
        self.assertEqual(use2._get_dependencies()[0], src2)
        self.assertEqual(src1.get_type(), pycml.VarTypes.Mapped)
        self.assertEqual(src2.get_type(), pycml.VarTypes.Mapped)
        self.assertEqual(targ1.get_type(), pycml.VarTypes.Computed)
        self.assertEqual(targ2.get_type(), pycml.VarTypes.Computed)
        
        # Also check error cases:
        # Where we try to create a connection but there's an existing conflicting variable
        p = self.CreateLr91Test()
        time2 = self.NewVariable(u'environment,time2', u'millisecond')
        t2def = self.NewAssign(time2.name, u'membrane,time')
        p.inputs = [time2, t2def]
        self.assertRaises(protocol.ProtocolError, p.modify_model)
        
        # The following are odd cases that do actually work!
        # Where we use the real source variable as target
        p = self.CreateLr91Test()
        time2 = self.NewVariable(u'fast_sodium_current,time2', u'millisecond')
        t2def = self.NewAssign(time2.name, u'fast_sodium_current_m_gate,time')
        p.inputs = [time2, t2def]
        p.modify_model()
        # The order in which we add assignments shouldn't matter
        p = self.CreateLr91Test()
        n1 = u'fast_sodium_current,stim1'
        n2 = u'slow_inward_current,stim2'
        p.inputs = [self.NewVariable(n1, u'microA_per_cm2'),
                    self.NewVariable(n2, u'microA_per_cm2'),
                    self.NewAssign(n2, n1),
                    self.NewAssign(n1, u'membrane,I_stim')]
        p.modify_model()
        # Likewise for the order in which they appear in the model
        p = self.CreateLr91Test()
        n2 = u'fast_sodium_current,stim1'
        n1 = u'slow_inward_current,stim2'
        p.inputs = [self.NewVariable(n1, u'microA_per_cm2'),
                    self.NewVariable(n2, u'microA_per_cm2'),
                    self.NewAssign(n2, n1),
                    self.NewAssign(n1, u'membrane,I_stim')]
        p.modify_model()
