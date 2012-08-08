/*

Copyright (C) University of Oxford, 2005-2009

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

*/

#ifndef _TESTMETADATACELLMLMODELS_HPP_
#define _TESTMETADATACELLMLMODELS_HPP_

#include <cxxtest/TestSuite.h>
#include "Exception.hpp"

#include <boost/shared_ptr.hpp>
#include "ZeroStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "AbstractCardiacCell.hpp"
#include "AbstractCvodeCell.hpp"
#include "AbstractCardiacCellWithModifiers.hpp"
#include "Modifiers.hpp"

#include "hund_rudy_2004_a.hpp"
#include "hund_rudy_2004_aOpt.hpp"
#include "hund_rudy_2004_aCvode.hpp"
#include "hund_rudy_2004_aCvodeOpt.hpp"
#include "mahajan_2008.hpp"
#include "mahajan_2008Opt.hpp"
#include "mahajan_2008Cvode.hpp"
#include "mahajan_2008CvodeOpt.hpp"
#include "Shannon2004.hpp"
//#include "Shannon2004Opt.hpp" // Opt files not made to save time as Shannon is in trunk.
#include "Shannon2004Cvode.hpp"
//#include "Shannon2004CvodeOpt.hpp" // Opt files not made to save time as Shannon is in trunk.
#include "ten_tusscher_model_2006_epi.hpp"
#include "ten_tusscher_model_2006_epiOpt.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "ten_tusscher_model_2006_epiCvodeOpt.hpp"
#include "grandi2010ss.hpp"
#include "grandi2010ssOpt.hpp"
#include "grandi2010ssCvode.hpp"
#include "grandi2010ssCvodeOpt.hpp"

class TestMetadataCellmlModels : public CxxTest::TestSuite
{
public:

    /**
     * This test just ensures that all cell models can be compiled and given metadata values.
     *
     * (i.e. that pyCML and the metadata tags generated valid Chaste cell models.)
     */
    void TestMetadataHasBeenCorrectlyTranslatedToChaste(void) throw(Exception)
    {
        // Regular stimulus and an ode solver object
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus);

        // 'Modifier' classes allow us to intervene in the right hand side of the ODE.
        boost::shared_ptr<AbstractModifier> a = boost::shared_ptr<AbstractModifier>(new FactorModifier(0.9));
        boost::shared_ptr<AbstractModifier> b = boost::shared_ptr<AbstractModifier>(new FactorModifier(0.8));
        boost::shared_ptr<AbstractModifier> c = boost::shared_ptr<AbstractModifier>(new FactorModifier(0.7));

        for (unsigned model_index = 1u; model_index < 6u ; model_index++)
        {
            AbstractCardiacCellWithModifiers<AbstractCardiacCell>* p_chaste_cell;
            AbstractCardiacCellWithModifiers<AbstractCardiacCell>* p_chaste_cell_opt;
            AbstractCardiacCellWithModifiers<AbstractCvodeCell>* p_cvode_cell;
            AbstractCardiacCellWithModifiers<AbstractCvodeCell>* p_cvode_cell_opt;

            std::cout << "model index = " << model_index << std::endl << std::flush;
            switch(model_index)
            {
                case 1u:
                {
                    p_chaste_cell = new CellShannon2004FromCellML(p_solver, p_stimulus);
                    p_chaste_cell_opt = new CellShannon2004FromCellML(p_solver, p_stimulus); // Opt files not made as Shannon is in trunk
                    p_cvode_cell = new CellShannon2004FromCellMLCvode(p_solver, p_stimulus);
                    p_cvode_cell_opt = new CellShannon2004FromCellMLCvode(p_solver, p_stimulus);  // Opt files not made as Shannon is in trunk
                    break;
                }
                case 2u:
                {
                    p_chaste_cell = new Cellten_tusscher_model_2006_epiFromCellML(p_solver, p_stimulus);
                    p_chaste_cell_opt = new Cellten_tusscher_model_2006_epiFromCellMLOpt(p_solver, p_stimulus);
                    p_cvode_cell = new Cellten_tusscher_model_2006_epiFromCellMLCvode(p_solver, p_stimulus);
                    p_cvode_cell_opt = new Cellten_tusscher_model_2006_epiFromCellMLCvodeOpt(p_solver, p_stimulus);
                    break;
                }
                case 3u:
                {
                    p_chaste_cell = new Cellmahajan_2008FromCellML(p_solver, p_stimulus);
                    p_chaste_cell_opt = new Cellmahajan_2008FromCellMLOpt(p_solver, p_stimulus);
                    p_cvode_cell = new Cellmahajan_2008FromCellMLCvode(p_solver, p_stimulus);
                    p_cvode_cell_opt = new Cellmahajan_2008FromCellMLCvodeOpt(p_solver, p_stimulus);
                    break;
                }
                case 4u:
                {
                    p_chaste_cell = new Cellhund_rudy_2004_aFromCellML(p_solver, p_stimulus);
                    p_chaste_cell_opt = new Cellhund_rudy_2004_aFromCellMLOpt(p_solver, p_stimulus);
                    p_cvode_cell = new Cellhund_rudy_2004_aFromCellMLCvode(p_solver, p_stimulus);
                    p_cvode_cell_opt = new Cellhund_rudy_2004_aFromCellMLCvodeOpt(p_solver, p_stimulus);
                    break;
                }
                case 5u:
                {
                    p_chaste_cell = new Cellgrandi2010ssFromCellML(p_solver, p_stimulus);
                    p_chaste_cell_opt = new Cellgrandi2010ssFromCellMLOpt(p_solver, p_stimulus);
                    p_cvode_cell = new Cellgrandi2010ssFromCellMLCvode(p_solver, p_stimulus);
                    p_cvode_cell_opt = new Cellgrandi2010ssFromCellMLCvodeOpt(p_solver, p_stimulus);
                    break;
                }
                default:
                {
                    EXCEPTION("Unrecognised cell model");
                }
            }
            p_chaste_cell->SetModifier("membrane_fast_sodium_current_conductance", a);
            p_chaste_cell->SetModifier("membrane_L_type_calcium_current_conductance", b);
            p_chaste_cell->SetModifier("membrane_rapid_delayed_rectifier_potassium_current_conductance", c);
            TS_ASSERT_EQUALS(p_chaste_cell->HasCellMLDefaultStimulus(), true);
            p_chaste_cell->UseCellMLDefaultStimulus();

            p_chaste_cell_opt->SetModifier("membrane_fast_sodium_current_conductance", a);
            p_chaste_cell_opt->SetModifier("membrane_L_type_calcium_current_conductance", b);
            p_chaste_cell_opt->SetModifier("membrane_rapid_delayed_rectifier_potassium_current_conductance", c);
            TS_ASSERT_EQUALS(p_chaste_cell_opt->HasCellMLDefaultStimulus(), true);
            p_chaste_cell_opt->UseCellMLDefaultStimulus();

            p_cvode_cell->SetModifier("membrane_fast_sodium_current_conductance", a);
            p_cvode_cell->SetModifier("membrane_L_type_calcium_current_conductance", b);
            p_cvode_cell->SetModifier("membrane_rapid_delayed_rectifier_potassium_current_conductance", c);
            TS_ASSERT_EQUALS(p_cvode_cell->HasCellMLDefaultStimulus(), true);
            p_cvode_cell->UseCellMLDefaultStimulus();

            p_cvode_cell_opt->SetModifier("membrane_fast_sodium_current_conductance",a);
            p_cvode_cell_opt->SetModifier("membrane_L_type_calcium_current_conductance", b);
            p_cvode_cell_opt->SetModifier("membrane_rapid_delayed_rectifier_potassium_current_conductance", c);
            TS_ASSERT_EQUALS(p_cvode_cell_opt->HasCellMLDefaultStimulus(), true);
            p_cvode_cell_opt->UseCellMLDefaultStimulus();

            delete p_chaste_cell;
            delete p_chaste_cell_opt;
            delete p_cvode_cell;
            delete p_cvode_cell_opt;
        }
    }

};


#endif //_TESTMETADATACELLMLMODELS_HPP_
