/*

Copyright (C) University of Oxford, 2005-2010

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

#ifndef TESTSTIMULUS_HPP_
#define TESTSTIMULUS_HPP_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <cfloat>
#include <cmath>

#include "AbstractStimulusFunction.hpp"
#include "SimpleStimulus.hpp"
#include "RegularStimulus.hpp"
#include "TimeStepper.hpp"
#include "ZeroStimulus.hpp"
#include "MultiStimulus.hpp"
#include "OutputFileHandler.hpp"

class TestStimulus : public CxxTest::TestSuite
{
public:

    void TestSimpleStimulus()
    {
        double magnitude_of_stimulus = 1.0;
        double duration_of_stimulus  = 0.5;  // ms
        double when = 100.0;

        SimpleStimulus initial_at_zero(magnitude_of_stimulus, duration_of_stimulus);

        TS_ASSERT_EQUALS(initial_at_zero.GetStimulus(0.0), magnitude_of_stimulus);
        TS_ASSERT_EQUALS(initial_at_zero.GetStimulus(duration_of_stimulus*(1-DBL_EPSILON)),
            magnitude_of_stimulus);
        TS_ASSERT_EQUALS(initial_at_zero.GetStimulus(duration_of_stimulus), magnitude_of_stimulus);

        //Made more sloppy
        TS_ASSERT_EQUALS(initial_at_zero.GetStimulus(duration_of_stimulus*(1+3*DBL_EPSILON)),
            0.0);

        SimpleStimulus initial_later(magnitude_of_stimulus, duration_of_stimulus, when);

        TS_ASSERT_EQUALS(initial_later.GetStimulus(when*(1-DBL_EPSILON)),
            0.0);
        TS_ASSERT_EQUALS(initial_later.GetStimulus(when), magnitude_of_stimulus);
        TS_ASSERT_EQUALS(initial_later.GetStimulus(when*(1+DBL_EPSILON)),
            magnitude_of_stimulus);
        TS_ASSERT_EQUALS(initial_later.GetStimulus((when+duration_of_stimulus)*(1-DBL_EPSILON)),
            magnitude_of_stimulus);
        TS_ASSERT_EQUALS(initial_later.GetStimulus(when+duration_of_stimulus), magnitude_of_stimulus);
        TS_ASSERT_EQUALS(initial_later.GetStimulus((when+duration_of_stimulus)*(1+2*DBL_EPSILON)),
            0.0);
    }


    void TestDelayedSimpleStimulus()
    {
        double magnitude_of_stimulus = 1.0;
        double duration_of_stimulus  = 0.002;  // ms
        double when = 0.06;

        SimpleStimulus initial(magnitude_of_stimulus, duration_of_stimulus, when);
        TS_ASSERT_EQUALS(initial.GetStimulus(0.062),  magnitude_of_stimulus);
    }

    void TestRegularStimulus()
    {
        double magnitude_of_stimulus = 1.0;
        double duration_of_stimulus = 0.5;  // ms
        double period = 1000.0; // 1s
        double when = 100.0;

        RegularStimulus regular_stimulus(magnitude_of_stimulus,
                                         duration_of_stimulus,
                                         period,
                                         when);

        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(0.0),
            0.0);

        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.0*(1-DBL_EPSILON)),
            0.0);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.0),
            magnitude_of_stimulus);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.0*(1+DBL_EPSILON)),
            magnitude_of_stimulus);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.5*(1-DBL_EPSILON)),
            magnitude_of_stimulus);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.5),
            magnitude_of_stimulus);

        TS_ASSERT_DELTA(regular_stimulus.GetStimulusPeriod(),period,1e-9);

        //Made more sloppy
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.5+(1000*DBL_EPSILON)),
            0.0);

        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(5100.0*(1-DBL_EPSILON)),
            0.0);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(5100.0),
            magnitude_of_stimulus);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(5100.0*(1+DBL_EPSILON)),
            magnitude_of_stimulus);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(5100.5*(1-DBL_EPSILON)),
            magnitude_of_stimulus);

        //An overeager floating point optimiser may fail the following test
        //by turning the stimulus off too early
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(5100.5),
                         magnitude_of_stimulus);

        //Made more sloppy
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(5100.5*(1+2*DBL_EPSILON)),
            0.0);

        // Test SetStimulusStartTime method, by shifting the stimulus slightly
        double delay = 10.0;
        regular_stimulus.SetStimulusStartTime(when + delay);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(5100.5 + delay),
                         magnitude_of_stimulus);
    }

    void TestRegularStimulusStopping()
    {
        double magnitude_of_stimulus = 1.0;
        double duration_of_stimulus  = 0.5;  // ms
        double period = 1000.0; // 1s
        double when = 100.0;
        double until = 3500.0;

        RegularStimulus regular_stimulus(magnitude_of_stimulus,
                                         duration_of_stimulus,
                                         period,
                                         when,
                                         until);

        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(0.0),
            0.0);

        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.0*(1-DBL_EPSILON)),
            0.0);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.0),
            magnitude_of_stimulus);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.0*(1+DBL_EPSILON)),
            magnitude_of_stimulus);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.5*(1-DBL_EPSILON)),
            magnitude_of_stimulus);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.5),
            magnitude_of_stimulus);

        //Made more sloppy
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(3100.0*(1-DBL_EPSILON)),
            0.0);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(3100.0),
            magnitude_of_stimulus);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(3100.0*(1+DBL_EPSILON)),
            magnitude_of_stimulus);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(3100.5*(1-DBL_EPSILON)),
            magnitude_of_stimulus);

        //Made more sloppy
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(100.5+(1000*DBL_EPSILON)),
            0.0);

        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(4100.0*(1-DBL_EPSILON)),
            0.0);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(4100.0),
            0.0);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(4100.0*(1+DBL_EPSILON)),
            0.0);
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(4100.5*(1-DBL_EPSILON)),
            0.0);

        //An overeager floating point optimiser may fail the following test
        //by turning the stimulus off too early
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(5100.5),
            0.0);

        //Made more sloppy
        TS_ASSERT_EQUALS(regular_stimulus.GetStimulus(5100.5*(1+2*DBL_EPSILON)),
            0.0);
    }

    //void TestBasicFmod() removed since the exact floating point behaviour
    //is too difficult to reproduce

    // SumStimulus is redundant as it is a special case of MultiStimulus, so we have refactored
    // the code by deleting SumStimulus. This test shows how to use MultiStimulus instead
    void TestSumStimulus()
    {
        MultiStimulus multi_stim;
        boost::shared_ptr<SimpleStimulus> p1(new SimpleStimulus(2,1,0));
        boost::shared_ptr<SimpleStimulus> p2(new SimpleStimulus(3,1,3));

        multi_stim.AddStimulus(p1);
        multi_stim.AddStimulus(p2);

        TimeStepper t(0,10,1);
        while (!t.IsTimeAtEnd())
        {
            TS_ASSERT_EQUALS( multi_stim.GetStimulus(t.GetTime()),
                              p1->GetStimulus(t.GetTime()) + p2->GetStimulus(t.GetTime()) );
            t.AdvanceOneTimeStep();
        }
    }

    void TestZeroStimulus()
    {
        ZeroStimulus zero_stim;
        TS_ASSERT_EQUALS( zero_stim.GetStimulus(1), 0);
    }

    void TestMultiStimulus()
    {
        MultiStimulus multi_stim;

        // No stimulus after creation.
        TS_ASSERT_EQUALS( multi_stim.GetStimulus(1.0), 0.0);

        boost::shared_ptr<SimpleStimulus> p_init_stim_a(new SimpleStimulus(2,1,0));
        boost::shared_ptr<SimpleStimulus> p_init_stim_b(new SimpleStimulus(3,1,30));
        // RegularStimulus result at a boundary point isn't the same for Default and IntelProduction build
        // (in fact it gives different answers to "fmod" within the IntelProduction test)
        boost::shared_ptr<RegularStimulus> p_regular_stim(new RegularStimulus(2.0, 1.0, 1.0/0.15, 1));

        multi_stim.AddStimulus(p_init_stim_a);
        multi_stim.AddStimulus(p_init_stim_b);
        multi_stim.AddStimulus(p_regular_stim);

        TimeStepper t(0,100,1);
        while (!t.IsTimeAtEnd())
        {
            double time=t.GetTime();
            // Stimulus equals to the sum of the individual stimuli
            TS_ASSERT_EQUALS( multi_stim.GetStimulus(time),
                              p_init_stim_a->GetStimulus(time)+
                              p_init_stim_b->GetStimulus(time)+
                              p_regular_stim->GetStimulus(time)
                            );
            t.AdvanceOneTimeStep();
        }

    }

    void TestComplicatedStimulus()
    {
        MultiStimulus multi_stim;

        boost::shared_ptr<RegularStimulus> pr1(new RegularStimulus(1,10,20,100,200));
        // First stimulus applies pulses for 10ms every 20ms between t=100 and t=200.
        boost::shared_ptr<RegularStimulus> pr2(new RegularStimulus(2,20,40,300,400));
        // Second stimulus applies pulses for 20ms every 40ms between t=300 and t=400.

        multi_stim.AddStimulus(pr1);
        multi_stim.AddStimulus(pr2);

        // Test on 0.5s so we avoid worrying about <= or < and stuff!
        for (double time=0.5; time<500.0; time=time+1)
        {
            if (time>=100 && time <=200)
            {   // first stimulus active when mod(t-100,20) < 10
                if (fmod(time-100,20) < 10)
                {
                    TS_ASSERT_DELTA(multi_stim.GetStimulus(time), 1.0, 1e-9);
                }
                else
                {
                    TS_ASSERT_DELTA(multi_stim.GetStimulus(time), 0.0, 1e-9);
                }
            }
            else if (time>=300 && time <=400)
            {   // second stimulus active when mod(t-300,40) < 20
                if (fmod(time-300,40) < 20)
                {
                    TS_ASSERT_DELTA(multi_stim.GetStimulus(time), 2.0, 1e-9);
                }
                else
                {
                    TS_ASSERT_DELTA(multi_stim.GetStimulus(time), 0.0, 1e-9);
                }
            }
            else
            {   //zero
                TS_ASSERT_DELTA(multi_stim.GetStimulus(time), 0.0, 1e-9);
            }
        }
    }

    void TestArchivingStimuli() throw(Exception)
    {
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "stimuli.arch";

        // Create and archive simulation time
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set up a zero stimulus
            AbstractStimulusFunction* const p_zero_stimulus = new ZeroStimulus;

            double magnitude = 1.0;
            double duration  = 0.5;  // ms
            double applied_time = 100.0;
            AbstractStimulusFunction* const p_simple_stimulus = new SimpleStimulus(magnitude, duration, applied_time);

            double magnitude_of_stimulus = 1.0;
            double duration_of_stimulus  = 0.5;  // ms
            double period = 1000.0; // 1s
            double when = 100.0;
            AbstractStimulusFunction* const p_regular_stimulus = new RegularStimulus(magnitude_of_stimulus,
                                         duration_of_stimulus,
                                         period,
                                         when);

            MultiStimulus* p_multiple_stimulus = new MultiStimulus;
            boost::shared_ptr<SimpleStimulus> p1(new SimpleStimulus(2,1,0));
            boost::shared_ptr<SimpleStimulus> p2(new SimpleStimulus(3,1,3));
            p_multiple_stimulus->AddStimulus(p1);
            p_multiple_stimulus->AddStimulus(p2);

            AbstractStimulusFunction* const p_multiple_stimulus2 = p_multiple_stimulus;   // make const now we've added stimuli

            // Should always archive a pointer
            output_arch << p_zero_stimulus;
            output_arch << p_simple_stimulus;
            output_arch << p_regular_stimulus;
            output_arch << p_multiple_stimulus2;

            // Change stimulus a bit
            delete p_zero_stimulus;
            delete p_simple_stimulus;
            delete p_regular_stimulus;
            delete p_multiple_stimulus;
        }

        // Restore
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Create a pointer
            AbstractStimulusFunction* p_zero;
            AbstractStimulusFunction* p_simple;
            AbstractStimulusFunction* p_regular;
            AbstractStimulusFunction* p_multiple;
            input_arch >> p_zero;
            input_arch >> p_simple;
            input_arch >> p_regular;
            input_arch >> p_multiple;

            // Check the stimulus now has the properties of the one archived above.
            TS_ASSERT_DELTA(p_zero->GetStimulus(0.0), 0.0, 1e-9);

            TS_ASSERT_DELTA(p_simple->GetStimulus(0.0), 0.0, 1e-9);
            TS_ASSERT_DELTA(p_simple->GetStimulus(100.1), 1.0, 1e-9);
            TS_ASSERT_DELTA(p_simple->GetStimulus(100.6), 0.0, 1e-9);

            TS_ASSERT_DELTA(p_regular->GetStimulus(0.0), 0.0, 1e-9);
            TS_ASSERT_DELTA(p_regular->GetStimulus(100.1), 1.0, 1e-9);
            TS_ASSERT_DELTA(p_regular->GetStimulus(100.6), 0.0, 1e-9);
            TS_ASSERT_DELTA(p_regular->GetStimulus(1100.1), 1.0, 1e-9);
            TS_ASSERT_DELTA(p_regular->GetStimulus(1100.6), 0.0, 1e-9);

            TimeStepper t(0,10,1);
            SimpleStimulus r1(2,1,0);
            SimpleStimulus r2(3,1,3);
            while (!t.IsTimeAtEnd())
            {
                TS_ASSERT_EQUALS( p_multiple->GetStimulus(t.GetTime()),
                                  r1.GetStimulus(t.GetTime())+
                                  r2.GetStimulus(t.GetTime()) );
                t.AdvanceOneTimeStep();
            }

            delete p_zero;
            delete p_simple;
            delete p_regular;
            delete p_multiple;
        }
    }

};

#endif /*TESTSTIMULUS_HPP_*/
