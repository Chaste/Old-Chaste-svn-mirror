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

#ifndef SENSITIVITYMODIFIERS_HPP_
#define SENSITIVITYMODIFIERS_HPP_

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

/**
 * This class specifies the calc() method which all modifiers must call to intervene on the
 * right hand side of a cell model.
 */
class AbstractSensitivityModifier
{
  public:
    /**
     * Default constructor.
     */
    AbstractSensitivityModifier(void){}

    /**
     * Default destructor.
     */
    virtual ~AbstractSensitivityModifier(){};

    /**
     * Pure virtual function which must be overriden in subclasses.
     */
    virtual double calc(double param, double time) = 0;
};

/**
 * This class allows modification of parameters by a scale factor.
 */
class FactorModifier : public AbstractSensitivityModifier
{
private:
    /** Factor to multiply parameter of interest by. */
    double mFactor;

public:
    /**
     * Constructor
     * @param factor  Scale factor to use, defaults to 1 (no effect)
     */
    FactorModifier(double factor=1)
        :mFactor(factor)
    {}

    virtual double calc(double param, double time)
    {
        return (param * mFactor);
    }
};

/**
 * This is just an example class to show how you might specify a custom
 * modifier to change a parameter through time. In this case it implements
 * a sin(time)*default_parameter factor modifier.
 */
class TimeModifier : public AbstractSensitivityModifier
{
public:
    /**
     * Constructor
     */
    TimeModifier(){}

    /**
     * @return The value that the parameter should take at this time.
     */
    virtual double calc(double param, double time)
    {
        return param * sin(time);
    }
};

/**
 * This class just returns a fixed value, regardless of the parameter's default or the time.
 */
class FixedModifier : public AbstractSensitivityModifier
{
private:
    /** Fixed value to clamp parameter at */
    double mValue;

public:
    /**
     * Constructor
     * @param value  The fixed value to use.
     */
    FixedModifier(double value)
      :mValue(value)
    {
    }

    /**
     * @return The fixed value (ignores inputs).
     */
    virtual double calc(double param, double time)
    {
        return mValue;
    }
};

/**
 * This class returns the parameter's default value and does not modify it.
 */
class DummyModifier : public AbstractSensitivityModifier
{
private:

public:
    /**
     * Default Constructor
     */
    DummyModifier(){}

    /**
     * This calculate does nothing and returns the parameter 'unharmed'.
     */
    virtual double calc(double param, double time)
    {
        return param;
    }
};


#endif  //SENSITIVITYMODIFIERS_HPP_
