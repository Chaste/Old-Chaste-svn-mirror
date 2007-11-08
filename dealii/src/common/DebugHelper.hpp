#ifndef DEBUGHELPER_HPP_
#define DEBUGHELPER_HPP_

/** This is a noddy file just to have some code (ie a source file that
 *  isn't templated in the Dealii folder, else the dealii library will 
 *  empty
 */

#include <iostream>
#include <string>

// eventually could add print methods for dealii classes, eg Print(Vector<double> vec) 
 
class DebugHelper
{
    static void Print(std::string message);
};


#endif /*DEBUGHELPER_HPP_*/
