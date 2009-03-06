#ifndef DEBUG_HPP_
#define DEBUG_HPP_

#include <iostream>
#include <cassert>

// A bunch of useful macros for debugging. Note, these should be removed from source
// code when committing.

/* Print the given message */
#define TRACE(stuff) std::cout << "DEBUG: " << stuff << std::endl << std::flush;

/* Print the name and value of the given variables */
#define PRINT_VARIABLE(var) std::cout << "DEBUG: " #var " = " << var << std::endl << std::flush;
#define PRINT_VARIABLES(var1,var2) std::cout << "DEBUG: " #var1 " = " << var1 << ", " #var2 " = " << var2 << std::endl << std::flush;
#define PRINT_3_VARIABLES(var1,var2,var3) std::cout << "DEBUG: " #var1 " = " << var1 << ", " #var2 " = " << var2 << ", " #var3 " = " << var3 << std::endl << std::flush;
#define PRINT_4_VARIABLES(var1,var2,var3,var4) std::cout << "DEBUG: " #var1 " = " << var1 << ", " #var2 " = " << var2 << ", " #var3 " = " << var3 << ", " #var4 " = " << var4 << std::endl << std::flush;

/* Quit (assert(0)) on the n-th time this line is reached, for the given n */ 
#define QUIT_AFTER_N_VISITS(n) { static unsigned counter=1; if(counter++==(n)) {TRACE("User-forced quit."); assert(0);} }

/* Print how many times this line has been reached, everytime it is reached */
#define HOW_MANY_TIMES_HERE(message) { static unsigned counter=1; std::cout << "DEBUG: Num times here ("<< message << "): " << counter++ << std::endl << std::flush; }

#endif /*DEBUG_HPP_*/
