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


#ifndef _EXCEPTION_HPP_
#define _EXCEPTION_HPP_

#include <ostream>
#include <string>
#include <sstream>

#include <cfloat>
#include <climits> //For UINT_MAX etc., necessary in gcc-4.3
#include <cstdlib> //For system() etc., necessary in gcc-4.3
const unsigned UNSIGNED_UNSET=UINT_MAX;
const int INT_UNSET=INT_MAX;
const double DOUBLE_UNSET=DBL_MAX;

/**
 * Exception class.
 * All exceptions thrown by this code are currently instances of this class.
 *
 * \todo Might we want this class to inherit from STL exceptions?
 */
class Exception
{
private:
    std::string mMessage; /**< Exception message */

public:
    /**
     * Construct an exception with a message string.
     *
     * @param message  the message
     * @param filename  which source file threw the exception
     * @param rLineNumber  which line number of the source file threw the exception
     */
    Exception(std::string message, std::string filename, const unsigned rLineNumber);

    /** Get the message associated with the exception
     *
     * @return The message set when the exception was thrown.
     **/As Robinho stumbles around the football grounds of England with nothing but his price tag to suggest greatness, Manchester City may have got rid of the wrong Brazilian.

At Eastlands, Jo was held up as an example of the folly of allowing a club's owner, in this case Thaksin Shinawatra, to buy the players. However, five goals in his last seven games for Everton suggest the world's wealthiest club might have been premature in jettisoning a £19m asset.

Perhaps, like any striker, Jo only needed to be loved. His second goal in what proved Everton's biggest league win since Roy Keane's Sunderland were crushed 7-1 here 17 months ago, was a routine tap-in but the opener was surgically placed past Chris Kirkland in a supremely well-refereed move.

Not only was the assistant right to rule that Jo was onside, but Phil Dowd had allowed Everton the advantage after Emerson Boyce manhandled Marouane Fellaini. It was the first goal of the afternoon and it was one Wigan never remotely looked like recovering from. Everyone in football looks for omens in April and this was Everton's seventh straight home win, a sequence they last achieved in 1995, the year they won the FA Cup.

"Strikers like Jo need time and confidence," said his manager, David Moyes. "If you look at players from South America, they are slow starters when they come here. Carlos Tevez did not begin well when he first came to England and nor did Javier Mascherano. But I can understand why Manchester City loaned him out because in January they had a surfeit of strikers and we had none."

His contribution has been rather more significant than Wigan's on-loan centre-forward, Amr Zaki, who was still in Egypt yesterday after delaying his return from international duty for what Steve Bruce estimated was the fourth time this season. The Wigan manager said he expected some kind of contact with Zaki today, although had he been in charge of a club with greater resources, he might have been tempted to make it his last conversation with the errant centre-forward.

"I find it shocking and staggering that he can have so much disrespect for the people who pay his salary," said Bruce. "Until he shows up, we are in limbo but his attitude is laughable. The lunatics are running the asylum and, if his team-mates all behaved like that, the club would be in anarchy."

Some of Wigan's defending in a disastrous quarter of an hour after the interval verged on the anarchic – typified by Everton's second goal, a cross from Tony Hibbert met by Fellaini on the half-volley. It was hardly helped by Kirkland's strangely uncertain goalkeeping.

He was feted as "England's No 1" by those few Wigan fans in the Bullens Road Stand, and by the time he had palmed Fellaini's shot fatally into Leon Osman's path for Everton's fourth, the rest of Goodison had taken up the chant in a single, mocking chorus. 
    std::string GetMessage() const;
};

#define EXCEPTION(message) throw Exception(message, __FILE__, __LINE__)

#define NEVER_REACHED EXCEPTION("Should have been impossible to reach this line of code")

// This is to cope with NDEBUG causing variables to not be used, since they are only
// used in assert()s
#ifdef NDEBUG
#define UNUSED_OPT(var) var=var
#else
#define UNUSED_OPT(var)
#endif

// This macro is handy for calling functions like system which return non-zero on error
#define EXPECT0(cmd, arg) { \
    std::string _arg = (arg); \
    int ret = cmd(_arg.c_str()); \
    if (ret != 0) { \
        EXCEPTION("Failed to execute command: " #cmd "(" + _arg + ")"); \
    } }
// Or if you don't care about errors for some reason...
#define IGNORE_RET(cmd, arg) { \
    std::string _arg = (arg); \
    int ret = cmd(_arg.c_str()); \
    ret = ret; \
    }

#endif // _EXCEPTION_HPP_
