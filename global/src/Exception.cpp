/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Exception.hpp"

Exception::Exception(std::string message,
                     std::string filename, const unsigned rLineNumber)
{
    std::stringstream line_number;
    line_number << rLineNumber;
    
    mMessage = std::string("\nChaste error: ") + filename + ":"  + line_number.str()  + ": " + message;
}


std::string Exception::GetMessage() const
{
    return mMessage;
}
