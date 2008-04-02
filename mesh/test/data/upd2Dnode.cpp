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

#include <iostream>
#include <math.h>

int main()
{
    int n = 0;
    
    std::cout << std::endl;
    
    while (n < 100)
    {
        std::cout << " " << 121+n << "\t" << 0.005+(n%10)*0.01 << "\t" << 0.005+floor(n/10)*0.01 << "\t0" << std::endl;
        
        n++;
    }
    
    return 0;
}
