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
