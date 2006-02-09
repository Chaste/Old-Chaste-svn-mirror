#include <iostream>

int main()
{
   int n = 0, t = 0, i = 0;

   std::cout << "200\t3\t0" << std::endl;

   while (n < 100)
   {
      if (t != 10)
      {
	// Top-left - bottom-right orientation

         std::cout << 2*n << "\t" << i << "\t" << i+1 << "\t" << i+11 << std::endl;
         std::cout << 2*n+1 << "\t" << i+1 << "\t" << i+12 << "\t" << i+11 << std::endl;

         // Bottom-left - top-right orientation

//         std::cout << 2*n << "\t" << i << "\t" << i+1 << "\t" << i+12 << std::endl;
//         std::cout << 2*n+1 << "\t" << i << "\t" << i+12 << "\t" << i+11 << std::endl;

         // Criss-cross shape - !!remember to change 200 to 400 in line above outside loop

	//         std::cout << 4*n << "\t" << i << "\t" << i+1 << "\t" << n+121 << std::endl;
	//         std::cout << 4*n+1 << "\t" << i+1 << "\t" << i+12 << "\t" << n+121 << std::endl;
	//std::cout << 4*n+2 << "\t" << i+12 << "\t" << i+11 << "\t" << n+121 << std::endl;
	//std::cout << 4*n+3 << "\t" << i+11 << "\t" << i << "\t" << n+121 << std::endl;

         t++;
         n++;
      }
      else
      {
         t = 0;
      }

      i++;
   }

   return 0;
}
