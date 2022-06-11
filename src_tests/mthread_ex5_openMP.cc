// from https://cplusplus.com/forum/general/253103/

// /usr/local/Cellar/gcc/11.3.0/bin/g++-11 -I/usr/local/include -fopenmp mthread_ex5_openMP.cc -o mthread_ex5_openMP
// /usr/local/Cellar/gcc@12/12.1.0/bin/g++-12 -I/usr/local/include -fopenmp mthread_ex5_openMP.cc -o mthread_ex5_openMP
// clang -Xclang -fopenmp -lomp -std=c++20 -lc++abi -lstdc++ -pthread -O3 mthread_ex5_openMP.cc -o mthread_ex5_openMP


#include <iostream>
#include <ctime>
#include "omp.h"
using namespace std;

int main()
{
   const int N = 100000000;
   static double A[N];
   for ( int i = 0, sign = 1; i < N; i++, sign = -sign )
   {
      A[i] = sign * 4.0 / ( 1.0 + 2.0 * i );
   }

   for ( int nThreads = 1; nThreads <= 8; nThreads++ )
   {
      clock_t start = clock();
      double sum = 0.0;

      omp_set_num_threads( nThreads );                     // OpenMP
      #pragma omp parallel for reduction( + : sum )        // OpenMP
      for ( int i = 0; i < N; i++ )
      {
         sum += A[i];
      }

      clock_t end = clock();

      cout << "Threads: " << nThreads << "     Sum: " << sum << "     Time: " << (double)( end - start ) / CLOCKS_PER_SEC << " s\n";
   }
}
