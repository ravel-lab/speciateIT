// from https://cplusplus.com/forum/general/253103/
//
// clang -std=c++20 -lc++abi -lstdc++ -pthread -O3 -DNDEBUG mthread_ex4.cc -o mthread_ex4

#include <iostream>
#include <numeric>
#include <future>
#include <ctime>
#include <chrono>
#include <string>

auto psum( const double* begin, const double* end, std::size_t seg_size )
{
    const std::size_t n =  end - begin ;
    if( n <= seg_size ) return std::accumulate( begin, end, 0.0 ) ;

    const double* next_seg_begin = begin + seg_size ;

    // compute the sum of the tail asynchronously (potentially in another thread)
    auto future_sum_tail = std::async( std::launch::async, psum, next_seg_begin, end, seg_size ) ;

    // compute the sum of the first segment in this thread
    const auto sum_head = std::accumulate( begin, begin+seg_size, 0.0 ) ;

    return sum_head + future_sum_tail.get() ;
}

// high_resolution_clock::time_point t1 = high_resolution_clock::now();
//  for(int n=0; n<100000000; ++n) sum += A[n];
//  high_resolution_clock::time_point t2 = high_resolution_clock::now();
//   duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

struct timer
{
    timer( std::string str ) : text(str) {}

    ~timer()
    {
        using namespace std::chrono ;
        const auto finish_p = std::clock() ;
        const auto finish_w = steady_clock::now() ;
        std::cout << '\n' << text << "\n---------\n"
                  << " processor: " << (finish_p-start_p) * 1000.0 / CLOCKS_PER_SEC << " milliseconds.\n"
                  << "wall clock: " << duration_cast<milliseconds>(finish_w-start_w).count() << " milliseconds.\n" ;
    }

    const std::string text ;
    const std::clock_t start_p = std::clock() ;
    const std::chrono::steady_clock::time_point start_w =  std::chrono::steady_clock::now() ;
};

int main()
{
    const std::size_t N = 8'000'000'000 ;
    static double seq[N];
    std::iota( seq, seq+N, 0 ) ;

    double sum ;
    {
        timer t { "linear" };
        sum = std::accumulate( seq, seq+N, 0.0 ) ;
    }
    std::cout << "sum: " << sum << '\n' ;

    {
        timer t { "async (12 segments)" };
        sum = psum( seq, seq+N, N/12 ) ;
    }
    std::cout << "sum: " << sum << '\n' ;
}

// linear
// ---------
//  processor: 95792.7 milliseconds.
// wall clock: 190165 milliseconds.
// sum: 1.74011e+17

// async (12 segments)
// ---------
//  processor: 217786 milliseconds.
// wall clock: 56725 milliseconds.
// sum: 1.74011e+17

// 190165 / 56725 = 3.352402
