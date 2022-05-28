//
// clang -std=c++11 -lc++abi -lstdc++ -pthread -O3 -DNDEBUG mthread_ex.cpp -o mthread_ex
//

#include <iostream>
#include <thread>

void f()
{
  std::cout << "Hello World\n";
}

int main()
{
  std::thread t(f);
  t.join();
}
