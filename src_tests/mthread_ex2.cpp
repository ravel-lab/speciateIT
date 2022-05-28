//
// clang -std=c++11 -lc++abi -lstdc++ -pthread -O3 -DNDEBUG mthread_ex2.cpp -o mthread_ex2
//

#include <iostream>
#include <thread>
#include <string>
#include <mutex>

static const int num_threads = 10;

std::mutex mtx;

//This function will be called from a thread
void call_from_thread(int tid) {
  mtx.lock();
  std::cout << "Launched by thread " << tid << std::endl;
  mtx.unlock();
}

int main() {
    std::thread t[num_threads];

    //Launch a group of threads
    for (int i = 0; i < num_threads; ++i) {
        t[i] = std::thread(call_from_thread, i);
    }

    mtx.lock();
    std::cout << "Launched from the main\n";
    mtx.unlock();

    //Join the threads with the main thread
    for (int i = 0; i < num_threads; ++i) {
        t[i].join();
    }

    return 0;
}
