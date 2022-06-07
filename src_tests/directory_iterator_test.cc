//
// from https://stackoverflow.com/questions/612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c/37494654#37494654
//

// compile with
// clang -std=c++17 -lc++abi -lstdc++ directory_iterator_test.cc -o directory_iterator_test

#include <string>
#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 2) {
      fprintf(stderr, "Usage: %s <pathname>\n", argv[0]);
      exit(EXIT_FAILURE);
    }

    string dir = argv[1];
    cout << "\n\n========================================\n" <<
      "Using directory_iterator() to list files in the directory " << dir << endl;
    for (const auto & entry : fs::directory_iterator(dir))
      cout << entry.path() << endl;

    cout << endl;

    return EXIT_SUCCESS;
}
