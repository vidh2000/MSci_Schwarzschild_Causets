#include <iostream>
#include <cstdlib>
using namespace std;

// Compile with g++ for it to work, then:
// ./a.out 2 4 --> will give Sum = 6 etc.

int main(int argc, char* argv[]) {
   for(int i = 1; i < argc; i++)
      cout << atoi(argv[i]) << endl;

    int a = atoi(argv[1]);
    int b = atoi(argv[2]);
    int sum = a+b;
    cout << "Sum = "<< sum << endl;
    return 0;
}