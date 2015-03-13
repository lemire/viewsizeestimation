/**
* the goal of this file is to test out Flajolet's theoretical result
*/

#include <fstream>
#include <iostream>
#include <cstdlib>
using namespace std;

int main(int params, char ** args) {
  if (params <= 2) {
    cout << "You must provide two parameters, number of distinct values and length of stream " << endl;
    return -1;
  } 
  //cout << "parsing using " << sizeof(long long) << " bytes" << endl;
  long long n = atoll(args[2]);
  long long k = atoll(args[1]);
  if(k > n) {
    cout << "The first parameter must, logically, be smaller than the second one " << endl;
    return -1;
  }
  cout << "unique" << endl; // name of dim is always "unique"
  for(long long i = 0; i < n ; ++i) {
   cout << (i % k) << ',' << endl;
  }
  return 0;
}
