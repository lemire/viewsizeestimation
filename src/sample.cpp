/**
 * This code will sample some rows out of a large text file.
 */
#include <algorithm>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <time.h>

using namespace std;
typedef unsigned int uint;

int main(int params, char **args) {
  float p(1.0); // Sampling a percentage of data
  uint n(0);    // Sampling a number of data
  map<int, bool> randomnumber;
  string samplefile("sample.out.txt"); // default output
  int r = 4567;                        // rand seed

  cout << "This program generates samples from an input data file" << endl;

  if (params <= 3) {
    cerr << "You must provide some parameters including  a sample size and a "
            "filename "
         << endl;
    cerr << "Usage: ./sample -n 1000 myfilen.txt" << endl;
    return -1;
  }

  string datafilename(args[params - 1]); // Data file name
  cout << "Data file name:  " << datafilename.c_str() << endl;

  ifstream filein(datafilename.c_str()); // Opening

  if (!filein) {
    cout << "Can't open the damned file. " << endl;
    return -10;
  }
  clock_t start, finish;
  start = clock();
  // Reading and checking the program parareters
  for (int i = 1; i < params - 1; ++i) {
    if (strcmp(args[i], "-n") == 0) {
      n = atoi(args[++i]);
      cout << "Number of facts in the sample you are generating: " << n << endl;
    }
    if (strcmp(args[i], "-p") == 0) {
      p = (float)atof(args[++i]);
      cout << "Percentage of sampling: " << p << "%" << endl;
      n = 0;
    }
    if (strcmp(args[i], "-o") == 0) {
      samplefile = args[++i];
      cout << "Name of the sample file: " << samplefile.c_str() << endl;
    }
    if (strcmp(args[i], "-r") == 0) {
      r = atoi(args[++i]);
      cout << "Random seed: " << r << endl;
    }
  }

  // srand(time(NULL));
  srand(r);

  uint N = count(istreambuf_iterator<char>(filein), istreambuf_iterator<char>(),
                 '\n'); // counting the number of lines
  filein.close();

  cout << "Number of facts: " << N << endl;
  uint ec = n;
  if (n == 0)
    ec = (uint)(p * N / 100); // if -n is not set

  while (randomnumber.size() < ec) {
    int randomvalue = 1 + (int)((double)rand() / ((double)RAND_MAX + 1) * N);
    randomnumber.insert(make_pair(randomvalue, true));
  }
  ofstream fileout(samplefile.c_str());

  if (!fileout) {
    cerr << "Cannot create the output file" << endl;
    return 11;
  }
  filein.open(datafilename.c_str(), ifstream::in); // open again fileib
  string line;
  uint linenumber(0);
  getline(filein, line); // attribute names
  fileout << line << endl;
  while (filein) {
    getline(filein, line);
    ++linenumber;
    map<int, bool>::iterator it = randomnumber.find(linenumber);
    if (it != randomnumber.end()) {
      fileout << line << endl;
      randomnumber.erase(it);
    }
    if (randomnumber.size() == 0)
      break;
  }
  filein.close();
  fileout.close();
  finish = clock();
  cout << "Size of the generated sample is: " << ec << endl;
  cout << "Time used to read data: "
       << (double)(finish - start) / CLOCKS_PER_SEC << " seconds or "
       << (finish - start) << " ticks " << endl;
  return 0;
}
