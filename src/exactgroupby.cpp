#include <algorithm>
#include <cassert> //hopefully cassert and cmath exists under vc++?
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <time.h>
#include <vector>

using namespace std;

enum { MAXLENGTHOFALINE = 256 };

typedef unsigned int uint;
typedef unsigned long int uint64;

inline void split(const string &str, vector<string> &tokens,
                  const string &delimiters = " ") {
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

template <class T> class StringConverter {
public:
  StringConverter() {}
  T toType(const char *input) const { return (T)input; }
};

template <> class StringConverter<int> {
public:
  StringConverter() {}
  int toType(const char *input) const { return atoi(input); }
};

template <> class StringConverter<char> {
public:
  StringConverter() {}
  char toType(const char *input) const { return (char)atoi(input); }
};

template <class T>
uint64 ExactGroupBy(istream &in, vector<string> &groupby, int &N,
                    const char delimiter = ',') {
  const uint n = (uint)groupby.size(); // number of dimensions
  StringConverter<T> sc;
  map<vector<T>, bool> histogram;

  string line; // Read a line of data
  vector<T> attribute_value(n);
  vector<string> attribute_list; // Pour stocker les attributs de dimension et
                                 // leur valeur de chaque faits
  vector<unsigned int> dimindice; // Pour stocker les indices des colonnes
                                  // correspondant aux attributs du group by
  getline(in, line);
  string delimitertostring(1, delimiter);
  split(line, attribute_list, delimitertostring);

  for (unsigned int j = 0; j < groupby.size(); ++j)
    for (unsigned int i = 0; i < attribute_list.size(); ++i)
      if (groupby[j] == attribute_list[i]) {
        dimindice.push_back(i);
        break;
      }

  sort(dimindice.begin(),
       dimindice.end()); // Ascending sort of the dimension indices
  assert(dimindice.size() == n);

  char chline[256];
  while (in.getline(chline, MAXLENGTHOFALINE, delimiter)) {
    ++N;
    uint col(0);
    bool endline(false);
    for (uint k = 0; k < dimindice.size(); ++k) {
      while (col < dimindice[k]) {
        ++col;
        if (col == attribute_list.size() - 1) {
          in.getline(chline, MAXLENGTHOFALINE);
          endline = true;
        } else
          in.getline(chline, MAXLENGTHOFALINE, delimiter);
      }
      attribute_value[k] = sc.toType(chline);
    }
    if (!endline)
      in.getline(chline, MAXLENGTHOFALINE);

    typename map<vector<T>, bool>::iterator it =
        histogram.find(attribute_value);
    if (it == histogram.end())
      histogram[attribute_value] = 1;
  }
  return (uint64)histogram.size();
}

int main(int params, char **args) {
  char delimiter(',');
  vector<string> groupby; // Group by attributes
  bool eightbits = false;
  bool thirtytwobits = false;
  bool verbose = false;
  int N(0); // Total number of facts

  if (params <= 3) {
    cout << " You must provide some parameters including a groupby and a "
            "filename "
         << endl;
    return -1;
  }

  string datafilename(args[params - 1]); // Data file name

  ifstream filein(datafilename.c_str()); // Opening

  // Reading and checking the program parareters
  for (int i = 1; i < params - 1; ++i) {
    if (strcmp(args[i], "-v") == 0) {
      verbose = true;
    }
    if (strcmp(args[i], "--groupby") == 0) {
      string sgroupby; // Group by attributes as a string

      if (params - i > 1) {
        sgroupby = args[++i];
        split(sgroupby, groupby, ",");
        if (verbose == true) {
          cout << "Group by attributes: ";
          for (uint l = 0; l < groupby.size(); ++l)
            cout << groupby[l] << " ";
          cout << endl;
        }
      } else {
        cerr << "--groupby expects a string" << endl;
        return -2;
      }
    }
    if (strcmp(args[i], "--8bits") == 0) {
      eightbits = true;
      if (verbose == true)
        cout << "Assuming attribute values fit in 8 bits" << endl;
    }
    if (strcmp(args[i], "--32bits") == 0) {
      thirtytwobits = true;
      if (verbose == true)
        cout << "Assuming attribute values fit in 32 bits" << endl;
    }
    if (strcmp(args[i], "-d") == 0) {
      string stringdelim = args[++i];
      delimiter = stringdelim.c_str()[0];
      if (verbose == true)
        cout << "Delimiter of group by values: " << delimiter << endl;
    }
  }
  if (verbose == true)
    cout << "Data file name:  " << datafilename.c_str() << endl;
  if (!filein) {
    cerr << "Can't open the damned file. " << endl;
    return -10;
  }

  int answer = 0; // the count-distinct estimate

  clock_t start, finish;
  start = clock();

  if (eightbits)
    answer = ExactGroupBy<char>(filein, groupby, N, delimiter);
  else if (thirtytwobits)
    answer = ExactGroupBy<int>(filein, groupby, N, delimiter);
  else {
    if (verbose == true)
      cout << "Assuming attribute values are strings???" << endl;
    if (verbose == true)
      cout << "Please try --8bits or --32bits instead, this is going to be slow"
           << endl;
    answer = ExactGroupBy<string>(filein, groupby, N, delimiter);
  }

  finish = clock();
  if (verbose == true) {
    cout << "Exact group by size: " << answer << endl;
    cout << "Time used to read data: "
         << (double)(finish - start) / CLOCKS_PER_SEC << " seconds or "
         << (finish - start) << " ticks " << endl;
    cout << "Nb of tuples: " << N << endl;
  } else {
    cout << answer << " " << (double)(finish - start) / CLOCKS_PER_SEC << endl;
  }
  return 0;
}
