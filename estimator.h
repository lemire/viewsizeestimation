#ifndef ESTIMATOR_H
#define ESTIMATOR_H

#include "randomhasher.h"
#include "trailingzeros.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctype.h>
#include <assert.h>
#include <sstream>
#include <algorithm>

using namespace std;


enum{MAXLENGTHOFALINE=256,stoverbose = false, verbose = false, readonly=true};

typedef unsigned int uint;



/* The upper 32 bits will be count, the lower are the traditional "number of trailing zeros" */

inline void split(const string& str, vector<string>& tokens, const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
  }
}

template <class T>
class StringConverter {
        public:
                StringConverter(){}
                T toType(const char * input) const {
                        return (T) input;
                }
};
template <>
class StringConverter<int> {
        public:
                StringConverter(){}
                int toType(const char * input) const {
                        return  atoi(input);
                }
};
template <>
class StringConverter<char> {
        public:
                StringConverter(){}
                char toType(const char * input) const {
                        return (char) atoi(input);
                }
};



#include "stogibbons.h"

template <class T>
uint64 GibbonsTirthapura(istream & in, vector<string> & groupby, int& N, int L=19, int storagelimit=512, char& delimiter = ',', int rng = rand_rng) 
{
  map< vector<T>, uint64> buffer;
  const uint n = groupby.size();//number of dimensions
  StringConverter<T> sc;

  // Estimating simpOuie functions on the union of data streams by Gibbons and Tirthapura
  TupleFactHash<T> cph(n,L,rng);
  uint t = 0;
  uint64 mask = (1 << (t)) - 1;   // which is initially 0
  string line; //Read a line of data
  vector<T> attribute_value(n); 
  vector<string> attribute_list;//Pour stocker les attributs de dimension et leur valeur de chaque faits
  vector<unsigned int> dimindice;//Pour stocker les indices des colonnes correspondant aux attributs du group by

  N=0;
  assert(buffer.size() == 0);  // if not, we needed to clear...

  //Récupérer les indices des colones des dimensions
  getline (in, line);//La première ligne du fichier de données contient les noms des dimensions

  string delimitertostring (1,delimiter);
  split (line,attribute_list,delimitertostring);
  
  for (unsigned int j=0; j < groupby.size(); ++j)
    for (unsigned int i=0; i < attribute_list.size (); ++i)
      if ( groupby[j] == attribute_list[i]) {
        dimindice.push_back (i);
        break;
      }
  sort(dimindice.begin(), dimindice.end());//Ascending sort of the dimension indices 

  assert(dimindice.size() == n);
  char chline[256];
  while (in.getline(chline,MAXLENGTHOFALINE,delimiter))
  {
    ++N;

    //cout << endl << "premier " << " " << chline << endl;
    uint col(0);
    bool endline(false);
    for(uint k = 0; k < dimindice.size(); ++k) 
    {
      while( col < dimindice[k])  
      {
        ++col;
        if (col == attribute_list.size()-1)
        {
          in.getline(chline, MAXLENGTHOFALINE);
          endline = true;
        }
        else 
          in.getline(chline, MAXLENGTHOFALINE,delimiter);
      }
      //cout << endl << "trouve " << k << " " << chline << endl;
      attribute_value[k] = sc.toType(chline);
    }
    if (!endline) in.getline(chline, MAXLENGTHOFALINE);

    uint64 oldhash = cph.slow(attribute_value); // hashing
    //buffer[attribute_value] = (1LL << 32) | trailingzeros(oldhash,L);  // count= 1 
    uint projection = oldhash & mask;
    if (projection == 0) 
    {
      typename map< vector<T>,uint64>::iterator it = buffer.find(attribute_value); 
      if (it == buffer.end()) 
      {
        //buffer[attribute_value] = (1LL << 32) | trailingzeros(oldhash,L);
        buffer[attribute_value] = trailingzeros(oldhash,L);//(1LL << 32) |
        while( (buffer.size() >(uint) storagelimit) ) 
        {
          t+= 1;
          mask = (1 << (t)) - 1;
          assert(t < (uint)L);
          for(typename map<vector<T>, uint64>::iterator j = buffer.begin(); j != buffer.end();) 
          //if( (j->second & 0xffffffff)  < t) buffer.erase(j++); 
          if( (j->second )  < t) buffer.erase(j++); // & 0xffffffff
          else ++j;
        }
      } //else // previously seen, so bump count by 1
      //buffer[attribute_value] += (1LL << 32);
    }
  }
  return buffer.size() * (1 << t);
}



//Probabilistic counting: several hash functions
template <class T>
uint64 ProbabilisticCounting (istream & in, vector<string> & groupby, int& N, int L=19, const char& delimiter = ',', int rng = rand_rng)
{
  const float phi = 0.77351;
  vector<bool> bitmap(L);
  for (uint i=0; i < L; ++i)
    bitmap[i] = 0;

  const uint n = groupby.size();//number of dimensions
  StringConverter<T> sc;

  // Estimating simple functions on the union of data streams by Gibbons and Tirthapura
  TupleFactHash<T> cph(n,L,rng);
  string line; //Read a line of data
  vector<T> attribute_value(n); 
  vector<string> attribute_list;//Pour stocker les attributs de dimension et leur valeur de chaque faits
  vector<unsigned int> dimindice;//Pour stocker les indices des colonnes correspondant aux attributs du group by

  N=0;
  //Récupérer les indices des colones des dimensions
  getline (in, line);//La première ligne du fichier de données contient les noms des dimensions
  string delimitertostring (1,delimiter);
  split (line, attribute_list,delimitertostring);
  
  for (unsigned int j=0; j < groupby.size(); ++j)
    for (unsigned int i=0; i < attribute_list.size (); ++i)
      if ( groupby[j] == attribute_list[i]) {
        dimindice.push_back (i);
        break;
      }

  sort(dimindice.begin(), dimindice.end());//Ascending sort of the dimension indices 
  assert(dimindice.size() == n);

  char chline[256]; 
  int index;
  while (in.getline(chline,MAXLENGTHOFALINE, delimiter))
  {
    ++N;
    //cout << endl << "premier " << " " << chline << endl;
    uint col(0);
    bool endline(false);
    for(uint k = 0; k < dimindice.size(); ++k) 
    {
      while( col < dimindice[k])
      {
        ++col;
        if (col == attribute_list.size()-1)
        {
          in.getline(chline, MAXLENGTHOFALINE);
          endline = true;
        }
        else 
          in.getline(chline, MAXLENGTHOFALINE,delimiter);
      }
      //cout << endl << "trouve " << k << " " << chline << endl;
      attribute_value[k] = sc.toType(chline);
    }
    if (!endline) in.getline(chline, MAXLENGTHOFALINE);
    
    uint64 hash = cph.slow(attribute_value); // hashing

    index = trailingzeros(hash,L);
    cout << endl << "index " << index << " hash " << hash << endl;
    if (bitmap[index]== 0) 
      bitmap[index]=1;
  }

  for (uint k=0; k <L; ++k)
    cout << bitmap[k] << "\t";
  cout << endl;

  //Left zero in bitmap
  uint i(0);
  uint64 ans(1);
  while ((bitmap[i] ==1) && (i < L))
  {
    i++;
    ans *=2;
  }

  //cout << endl << " i " << i << endl;
  return ans/phi;
}


//Probabilistic loglogcounting
template <class T>
uint64 LogLogCounting (istream & in, vector<string> & groupby, int& N, int L=32, int m=64, const char& delimiter = ',', int rng = rand_rng)
{
  const float pi=3.14;
  const double log2= log ((double) 2);

  vector<uint> bucket(m,0);// initialized to zero by constructor
  const uint kbits = (uint) round(log ((double) m)/ log2); //k-bits to extract from the hashed values
  //cout << kbits << " " << (1 << kbits)<< " " << m << endl;
  assert(m == 1 << kbits);
  //cout << endl << "valeur de k " << kbits << endl;
  assert (kbits < (uint) L);
  //assert (L <= 32);
  //const double alpha = (double) (0.79402- (2*pi*pi + log2 * log2)/(24*m)); //for checking  alpha(64) must be equal to 0.783
  const double alpha = (double) (0.3970- (2*pi*pi + log2 * log2)/(48*m)); //for checking  alpha(64) must be equal to 0.783
  const uint n = groupby.size();//number of dimensions
  StringConverter<T> sc;
  /*for (int j=0; j < m; ++j)
    bucket[j]=0;*/

  // Estimating simple functions on the union of data streams by Gibbons and Tirthapura
  TupleFactHash<T> cph(n,L,rng);
  string line; //Read a line of data
  vector<T> attribute_value(n); 
  vector<string> attribute_list;//Pour stocker les attributs de dimension et leur valeur de chaque faits
  vector<uint> dimindice;//Pour stocker les indices des colonnes correspondant aux attributs du group by

  N=0;
  //Récupérer les indices des colones des dimensions
  getline (in, line);//La première ligne du fichier de donnès contient les noms des dimensions
  string delimitertostring (1,delimiter);
  split (line,attribute_list,delimitertostring);
  
  for (uint j=0; j < groupby.size(); ++j)
    for (uint i=0; i < attribute_list.size (); ++i)
      if ( groupby[j] == attribute_list[i]) {
        dimindice.push_back (i);
        break;
      }

  sort(dimindice.begin(), dimindice.end());//Ascending sort of the dimension indices 
  assert(dimindice.size() == n);

  char chline[256]; 
  while (in.getline(chline,MAXLENGTHOFALINE,delimiter))
  {
    ++N;
    uint col(0);
    bool endline(false);
    for(uint k = 0; k < dimindice.size(); ++k) 
    {
      while( col < dimindice[k])
      {
        ++col;
        if (col == attribute_list.size()-1)
        {
          in.getline(chline, MAXLENGTHOFALINE);
          endline = true;
        }
        else 
          in.getline(chline, MAXLENGTHOFALINE,delimiter);
      }
      //cout << endl << "trouve " << k << " " << chline << endl;
      attribute_value[k] = sc.toType(chline);
    }
    if (!endline) in.getline(chline, MAXLENGTHOFALINE);
    
    const uint64 hash = cph.slow(attribute_value); // hashing
    //cout << "hash " << hash << " " << k << endl;
    const uint firstkbits = hash >> (uint)(L-kbits); //Extracting first k-bits b0b1...bk-1
    //cout << "kbits " << firstkbits << " " << k << endl;
    uint leftposition = lefttrailingzeros((hash << kbits),L) + 1;//the position of the left 1-bit in bk... (rank from 1)
    bucket[firstkbits] = max (bucket[firstkbits],leftposition);
  }
  double s(0);
  for (uint j=0; j < bucket.size(); ++j)
  {
    s += bucket[j];
    //cout << " " << j << " " << bucket[j] ;
  }
	//cout << "gamma "<<  -1.0/m << " is " << gamma(-1.0/m) << endl;
	/*cout << "gamma "<<  5.0/2 << " is " << exp(gamma((double) (5.0/2))) << endl;
	cout << "gamma "<<  10 << " is " << exp(gamma(10)) << endl;
	//cout << "gamma "<<  -3.0/2 << " is " << gamma((double)(-3/2)) << endl;
	double xg=-1.5;
	cout << "gamma xg " << xg << " is " << (double) ( pi / (double)( exp(gamma(1-xg)) * sin(pi*xg))) << endl;
	xg=-1.0/m;
	cout << "gamma -1/m " << xg << " is " << (double) ( pi / (double)( exp(gamma(1-xg)) * sin(pi*xg))) << endl;
	cout << "Flajolet  -1/m" << xg << " is " << -1*m*gamma((double)(1-xg)) << endl;
	cout << "gamma "<<  2 << " is " << exp(gamma(2)) << endl;
	cout << "gamma "<<  3 << " is " << exp(gamma(3)) << " -- "<< log ((double) 2) << endl;
	//cout << "gamma "<<  << " is " << exp(gamma(3)) << " -- "<< log ((double) 2) << endl;*/
	
 
	
	int ans=0;
	if (m >= 64){
		ans = (int)(alpha * m * pow (2, (double) (s/m)));
	}else{
		//cout << "power part " << (pow(2,1.0/m)-1)/log2<<endl;
		//cout << "gamma part " << -1*pi/(exp(gamma(1+(1.0/m)))*sin(pi/m)) << endl;
		//cout << "product part "<< (-1*pi/(exp(gamma(1+(1.0/m)))*sin(pi/m)))*((pow(2,1.0/m)-1)/log2)<<endl;
		//cout << "alpham "<< pow ((-1*pi/(exp(gamma(1+(1.0/m)))*sin(pi/m)))*((pow(2,1.0/m)-1)/log2),-1*m) <<endl;
		//cout << "gamma " << gamma(1+(1.0/m)) << endl;
		double alpham = pow ((-1*pi/(exp(lgamma(1+(1.0/m)))*sin(pi/m)))*((pow(2,1.0/m)-1)/log2),-1*m);
		ans=(int) (alpham * m * pow (2, (double) (s/m)));
	}
  return ans;
}


//This is the stochastic probalistic counting algorithm
//Only one hash function


void display(vector<bool> & vec) {
  //cout << " hi " << endl;
  for(uint k = 0; k < vec.size() ; ++k) 
    cout << (int) vec[k] << " ";
  cout << endl;
}

template <class T>
uint64 StoProbabilisticCounting (istream & in, vector<string> & groupby, int& N, int L=19, int nbofbitmap=64, const char& delimiter = ',', int rng = rand_rng)
{
  const float phi = 0.77351; //The magic constant
  vector <vector<bool> > bitmap(nbofbitmap);

  for (int j=0; j < nbofbitmap; ++j)
    for (int i=0; i < L; ++i)
      bitmap[j].push_back(0);

  const uint n = groupby.size();//number of dimensions
  StringConverter<T> sc;

  // Estimating simple functions on the union of data streams by Gibbons and Tirthapura
  TupleFactHash<T> cph(n,L,rng);
  string line; //Read a line of data
  vector<T> attribute_value(n); 
  vector<string> attribute_list;//Pour stocker les attributs de dimension et leur valeur de chaque faits
  vector<uint> dimindice;//Pour stocker les indices des colonnes correspondant aux attributs du group by

  N=0;
  //Récupérer les indices des colones des dimensions
  getline (in, line);//La première ligne du fichier de donnès contient les noms des dimensions
  string delimitertostring (1,delimiter);
  split (line,attribute_list,delimitertostring);
  
  for (uint j=0; j < groupby.size(); ++j)
    for (uint i=0; i < attribute_list.size (); ++i)
      if ( groupby[j] == attribute_list[i]) {
        dimindice.push_back (i);
        break;
      }

  sort(dimindice.begin(), dimindice.end());//Ascending sort of the dimension indices 
  assert(dimindice.size() == n);

  char chline[256];
  uint index;
  uint alpha;
  while (in.getline(chline,MAXLENGTHOFALINE,delimiter))
  {
    ++N;
    uint col(0);
    bool endline(false);
    for(uint k = 0; k < dimindice.size(); ++k) 
    {
      while( col < dimindice[k])
      {
        ++col;
        if (col == attribute_list.size()-1)
        {
          in.getline(chline, MAXLENGTHOFALINE);
          endline = true;
        }
        else 
          in.getline(chline, MAXLENGTHOFALINE, delimiter);
      }
      attribute_value[k] = sc.toType(chline);
    }
    if (!endline) in.getline(chline, MAXLENGTHOFALINE);
    
    uint64 hash = cph.slow(attribute_value); // hashing
    alpha = hash % nbofbitmap;
    index = trailingzeros((hash / nbofbitmap),L);
//cout << "alpha = "<< alpha << " projection = " << (hash / nbofbitmap) << endl;
    if (bitmap[alpha][index]== 0)  bitmap[alpha][index]=1;
  }
  double A(0);
  
  for (int j=0; j < nbofbitmap; ++j)
  {
    int r(0);
    /*if(stoverbose) {
      for(int x = 0; x < L;++x) 
       cout << bitmap[j][x] << " ";
      cout << endl;
    }*/
    while ( (bitmap[j][r] == 1) && (r < L))
    {
       r++;
    }
    
    /*int leftmost1 = 0;
    while ( r < L )
    {
       r++;
       if(bitmap[j][r] == 1) leftmost1 = r;
    }*/
     //if(stoverbose) cout << leftmost1 << endl;
    A += r;
  }
if(stoverbose) cout << "nbofbitmap"  << nbofbitmap <<endl;
if(stoverbose) cout << "A"  << A<<endl;  
  int ans = (uint) ((nbofbitmap/phi) * ((pow(2, (A / nbofbitmap)))));
  if(nbofbitmap <= 32){
	  ans = (uint) (ans/(1+0.31/nbofbitmap));
  }
  return ans;
}

//To estimate the time needed to read all the facts
template <class T>
int ReadFactsFromFile (istream & in, vector<string> & groupby, int& N, const char& delimiter = ',')
{
  const uint n = groupby.size();//number of dimensions
  StringConverter<T> sc;
  string line; //Read a line of data
  vector<T> attribute_value(n); 
  vector<string> attribute_list;//Pour stocker les attributs de dimension et leur valeur de chaque faits
  vector<unsigned int> dimindice;//Pour stocker les indices des colonnes correspondant aux attributs du group by

  N=0;

  //Récupérer les indices des colones des dimensions
  getline (in, line);//La première ligne du fichier de données contient les noms des dimensions
  string delimitertostring (1,delimiter);
  split (line,attribute_list,delimitertostring);
  
  for (unsigned int j=0; j < groupby.size(); ++j)
    for (unsigned int i=0; i < attribute_list.size (); ++i)
      if ( groupby[j] == attribute_list[i]) {
        dimindice.push_back (i);
        break;
      }

  sort(dimindice.begin(), dimindice.end());//Ascending sort of the dimension indices 
  assert(dimindice.size() == n);
  char chline[256]; 
  while (in.getline(chline,MAXLENGTHOFALINE,delimiter))
  {
    ++N;
    uint col(0);
    bool endline(false);
    for(uint k = 0; k < dimindice.size(); ++k) 
    {
      while( col < dimindice[k])  
      {
        ++col;
        if (col == attribute_list.size()-1)
        {
          in.getline(chline, MAXLENGTHOFALINE);
          endline = true;
        }
        else 
          in.getline(chline, MAXLENGTHOFALINE, delimiter);
      }
      attribute_value[k] = sc.toType(chline);
      //cout << endl << "trouve " << k << " " << chline << endl;
    }
    if (!endline) in.getline(chline, MAXLENGTHOFALINE);
  }
  return N;
}

//To compute the time needed to read  a data file and hash all its facts
template <class T>
uint64 ReadHashFact(istream & in, vector<string> & groupby, int& N, int L=19, char& delimiter = ',', int rng = rand_rng) 
{
  //map< vector<T>, uint64> buffer;
  const uint n = groupby.size();//number of dimensions
  StringConverter<T> sc;

  // Estimating simple functions on the union of data streams by Gibbons and Tirthapura
  TupleFactHash<T> cph(n,L,rng);
  //uint t = 0;
  //uint mask = (1 << (t)) - 1;   // which is initially 0
  string line; //Read a line of data
  vector<T> attribute_value(n); 
  vector<string> attribute_list;//Pour stocker les attributs de dimension et leur valeur de chaque faits
  vector<unsigned int> dimindice;//Pour stocker les indices des colonnes correspondant aux attributs du group by

  N=0;
  //assert(buffer.size() == 0);  // if not, we needed to clear...

  //Récupérer les indices des colones des dimensions
  getline (in, line);//La première ligne du fichier de données contient les noms des dimensions

  string delimitertostring (1,delimiter);
  split (line,attribute_list,delimitertostring);
  
  for (unsigned int j=0; j < groupby.size(); ++j)
    for (unsigned int i=0; i < attribute_list.size (); ++i)
      if ( groupby[j] == attribute_list[i]) {
        dimindice.push_back (i);
        break;
      }
  sort(dimindice.begin(), dimindice.end());//Ascending sort of the dimension indices 

  assert(dimindice.size() == n);
  char chline[256];
  while (in.getline(chline,MAXLENGTHOFALINE,delimiter))
  {
    ++N;

    //cout << endl << "premier " << " " << chline << endl;
    uint col(0);
    bool endline(false);
    for(uint k = 0; k < dimindice.size(); ++k) 
    {
      while( col < dimindice[k])  
      {
        ++col;
        if (col == attribute_list.size()-1)
        {
          in.getline(chline, MAXLENGTHOFALINE);
          endline = true;
        }
        else 
          in.getline(chline, MAXLENGTHOFALINE,delimiter);
      }
      //cout << endl << "trouve " << k << " " << chline << endl;
      attribute_value[k] = sc.toType(chline);
    }
    if (!endline) in.getline(chline, MAXLENGTHOFALINE);

    cph.slow(attribute_value); // hashing
 }
    //buffer[attribute_value] = (1LL << 32) | trailingzeros(oldhash,L);  // count= 1 
  return N;
}


//Adaptive counting
//by Min Cai et al.

template <class T>
uint64 AdaptiveCounting (istream & in, vector<string> & groupby, int& N, int L=32, int m=64, const char& delimiter = ',', int rng = rand_rng)
{
  const float pi=3.14;
  const double log2= log ((double) 2);

  vector<uint> bucket(m,0);// initialized to zero by constructor
  const uint kbits = (uint) round(log ((double) m)/ log2); //k-bits to extract from the hashed values
  assert(m == 1 << kbits);
  assert (kbits < (uint) L);
  const double alpha = (double) (0.3970- (2*pi*pi + log2 * log2)/(48*m)); //alpha(64) must be equal to 0.783
  const uint n = groupby.size();//number of dimensions
  StringConverter<T> sc;
  TupleFactHash<T> cph(n,L,rng);
  string line; //Read a data line
  vector<T> attribute_value(n); 
  vector<string> attribute_list;//Pour stocker les attributs de dimension et leur valeur de chaque faits
  vector<uint> dimindice;//Pour stocker les indices des colonnes correspondant aux attributs du group by

  N=0;
  //get dimesions' indices
  getline (in, line);//The firt line contains the name of each dimensions
  string delimitertostring (1,delimiter);
  split (line,attribute_list,delimitertostring);
  
  for (uint j=0; j < groupby.size(); ++j)
    for (uint i=0; i < attribute_list.size (); ++i)
      if ( groupby[j] == attribute_list[i]) {
        dimindice.push_back (i);
        break;
      }

  sort(dimindice.begin(), dimindice.end());//Ascending sort of the dimension indices 
  assert(dimindice.size() == n);

  char chline[256]; 
  while (in.getline(chline,MAXLENGTHOFALINE,delimiter))
  {
    ++N;
    uint col(0);
    bool endline(false);
    for(uint k = 0; k < dimindice.size(); ++k) 
    {
      while( col < dimindice[k])
      {
        ++col;
        if (col == attribute_list.size()-1)
        {
          in.getline(chline, MAXLENGTHOFALINE);
          endline = true;
        }
        else 
          in.getline(chline, MAXLENGTHOFALINE,delimiter);
      }
      //cout << endl << "trouve " << k << " " << chline << endl;
      attribute_value[k] = sc.toType(chline);
    }
    if (!endline) in.getline(chline, MAXLENGTHOFALINE);
    
    const uint64 hash = cph.slow(attribute_value); // hashing
    const uint firstkbits = hash >> (uint)(L-kbits); //Extracting first k-bits b0b1...bk-1
    uint leftposition = lefttrailingzeros((hash << kbits),L) + 1;//the position of the left 1-bit in bk... (rank from 1)
    bucket[firstkbits] = max (bucket[firstkbits],leftposition);
  }
  double s(0);
  int empty(0);
  for (uint j=0; j < bucket.size(); ++j)
  {
    s += bucket[j];
    if (bucket[j]==0) empty++;
  }
	
	int ans=0;
	if (((double)empty/m) >= 0.051)
	{
		ans = (int) (-m *log(((double)empty/m)));
		//cout << "empty double " << ((double)empty/m) << endl;
	}else{
		if (m >= 64){
			ans = (int)(alpha * m * pow (2, (double) (s/m)));
		}else{
			double alpham = pow ((-1*pi/(exp(lgamma(1+(1.0/m)))*sin(pi/m)))*((pow(2,1.0/m)-1)/log2),-1*m);
			ans=(int) (alpham * m * pow (2, (double) (s/m)));
		}
	}
  return ans;
}

#endif
