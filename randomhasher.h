#ifndef RANDOMHASHER
#define RANDOMHASHER
#include <vector>
#include <deque>
#include <map>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <assert.h>
#include <ext/hash_map>
#include <ext/hash_fun.h>
#include "mersennetwister.h"
extern "C" {
#include "randomlib.h"
}
using namespace std;
using __gnu_cxx::hash_map;  // gcc and SGI extension.
using __gnu_cxx::hash;  // gcc and SGI extension.
// some version to be blessed as unordered_map, someday.
using namespace std;
//size_t k = hash< const char*>("fsfd");
namespace  __gnu_cxx                                                                                
{                                                                                             
  template<> struct hash<string>                                                       
  {                                                                                           
    size_t operator()( const std::string& x ) const                                           
    {                                                                                         
      return hash< const char* >()( x.c_str() );                                              
    }                                                                                         
  };       

} 


//Type definition
typedef unsigned int uint;
typedef unsigned long long uint64;
typedef unsigned long uint32;


// kinds of random number generators known here
enum {rand_rng=1, mt_rng=2, noisefile_rng=3, mars_rng=4};

//const bool verboseusemt(false);

class mersenneRNG 
{
public:
	mersenneRNG(int nn) : mtr(25),n(nn) {};
	uint32 operator()() 
	{ 
		uint32 ans = mtr.randInt(n);
		return ans;
	} 
	void seed(uint32 seedval) 
	{ 
		mtr.seed(seedval);
	}
	void seed() 
	{ 
		mtr.seed();
	}
	uint32 rand_max() 
	{ 
		return n;
	}

private:
	MTRand mtr;
	int n;
};

mersenneRNG randomSource(0xffffffff);  // global var


// use a file from random.org
class noisefileRNG 
{
public:
	noisefileRNG(uint32 nn) : noisefile("10megs.001+"), n(nn)  {};
	uint32 operator()() 
	{
		uint32 ans(0);
		for (uint i=0; i < sizeof(ans); ++i) 
		{
			if (noisefile.eof()) cerr << "out of random bytes" << endl;
			unsigned char somebits = (unsigned char) noisefile.get();
			ans = (ans << 8) | somebits;
		}
		ans %= n;
		return ans;
	}
	// seed by seeking in file: far, but not too far.  Stay in first 2MB. But amplify seeds
	void seed(uint32 seedval) 
	{   
		int smallseed = (seedval * 9876543) % 2000000;
		noisefile.seekg(smallseed, ios_base::beg);
	}
	void seed() 
	{ 
		seed(0);
	}
	uint32 rand_max() 
	{ 
		return n;
	}
private:
	ifstream noisefile;
	int n;
};


noisefileRNG randomSourceFile(0xffffffff);  // global var


// Marsaglia and Zaman RNG
// note: there can be only one instance of this random number generator
/**class mzRNG 
{
public:
	mzRNG(uint32 nn) : n(nn) 
	{
		if (instanceCtr++ != 0) 
			cerr << "more than one instance of MZ RNG is a problem!" << endl;
		if (nn > 0x7fffffff) 
			cerr << "warning in mz: truncating randmax" << endl;
		//seed();
	};
  uint32 operator()() 
  { 
	  uint32 ans = RandomInt(0,n);
	  // cout << "gave Marsaglia " << n << "got back " << ans << endl;
	  return ans;
  } 
  void seed(uint32 seedval) 
  { 
	  unsigned short part1 = (seedval % 31329);
	  unsigned short part2 = (seedval >> 15) % 30082;
	  RandomInitialise(part1,part2);
  }
  void seed() 
  { 
	  seed(25);
  }
  uint32 rand_max() 
  { 
	  return n;
  }
private:
	uint32 n;
	static int instanceCtr;
};

int mzRNG::instanceCtr = 0;
mzRNG randomSourceMZ(0x7fffffff);  // global var
**/
/**
* This simply hash objects to 64 (by TupleFactHash below!!!
*/

template <class T>
class RandomHasher {
public:
	RandomHasher(): wordsize(-1), mMap(), logrm(-1), rngKind(rand_rng) {}
	
	RandomHasher(int mywordsize, int rng = rand_rng): wordsize(mywordsize), mMap(),logrm(0),rngKind(rng)
	{
		cout << "[WARNING] USING SLOW RANDOM HASHER, consider switching to integers."<< endl<< endl<< endl;
    unsigned int rm;
		switch (rngKind)
		{
			case rand_rng: rm = RAND_MAX; break;
			case mt_rng: rm = randomSource.rand_max(); break;
			//case mars_rng: rm = randomSourceMZ.rand_max(); break;
			case noisefile_rng: rm = randomSourceFile.rand_max(); break;
			default: cerr << "rngKind = " << rngKind << " unknown" << endl; rm = 57;
		}
		
		while(rm > 0) 
		{
			rm = rm >> 1;
			++logrm; // Nombre de bits décalés log rm
		}
	}
	
	uint64 hash(T t) 
	{
		typename map<T,uint64>::iterator p = mMap.find(t);
		if( p != mMap.end())
			return p->second;
		return mMap[t] = generateRandomInteger();
	}
	
	void forcehash(T t, uint64 hashvalue) 
	{
		mMap[t] = hashvalue;
	}
	
	uint64 generateRandomInteger() 
	{
		int bits = wordsize;
		uint64 answer(0);
		while(bits != 0) 
		{
			const int offset = bits - logrm;
			uint64 randombits = 0;
			switch(rngKind) 
			{
				case rand_rng: randombits = rand(); break;
				case mt_rng: randombits = randomSource(); break;
				//case mars_rng: randombits = randomSourceMZ(); break;
				case noisefile_rng : randombits = randomSourceFile(); break;
				default: cerr << "rngKind=" << rngKind << " unk in hash" << endl;
			}
			if(offset < 0) 
			{
				answer = answer | (randombits >> (-offset));
				break;
			} else 
			{
				answer = answer | (randombits << offset);
				bits = bits - logrm;
			}
		}
		return answer;
	}
	
	int getWordSize() const
	{
		return wordsize;
	}

private:
	int wordsize;
	map<T,uint64> mMap;
	int logrm;
	int rngKind;
};


template <>
class RandomHasher<string> {
public:
	RandomHasher<string>(): wordsize(-1), mMap(), logrm(-1), rngKind(rand_rng) {}
	
	RandomHasher<string>(int mywordsize, int rng = rand_rng): wordsize(mywordsize), mMap(),logrm(0),rngKind(rng)
	{
		unsigned int rm;
		switch (rngKind)
		{
			case rand_rng: rm = RAND_MAX; break;
			case mt_rng: rm = randomSource.rand_max(); break;
			//case mars_rng: rm = randomSourceMZ.rand_max(); break;
			case noisefile_rng: rm = randomSourceFile.rand_max(); break;
			default: cerr << "rngKind = " << rngKind << " unknown" << endl; rm = 57;
		}
		//cout << RAND_MAX << endl;
		while(rm > 0) 
		{
			rm = rm >> 1;
			++logrm; // Nombre de bits décalés log rm
		}
		//cout << logrm << endl;
	}
	
	uint64 hash(string t) 
	{
		hash_map<string,uint64>::iterator p = mMap.find(t);
		if( p != mMap.end())
			return p->second;
		return mMap[t] = generateRandomInteger();
	}
	

	
	uint64 generateRandomInteger() 
	{
		int bits = wordsize;
		uint64 answer(0);
		while(bits != 0) 
		{
			const int offset = bits - logrm;
			uint64 randombits = 0;
			switch(rngKind) 
			{
				case rand_rng: randombits = rand(); break;
				case mt_rng: randombits = randomSource(); break;
				//case mars_rng: randombits = randomSourceMZ(); break;
				case noisefile_rng : randombits = randomSourceFile(); break;
				default: cerr << "rngKind=" << rngKind << " unk in hash" << endl;
			}
			if(offset < 0) 
			{
				answer = answer | (randombits >> (-offset));
				break;
			} else 
			{
				answer = answer | (randombits << offset);
				bits = bits - logrm;
			}
		}
		return answer;
	}
	
	int getWordSize() const
	{
		return wordsize;
	}

private:
	int wordsize;
	hash_map<string,uint64> mMap;
	int logrm;
	int rngKind;
};


template <>
class RandomHasher<int> {
public:
	RandomHasher<int>(): wordsize(-1), mMap(), logrm(-1), rngKind(rand_rng) {}
	
	RandomHasher<int>(int mywordsize, int rng = rand_rng): wordsize(mywordsize), mMap(),logrm(0),rngKind(rng)
	{
		unsigned int rm;
		switch (rngKind)
		{
			case rand_rng: rm = RAND_MAX; break;
			case mt_rng: rm = randomSource.rand_max(); break;
			//case mars_rng: rm = randomSourceMZ.rand_max(); break;
			case noisefile_rng: rm = randomSourceFile.rand_max(); break;
			default: cerr << "rngKind = " << rngKind << " unknown" << endl; rm = 57;
		}
		
		while(rm > 0) 
		{
			rm = rm >> 1;
			++logrm; // Nombre de bits décalés log rm
		}
	}
	
	uint64 hash(int t) 
	{
		hash_map<int,uint64>::iterator p = mMap.find(t);
		if( p != mMap.end())
			return p->second;
		return mMap[t] = generateRandomInteger();
	}
	

	
	uint64 generateRandomInteger() 
	{
		int bits = wordsize;
		uint64 answer(0);
		while(bits != 0) 
		{
			const int offset = bits - logrm;
			uint64 randombits = 0;
			switch(rngKind) 
			{
				case rand_rng: randombits = rand(); break;
				case mt_rng: randombits = randomSource(); break;
				//case mars_rng: randombits = randomSourceMZ(); break;
				case noisefile_rng : randombits = randomSourceFile(); break;
				default: cerr << "rngKind=" << rngKind << " unk in hash" << endl;
			}
			if(offset < 0) 
			{
				answer = answer | (randombits >> (-offset));
				break;
			} else 
			{
				answer = answer | (randombits << offset);
				bits = bits - logrm;
			}
		}
		return answer;
	}
	
	int getWordSize() const
	{
		return wordsize;
	}

private:
	int wordsize;
	hash_map<int,uint64> mMap;
	int logrm;
	int rngKind;
};

template <>
class RandomHasher<char> {
public:
	RandomHasher<char>(): wordsize(-1), mMap(), logrm(-1), rngKind(rand_rng) {}
	
	RandomHasher<char>(int mywordsize, int rng = rand_rng): wordsize(mywordsize), mMap(),logrm(0),rngKind(rng)
	{
		unsigned int rm;
		switch (rngKind)
		{
			case rand_rng: rm = RAND_MAX; break;
			case mt_rng: rm = randomSource.rand_max(); break;
			//case mars_rng: rm = randomSourceMZ.rand_max(); break;
			case noisefile_rng: rm = randomSourceFile.rand_max(); break;
			default: cerr << "rngKind = " << rngKind << " unknown" << endl; rm = 57;
		}
		
		while(rm > 0) 
		{
			rm = rm >> 1;
			++logrm; // Nombre de bits décalés log rm
		}
    for(uint i = 0; i < 256; ++i) 
      mMap[i] = generateRandomInteger();
	}
	
	uint64 hash(char t) 
	{
		return mMap[(unsigned char) t];
	}
	
	
	uint64 generateRandomInteger() 
	{
		int bits = wordsize;
		uint64 answer(0);
		while(bits != 0) 
		{
			const int offset = bits - logrm;
			uint64 randombits = 0;
			switch(rngKind) 
			{
				case rand_rng: randombits = rand(); break;
				case mt_rng: randombits = randomSource(); break;
				//case mars_rng: randombits = randomSourceMZ(); break;
				case noisefile_rng : randombits = randomSourceFile(); break;
				default: cerr << "rngKind=" << rngKind << " unk in hash" << endl;
			}
			if(offset < 0) 
			{
				answer = answer | (randombits >> (-offset));
				break;
			} else 
			{
				answer = answer | (randombits << offset);
				bits = bits - logrm;
			}
		}
		return answer;
	}
	
	int getWordSize() const
	{
		return wordsize;
	}

private:
	int wordsize;
  uint64 mMap[256];
	//map<T,uint64> mMap;
	int logrm;
	int rngKind;
};


/**
* This is is a tuple-wise independent hash (should be very good, but slow).
*/
template <class T>
class TupleFactHash {
public:
	TupleFactHash(int myn, int mywordsize=19, int rng = rand_rng) : n(myn), wordsize(mywordsize), rh(0){
      //if(verboseusemt) cout << "[TupleWiseIndependent] rng = "<< rng << endl;
		 //cout << "[TupleWiseIndependent] rng = "<< rng << endl;
      for(int k = 0; k < n; ++k)
        rh.push_back(RandomHasher<T>(wordsize, rng));
    }
    uint64 slow(const vector<T>& chars)  {
      uint64 k = 0;
      for(int j = 0; j < n; ++j)
        k = k  ^ rh[j].hash(chars[j]);
      return k;     
    }

    uint64 slow(const deque<T>& chars)  {
      uint64 k = 0;
      for(int j = 0; j < n; ++j) 
        k = k  ^ rh[j].hash(chars[j]);
      return k;     
    }

    inline  uint64 semiFast(uint64 oldhash, T inchar, int pos) {
      return oldhash ^ rh[pos].hash(inchar);
    }

    int hashSize() const {
      return wordsize;
    }
    
    string __str__() const {
        stringstream strs;
        strs << "TupleFactHash" <<n <<"_"<<wordsize;
        return strs.str() ;
    }
    const int n, wordsize;

    vector<RandomHasher<T> > getCopyOfHashers() const {return rh;}
    
    void useACopyOfTheseHashers(vector<RandomHasher<T> > & myhasher) {
      for(uint k = 0; k < myhasher.size(); ++k)
        assert(myhasher[k].getWordSize()==rh[k].getWordSize());
      rh=myhasher;
    }
private: 
    vector<RandomHasher<T> > rh;
};


/*template <class T>
class RandomTupleHasher 
{
public:
	RandomTupleHasher(int myn, int mywordsize=19, int rng = rand_rng) : n(myn), wordsize(mywordsize), rh(mywordsize, rng){}

    //  convert to a deque if necessary, I guess
    //uint64 slow(const vector<T>& chars)  {
    //}


    uint64 slow(const deque<T>& chars)  
	{
		return rh.hash(chars);     // I think that == for deques will do the right thing.
                                 // Stroustrup says they're like vectors and vectors
                                 // have == overloaded in a nice way.
	}

    int hashSize() const 
	{
      return wordsize;
    }
    
    string __str__() const {
        stringstream strs;
        strs << "CPPRandomTuple" <<n <<"_"<<wordsize;
        return strs.str() ;
    }

    const int n, wordsize;
    
 private: 
    RandomHasher<deque<T> > rh;
};*/
#endif

