#include "estimator.h"


void unittesting(int L = 55) {
    RandomHasher<string> rh(L);
    int times = 10000;
    int counter = 0;
    for(int k = 0; k < times ; ++k) {
        if(rh.generateRandomInteger()>pow(2,L-1))
          ++counter; 
    }
    cout << "percentage of random numbers above average = " << ((double) counter / times) << endl;
}

int main(int params, char ** args) 
{
	enum {gibbons,counting,readfile,loglog,readhash,stogibbons,stosumgibbons,stoweightedsumgibbons,baryossef,adaptive};  // methods
	int L=19;
	int s=256;
	int r=time(0);//4567;  // rand seed
	int stogibbonscount = -1;
	int rng = rand_rng;
	char delimiter = ',';
	vector<string> groupby; //Group by attributes
	bool eightbits = false;
	bool thirtytwobits = false;
	bool verbose = false;
	int method = gibbons;
  	bool richeroutput = false;

	if (params <= 3) {
		cout << "You must provide some parameters including a groupby and a filename " << endl;
		return -1;
	}
	
	string datafilename (args[params-1]); //Data file name
	ifstream myfilein(datafilename.c_str()); //Opening
  	istream &filein =   (datafilename == "stdin") ? cin : myfilein;

	//Reading and checking the program parareters
	for(int i = 1; i < params - 1; ++i)
	{
		if(strcmp(args[i],"-v") == 0)
		{
			verbose = true;
		}
		if(strcmp(args[i],"--rich") == 0)
		{
			richeroutput = true;
		}
		if(strcmp(args[i],"--groupby") == 0)
		{
			string sgroupby;//Group by attributes as a string
			
			if( params - i > 1 ) {
				sgroupby = args[++i];
				split(sgroupby, groupby, ",");
				if (verbose==true) cout << "Group by attributes: ";
				if (verbose==true)
				{
					for (uint l=0; l < groupby.size(); ++l)
					cout << groupby[l] << " ";
					cout << endl;
				}
			} else
			{
				cerr << "--groupby expects a string" << endl;//Voir comment verifier une cha�e en c++
				return -2;
			}
		}
		if(strcmp(args[i],"-d") == 0)
		{
			string stringdelim = args[++i];
			delimiter = stringdelim.c_str()[0];
			if (verbose==true) cout << "Delimiter of group by values: " << delimiter <<endl;
		}
		if(strcmp(args[i],"-s") == 0)
		{
			s = atoi(args[++i]);
			if (verbose==true) cout << "Memory size or number of bitmaps: " << s <<endl;
		}

		if(strcmp(args[i],"-l") == 0)
		{
			L = atoi(args[++i]);
			if (verbose==true) cout << "Number of bits for hashing: " << L <<endl;
		}

		if(strcmp(args[i],"--8bits") == 0)
		{
			eightbits = true;
			if (verbose==true) cout << "Assuming attribute values fit in 8 bits" << endl;
		}
		if(strcmp(args[i],"--32bits") == 0)
		{
			thirtytwobits = true;
			if (verbose==true) cout << "Assuming attribute values fit in 32 bits" << endl;
		}
		if(strcmp(args[i],"--gibbons")==0)
		{
			method = gibbons;
			if (verbose==true) cout << "Counting by the Gibbons' method" << endl;
		}
		if(strcmp(args[i],"--baryossef")==0)
		{
			method = baryossef;
			if (verbose==true) cout << "Counting by the Bar-Yossef' method" << endl;
		}
		if(strcmp(args[i],"--adaptive")==0)
		{
			method = adaptive;
			if (verbose==true) cout << "Counting by the Adaptive counting' method" << endl;
		}
		if(strcmp(args[i],"--stogibbons")==0)
		{
			method = stogibbons;
			stogibbonscount = atoi(args[++i]);
			if (verbose==true) cout << "Counting by the Sto. Gibbons' method" << endl;
		}
		if(strcmp(args[i],"--stosumgibbons")==0)
		{
			method = stosumgibbons;
			stogibbonscount = atoi(args[++i]);
			if (verbose==true) cout << "Counting by the Sto. Sum Gibbons' method" << endl;
		}
		if(strcmp(args[i],"--stoweightedsumgibbons")==0)
		{
			method = stoweightedsumgibbons;
			stogibbonscount = atoi(args[++i]);
			if (verbose==true) cout << "Counting by the Sto. Weighted Sum Gibbons' method" << endl;
		}
		if(strcmp(args[i],"--counting")==0)
		{
			method = counting;
			if (verbose==true) cout << "Counting by the probalistic counting" << endl;
		}
		if(strcmp(args[i],"-r") == 0 )
                {
			r= atoi(args[++i]);
			if (verbose==true) cout << "Random seed: " << r <<endl;
		}
		if(strcmp(args[i],"--loglog")==0)
		{
			method = loglog;
			if (verbose==true) cout << "Counting by the loglog counting" << endl;
		}

		if(strcmp(args[i],"--readfile")==0)
		{
			method = readfile;
			if (verbose==true) cout << "Just reading..." << endl;
		}
		if(strcmp(args[i],"--readhash")==0)
		{
			method = readhash;
			if (verbose==true) cout << "Just reading and hashing..." << endl;
		}
	}
	// seed random number generators
	srand(r); 
	// extra rngs seeded on demand
	randomSource.seed(r);
	//cout << rand() << " " << rand() << endl;
	if (!filein) {
		cerr << "Can't open the damned file. " <<endl;
		return -10;
	}
	
	uint64 answer = 0;  // the count-distinct estimate
	int N(0);
	
	clock_t start,finish;
	start = clock();
	switch (method){
	case gibbons:
		if(eightbits)
			answer = GibbonsTirthapura<char>(filein,groupby,N,L,s,delimiter,rng);
		else if(thirtytwobits)
			answer = GibbonsTirthapura<int>(filein,groupby,N,L,s,delimiter,rng);
			else {
				if (verbose==true) cout << "Assuming attribute values are strings???" << endl;
				if (verbose==true) cout << "Please try --8bits or --32bits instead, this is going to be slow"<<endl;
				answer = GibbonsTirthapura<string>(filein, groupby,N,L,s,delimiter,rng);
			}
			break;
			
	case stogibbons:
		  if(eightbits)
			  answer = StoGibbonsTirthapura<char>(filein,groupby,N,L,s,stogibbonscount,delimiter,rng);
		  else if(thirtytwobits)
			  answer = StoGibbonsTirthapura<int>(filein,groupby,N,L,s,stogibbonscount,delimiter,rng);
		  else {
			  if (verbose==true) cout << "Assuming attribute values are strings???" << endl;
			  if (verbose==true) cout << "Please try --8bits or --32bits instead, this is going to be slow"<<endl;
			  answer = StoGibbonsTirthapura<string>(filein, groupby,N,L,s,stogibbonscount,delimiter,rng);
		  }
		  break;
		  
	case stosumgibbons:
		  if(eightbits)
			  answer = StoSumGibbonsTirthapura<char>(filein,groupby,N,L,s,stogibbonscount,delimiter,rng);
		  else if(thirtytwobits)
			  answer = StoSumGibbonsTirthapura<int>(filein,groupby,N,L,s,stogibbonscount,delimiter,rng);
		  else {
			  if (verbose==true) cout << "Assuming attribute values are strings???" << endl;
			  if (verbose==true) cout << "Please try --8bits or --32bits instead, this is going to be slow"<<endl;
			  answer = StoSumGibbonsTirthapura<string>(filein, groupby,N,L,s,stogibbonscount,delimiter,rng);
		  }
		  break;
	case stoweightedsumgibbons:
		if(eightbits)
			answer = StoWeightedSumGibbonsTirthapura<char>(filein,groupby,N,L,s,stogibbonscount,delimiter,rng);
		else if(thirtytwobits)
			answer = StoWeightedSumGibbonsTirthapura<int>(filein,groupby,N,L,s,stogibbonscount,delimiter,rng);
			else {
				if (verbose==true) cout << "Assuming attribute values are strings???" << endl;
				if (verbose==true) cout << "Please try --8bits or --32bits instead, this is going to be slow"<<endl;
				answer = StoWeightedSumGibbonsTirthapura<string>(filein, groupby,N,L,s,stogibbonscount,delimiter,rng);
			}
			break;

	case counting:
		  if(eightbits)
			  answer = StoProbabilisticCounting<char>(filein,groupby,N,L,s,delimiter,rng);
		  else if(thirtytwobits)
			  answer = StoProbabilisticCounting<int>(filein,groupby,N,L,s,delimiter,rng);
		  else {
			  if (verbose==true) cout << "Assuming attribute values are strings???" << endl;
			  if (verbose==true) cout << "Please try --8bits or --32bits instead, this is going to be slow"<<endl;
			  answer = StoProbabilisticCounting<string>(filein,groupby,N,L,s,delimiter,rng);
		  }
		  break;
     
	case baryossef:
		if(eightbits)
			answer = BarYossef<char>(filein,groupby,N,L,s,delimiter,rng);
		else if(thirtytwobits)
			answer = BarYossef<int>(filein,groupby,N,L,s,delimiter,rng);
			else {
				if (verbose==true) cout << "Assuming attribute values are strings???" << endl;
				if (verbose==true) cout << "Please try --8bits or --32bits instead, this is going to be slow"<<endl;
				answer = BarYossef<string>(filein,groupby,N,L,s,delimiter,rng);
			}
			break;
		
	case loglog:
		  if(eightbits)
			  answer = LogLogCounting<char>(filein,groupby,N,L,s,delimiter,rng);
		  else if(thirtytwobits)
			  answer = LogLogCounting<int>(filein,groupby,N,L,s,delimiter,rng);
		  else {
			  if (verbose==true) cout << "Assuming attribute values are strings???" << endl;
			  if (verbose==true) cout << "Please try --8bits or --32bits instead, this is going to be slow"<<endl;
			  answer = LogLogCounting<string>(filein,groupby,N,L,s,delimiter,rng);
		  }
		  break;
		  
	case adaptive:
		  if(eightbits)
			  answer = AdaptiveCounting<char>(filein,groupby,N,L,s,delimiter,rng);
		  else if(thirtytwobits)
			  answer = AdaptiveCounting<int>(filein,groupby,N,L,s,delimiter,rng);
		  else {
			  if (verbose==true) cout << "Assuming attribute values are strings???" << endl;
			  if (verbose==true) cout << "Please try --8bits or --32bits instead, this is going to be slow"<<endl;
			  answer = AdaptiveCounting<string>(filein,groupby,N,L,s,delimiter,rng);
		  }
		  break;
	
	case readfile:
		  if(eightbits)
			  answer = ReadFactsFromFile<char>(filein,groupby,N,delimiter);
		  else if(thirtytwobits)
			  answer = ReadFactsFromFile<int>(filein,groupby,N,delimiter);
		  else {
			  if (verbose==true) cout << "Assuming attribute values are strings???" << endl;
			  if (verbose==true) cout << "Please try --8bits or --32bits instead, this is going to be slow"<<endl;
			  answer = ReadFactsFromFile<string>(filein,groupby,N,delimiter);
		  }
		  break;
		
	 case readhash:
		  if(eightbits)
			  answer = ReadHashFact<char>(filein,groupby,N,L,delimiter,rng);
		  else if(thirtytwobits)
			  answer = ReadHashFact<int>(filein,groupby,N,L,delimiter,rng);
		  else {
			  if (verbose==true) cout << "Assuming attribute values are strings???" << endl;
			  if (verbose==true) cout << "Please try --8bits or --32bits instead, this is going to be slow"<<endl;
			  answer = ReadHashFact<string>(filein,groupby,N,L,delimiter,rng);
		  }
		  break;
	  
	default : cerr << "Unknown method " << method << endl;
	
	}
	finish = clock();
	if (verbose==true)
	{
		cout << "Estimated view size: " << answer << endl;
		cout << "Time used to read data: " << (double)(finish - start) / CLOCKS_PER_SEC << " seconds or " << (finish - start)<< " ticks "<<endl;
		cout << "Nb of tuples: " << N << endl;
	} if (richeroutput == true) {
    cout << answer<< " " << N << " " << (double)(finish - start) / CLOCKS_PER_SEC << endl;
  } else
	{
		//cout << s << "\t" << r << "\t" << answer << "\t" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
		cout << answer << " " << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	}
  myfilein.close();
	return 0;
}
