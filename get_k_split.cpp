#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib> 
#include <string> 
#include <algorithm>			    					// std::sort
#include <math.h>
#include <iomanip>
#include <ctime>
#include <time.h>
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/gamma_distribution.hpp"
#include "boost/random/normal_distribution.hpp"
#include <boost/math/distributions/beta.hpp>
#include <boost/random.hpp>
#include <boost/math/special_functions/digamma.hpp> 				//needed for computing the modified Dirichlet bound
#include <boost/math/special_functions/gamma.hpp>
#include <boost/foreach.hpp>
#include "omp.h"
#include "infiles.h"

//#include "boost/program_options.hpp" 
//#include <iterator>

using namespace std;
using namespace boost::math::policies;
using namespace boost::math;



int read()
{

//cout<<prob_file<<endl;
	ifstream inFile;
	string line;
	int K, i;

	inFile.open(vb_file, ifstream::in);
	if (!inFile) {
    		cerr << "Unable to open "<< vb_file<<" file"<<endl;
    		exit(1);   // call system to stop
	}
	//cout<< "Reading "<< vb_file<<" file"<<endl;

	//getline( inFile, line );
	inFile.ignore(3);	
	inFile >> i;
	inFile >> K;
	K++;
	K++;
	cout<< "number of components = "<< K<<endl;
	inFile.close();

return(K);
}

