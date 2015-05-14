#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib> 
#include <string> 
#include <sstream>
#include <algorithm>			    					// std::sort
#include <math.h>
#include <numeric>
#include <iomanip>
#include <ctime>
#include <time.h>
#include <boost/foreach.hpp>
#include "omp.h"
#include "get_k_split.h" 
#include "infiles.h" 
//#include "gen_dir.h" 
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/gamma_distribution.hpp"
#include "boost/random/normal_distribution.hpp"
#include <boost/math/distributions/beta.hpp>
#include <boost/random.hpp>
#include <boost/math/special_functions/digamma.hpp> 				//needed for computing the modified Dirichlet bound
#include <boost/math/special_functions/gamma.hpp>
//#include <boost/random/uniform_01.hpp>
//#include "omp.h"

using namespace std;
using namespace boost::math::policies;
using namespace boost::math;



//****************************************************************************************************************************************
	//sample A
	vector< vector <double> > x_1;					// double x[n][p]:	pdfvalues
	vector< int > y_1; 							// int y[n]:		number of mappings per read
	vector< vector <int> > c_1; 						// int c_new[n][p]:	relabelled index of mapping components per read.
	//sample B
	vector< vector <double> > x_2; 						// double x[n][p]:	pdfvalues
	vector< int > y_2; 							// int y[n]:		number of mappings per read
	vector< vector <int> > c_2; 						// int c_new[n][p]:	relabelled index of mapping components per read.

	vector< int > r;
	vector< int > s;

	vector< int > z_alloc;							//z = (z_1,...,z_r)
	vector< int > xi_alloc;							//xi = (xi_1,...,xi_s)

	//double * theta;
	//double * w;
	vector<double> theta (K,0.0);
	vector<double> w (K,0.0);	
	vector<double> sum_1 (K,0.0);
	vector<double> sum_2 (K,0.0);

	vector<double> sum_1_new (K,0.0);
	vector<double> sum_2_new (K,0.0);

	vector<int> cum_sum_1 (n_1,0);
	vector<int> cum_sum_2 (n_2,0);
	double p_alloc, p_alloc_new, myStheta, mySw;
	double log_likelihood, log_likelihood_new, gPrior;
//	gPrior denotes the probability of being DE, and it follows: gPrior ~ Beta(0.5,0.5), that is, the Jefreys Prior for a Bernoulli trial
	vector<int> c_vec (K,1);	//c_vector
	vector<int> c_new (K,0);	//c_vector
	vector<int> permutation (K,0);
	vector<int> inv_permutation (K,0);
	vector<int> permutation_new (K,0);
	vector<int> inv_permutation_new (K,0);

	//double * u;
	//double * v;
	vector<double> u (K,0.0);	
	vector<double> v (K,0.0);	

	double u_sum_c1, ty1, ty2;
	double u_sum_c1_new;	

	vector<double> alpha_prior (K,1.0);	//prior for u
	vector<double> gamma_prior (K,1.0);	//prior for v
	vector<double> v_new (K,0.0);
	
	int nMIX = 5, mixID, kRest, myClustID;
	vector<double> beta_paramNEW (nMIX,1.0); //beta parameters of mixture

	const double threshold = pow(10.0,-80.0);
	void gen_dir(vector<double>&, vector<double>&, int);
	void dir(vector<double>&, int);
	double allocate(int , int );
	double allocate_new(int , int );
	double allocateCollapsed(int , int );
	double unifRand();
	void log_like();
	void log_like_new();
	void permutation_current();
	void permutation_current_collapsed();
	void permutation_prop();
	void permutation_prop_collapsed();
	double bp(int);
	double dp(int);
	double  beta_rand( double , double );
	double log_beta(double , double , double );
	double log_f_c_old();
	double log_f_c_new();
	double birth_accept_prob(double , double ,double , double ,int ,int ,double );
	double birth_mh(double , double , int ,int );
//****************************************************************************************************************************************
	int n_obs, c_sum, c_sum_new;								// number of components and number of observations
//****************************************************************************************************************************************

	// temp vectors
	vector<double> z_prob (K,0.0);
	vector<double> alpha_gd(K ,0.0);
	vector<double> beta_gd(K ,0.0);
	vector<double> w_temp (K,0.0);	
	vector<double> u_temp (K,0.0);	
	vector<double> u_temp1 (K,0.0);	

	vector<double> papasA;	
	vector<double> papasB;

int main()
{
	int iterNEW, myBurn;
	myBurn = burn;
	//if(K > 200){
	//	myBurn = 3*burn;
	//}

	beta_paramNEW[0] = 1.0;
	beta_paramNEW[1] = 10.0;
	beta_paramNEW[2] = 100.0;
	beta_paramNEW[3] = 250.0;
	beta_paramNEW[4] = 500.0;



	iterNEW = iterations + myBurn;
	//if(K > 200){
	//	iterNEW = 5*iterations + myBurn;
	//}




	ifstream inFile, trFile;
	int i, j, num, k, h, place = 0, step = 0, iter, nan_number = 0, cluster, tid;
	double temp, expZero;
	expZero = threshold;
	string line, sample;
	string result;
	ostringstream convert;
	//tid: the first number of prob line denotes the cluster of each row

        inFile.open(tr_info_file.data(), ifstream::in);
//	inFile.ignore(3);       
        inFile.ignore(1000,' ');       
        inFile >> kRest;
        cout<< "kRest = "<< kRest<<endl;
	//alpha_prior[K-1] = kRest + 0.0;
        inFile >> temp;
        inFile >> cluster;
        cout<< "cluster = "<< cluster<<endl;
	myClustID = cluster;
        inFile.close();
	
	srand(time(NULL)) ;
	//clock_t start, end;
	//double index[K];

	inFile.open("conditionA.weights", ifstream::in);
	if (!inFile) {
		cerr << "Unable to open conditionA.weights file\n";
		exit(1);   // call system to stop
	}
	double malakia;
	for(h = 0; h < n_1; h++){
		inFile >> malakia;
		papasA.push_back(malakia);
	}
	inFile.close();


	inFile.open("conditionB.weights", ifstream::in);
	if (!inFile) {
		cerr << "Unable to open conditionB.weights file\n";
		exit(1);   // call system to stop
	}

	for(h = 0; h < n_2; h++){
		inFile >> malakia;
		papasB.push_back(malakia);
	}
	inFile.close();		


	for(h = 0; h < n_1; h++){
		cout << papasA[h]<< " ";
	}
	cout << endl;

	for(h = 0; h < n_2; h++){
		cout << papasB[h]<< " ";
	}
	cout << endl;
	//exit(1);


	//read conditionA files
	cout<<"Reading Condition A"<<endl;
	for (h = 0; h < n_1; h++){
		convert << h + 1;
		sample = prob_file_1.data();
		sample += "_";
		sample += convert.str();
		sample += ".panos";
		//cout<<"name = "<< sample.data()<<endl;
		convert.str("");

		//for (i = 0; i < K; ++i) index[i] = 1.0;
		inFile.open(sample.data(), ifstream::in);
		cout<< "              Reading "<< sample.data()<<" file: ";
		if (!inFile) {
			cerr << "Unable to open file .panos file\n";
			exit(1);   // call system to stop
		}
		i = -1;
		while (i<0){
			getline( inFile, line );
			//cout <<line <<  endl;
			if (line[0] != '#'){
				i = 0;	
				inFile.seekg(place);	// this sets the ifstream to the previous line
			}
			place = inFile.tellg(); // this reads the current position of the ifstream
		}
		i = 0;
		inFile >> tid;
		if (tid == cluster){
			inFile >> num;
			y_1.push_back(num);
			xi_alloc.push_back(num);
			vector<double> row(num);
			vector<int> c_row(num);
			temp = 0.0;
			for ( j = 0; j < num; ++j){
				inFile >> c_row[j] >> row[j];
				row[j] = exp(row[j]);
				if (row[j] < threshold){row[j] = threshold;nan_number++;}	
			}
			x_1.push_back(row);
			c_1.push_back(c_row);	
			inFile.ignore(10000000,'\n');
			i = i+1;		
		}else{
			inFile.ignore(10000000,'\n');
			//cout << "line ignored"<< endl;
		}

		while( inFile.good() ){

			inFile >> tid;
			if (tid == cluster){
				inFile >> num;
				y_1.push_back(num);
				xi_alloc.push_back(num);
				vector<double> row(num);
				vector<int> c_row(num);
				//cout << "i = "<< i + 1 << ", nMap = "<< num << "  ";
				for ( j = 0; j < num; ++j){
					inFile >> c_row[j] >> row[j];
					row[j] = exp(row[j]);
					if (row[j] < threshold){row[j] = threshold;nan_number++;}
					//cout << c_row[j] << ", " << row[j]<< " ";	
					
				}
				//cout << endl;
				x_1.push_back(row);
				x_1[i][num - 1] = expZero;
				c_1.push_back(c_row);
				inFile.ignore(10000000,'\n');
				i = i+1;
			}else{
				inFile.ignore(10000000,'\n');
			}
		}

		inFile.close();
		n_obs = i;
		r.push_back(n_obs);
		cout << "replicate "<< h + 1 << ", r["<<  h<<"] = "<< r[h] <<", total reads = "<< x_1.size()  <<endl;
		step += n_obs;
		cum_sum_1[h] = step;

	}
	if(nan_number>0){
		cout<<"WARNING: found "<< nan_number <<" alignments with weird pdf value."<<endl;
		nan_number = 0;
	}
	cout<<"Reading Condition B"<<endl;
	step = 0;
	for (h = 0; h < n_2; h++){
		convert << h + 1;
		sample = prob_file_2.data();
		sample += "_";
		sample += convert.str();
		sample += ".panos";
		convert.str("");

		//for (i = 0; i < K; ++i) index[i] = 1.0;
		inFile.open(sample.data(), ifstream::in);
		cout<< "              Reading "<< sample.data()<<" file: ";
		if (!inFile) {
			cerr << "Unable to open file .panos file\n";
			exit(1);   // call system to stop
		}

		i = -1;
		while (i<0){
			getline( inFile, line );
			if (line[0] != '#'){
				i = 0;	
				inFile.seekg(place);	// this sets the ifstream to the previous line
			}
			place = inFile.tellg(); // this reads the current position of the ifstream
		}
		i = 0;
		inFile >> tid; 
		if (tid == cluster){
			inFile >> num;
			y_2.push_back(num);
			z_alloc.push_back(num);
			vector<double> row(num);
			vector<int> c_row(num);
			temp = 0.0;
			for ( j = 0; j < num; ++j){
				inFile >> c_row[j] >> row[j];
				row[j] = exp(row[j]);
				if (row[j] < threshold){row[j] = threshold;nan_number++;}		
			}
			x_2.push_back(row);
			c_2.push_back(c_row);	
			inFile.ignore(10000000,'\n');
			i = i+1;
			}else{
				inFile.ignore(10000000,'\n');
			}
		

		while( inFile.good() ){

			inFile >> tid; 
			if (tid == cluster){
				inFile >> num;
				y_2.push_back(num);
				z_alloc.push_back(num);
				vector<double> row(num);
				vector<int> c_row(num);
				for ( j = 0; j < num; ++j){
					inFile >> c_row[j] >> row[j];
					row[j] = exp(row[j]);
					if (row[j] < threshold){row[j] = threshold;nan_number++;}	
				}
				x_2.push_back(row);
				x_2[i][num - 1] = expZero;
				c_2.push_back(c_row);
				inFile.ignore(10000000,'\n');
				i = i+1;
			}else{
				inFile.ignore(10000000,'\n');
			}
		}

		inFile.close();
		n_obs = i ;
		s.push_back(n_obs) ;
		cout << "replicate "<< h +1 << ", s["<<  h<<"] = "<< s[h] <<", total reads = "<< x_2.size()  <<endl;
		step += n_obs;
		cum_sum_2[h] = step;
	}
	if(nan_number>0){
		cout<<"WARNING: found "<< nan_number <<" alignments with weird pdf value."<<endl;
	}
//*******************************************************************************************


	int  initial_state;
	double t_start, t_end;
	vector<double> a_pars, b_pars;
	int chain, switchCondition = 0;
	t_start = omp_get_wtime( );

//	ofstream de_file ( output_file.data() );
//	ofstream theta_out ( theta_file.data() );
//	ofstream w_out ( w_file.data() );
//	ofstream gThan ( xgty_file.data() );

	vector<double> thetaMEAN (K,0.0);
	vector<double> wMEAN (K,0.0);	
	vector<double> stateMEAN (K,0.0);	
	int nIterations = 0;
	string myString;
	int numberChains = nChains;
	if (x_2.size() + x_1.size()> 2000000){
			numberChains = 4;
	}
	if (x_2.size() + x_1.size()> 10000000){
			numberChains = 2;
	}
	if (x_2.size() + x_1.size()> 30000000){
			numberChains = 2;
			iterNEW = (iterations - 1)/5 + 1;
			myBurn = 100;
	}


	//if (K > 200){ numberChains = 2; }

        ofstream ar_out ( ar_file.data() );

	for(chain = 0; chain < numberChains; chain++){ 
		

		convert << chain + 1;
		myString = output_file.data();
		myString += "_";
		myString += convert.str();
		cout<<"name = "<< myString.data()<<endl;
		convert.str("");
		ofstream de_file ( myString.data() );

		if (chain == numberChains/2){
			cout << "swap:"<<endl;
			x_1.swap(x_2);
			y_1.swap(y_2);
			c_1.swap(c_2);
			z_alloc.swap(xi_alloc);
			cum_sum_1.swap(cum_sum_2);
			papasA.swap(papasB);
			r.swap(s);
			switchCondition = 1;

		}


		convert << chain + 1;
		myString = theta_file.data();
		myString += "_";
		myString += convert.str();
		cout<<"name = "<< myString.data()<<endl;
		convert.str("");
		ofstream theta_out ( myString.data() );

		convert << chain + 1;
		myString = w_file.data();
		myString += "_";
		myString += convert.str();
		cout<<"name = "<< myString.data()<<endl;
		convert.str("");
		ofstream w_out ( myString.data() );

		//if( chain < numberChains/2){
		//	initial_state = 0;
		//}else{
		//	initial_state = 1;
		//}
		if ((chain+1) % 2 == 0){
			initial_state = 1;
		}else{
			initial_state = 0;
		}
		


		//if (K > 200){ initial_state = 0; }


	//**************************************************************************************************
	//	Algorithm begins by default from full DE. Otherwise:					**
	//												**
		//this corresponds to initialization from zero.						**
		if (initial_state == 0){						//		**
			for (k = 1; k< K - 1; k++){					//		**
				c_vec[k] = 0;						//		**
			}
			c_vec[0] = 1;							//		**
			c_vec[K-1] = 1;
		}else{
			for (k = 1; k< K - 1; k++){					//		**
				c_vec[k] = 1;						//		**
			}
		}


		//initialization
		iter = 0;
		c_sum =accumulate(c_vec.begin(),c_vec.end(),0);
		//gamma_prior[c_sum - 1] = kRest + 0.0;
		cout << "Initializing chain " << chain+1 << " from nDE = " << c_sum <<"... " << endl;

	/*
	##################################################################################################################################
	*/
		dir(alpha_prior,K);	//generates v
		for (k = 0; k<K; k++){
			u[k] = v[k];	//set u equal to v
		}
		permutation_current();	//this computes the permutation and sets theta as a function of u. It should be called AFTER u.
		//computing w
		if (c_sum>0){
			dir(gamma_prior,c_sum);
			for (k = 0; k < K - c_sum; k++){
				w_temp[k] = u[k];
			}
			for (k = K-c_sum; k < K; k++){
				w_temp[k] = v[k - K + c_sum]*u_sum_c1;
			}
			for (k = 0; k < K; k++){
				w[k] = w_temp[inv_permutation[k]];
			}
		}else{
			for (k = 0; k< K; k++){
				w[k] = theta[k];
			}
		}


	/*
	##################################################################################################################################
	*/

		int sample1, sample2;
		double cum_sum1, cum_sum2, tot_sum1, tot_sum2, uni, rj_rate = 0.0;
	/*			Burn-in period
	##################################################################################################################################
	*/

		tot_sum1 = r[0]/papasA[0] + s[0]/papasB[0] + K;
		tot_sum2 = r[0]/papasA[0] + K;
		sample1 = 1;
		sample2 = 1; 
		ty1 = r[0]/papasA[0] - r[0];
		ty2 = s[0]/papasB[0] - s[0];

		for(iter = 1; iter<burn;iter++){
			//generate random integers in 1,...,n_1 and 1,...,n_2
			//sample1 = 1 + rand() % n_1;
			//sample2 = 1 + rand() % n_2; 
			allocate(sample1, sample2);	// generate allocations 
			// compute GD hyperparameters
			//tot_sum1 = r[sample1 - 1]/papasA[sample1 - 1] + s[sample2 - 1]/papasB[sample2 - 1] + K;
			//tot_sum2 = r[sample1 - 1]/papasA[sample1 - 1] + K;
			cum_sum1 = 0.0;
			cum_sum2 = 0.0;
			for (k = 0; k < K - c_sum; k++){
				alpha_gd[k] = alpha_prior[permutation[k]] + sum_1[permutation[k]] + sum_2[permutation[k]];
				beta_gd[k] = tot_sum1 - (alpha_gd[k] + cum_sum1);
				cum_sum1 += alpha_gd[k];
				cum_sum2 += alpha_gd[k] - sum_2[permutation[k]];
			}
			if (c_sum>0){
				for (k = K - c_sum; k < K ; k++){
					alpha_gd[k] = alpha_prior[permutation[k]] + sum_1[permutation[k]];
					beta_gd[k] = tot_sum2 - (alpha_gd[k] + cum_sum2);
					cum_sum2 += alpha_gd[k];
					w_temp[k - K + c_sum] = gamma_prior[permutation[k -K + c_sum]] + sum_2[permutation[k]];			
				}	
			}

			gen_dir(alpha_gd,beta_gd,K);	//generate u
			permutation_current();
			log_like();
			if (c_sum>0){
				dir(w_temp,c_sum);	//generate v
				for (k = 0; k < K - c_sum; k++){
					w_temp[k] = u[k];
				}
				for (k = K-c_sum; k < K; k++){
					w_temp[k] = v[k - K + c_sum]*u_sum_c1;
				}
				for (k = 0; k < K; k++){
					w[k] = w_temp[inv_permutation[k]];
				}
			}else{
				for (k = 0; k< K; k++){
					w[k] = theta[k];
				}
			}
		}
		t_end = omp_get_wtime( );
		cout << "Burn-in done. Time: "<< setprecision(5)<< t_end - t_start<<" sec."<< "\n\n";

	/*			Sampler
	##################################################################################################################################
	*/

		vector<double> xgty (K,0.0);
		vector<double> ygtx (K,0.0);
		vector<double> cDist (2,0.0);
		tot_sum1 = r[0]/papasA[0] + s[0]/papasB[0] + K;
		tot_sum2 = r[0]/papasA[0] + K;
		vector<int> componentsIndex(K-2);

		
		//cout << tot_sum1 << ", " << tot_sum2 <<endl;
		sample1 = 1;
		sample2 = 1; 
		gPrior = beta_rand(0.5,0.5);
		for(iter = 1; iter<iterNEW;iter++){
			//cout << iter <<endl;
			allocateCollapsed(sample1, sample2);	// generate allocations (and compute c_sum) 
			permutation_current_collapsed();

			if (K < 20){
				for (k = 1; k < K - 1; k++){
					//k = 1 + rand() % (K-2);
					for (i = 0; i < K; i++){
						c_new[i] = c_vec[i];
					}
					uni = unifRand();
					c_new[k] = 0;
					cDist[0] = log_f_c_new() + log(1.0-gPrior);
					c_new[k] = 1;
					cDist[1] = log_f_c_new() + log(gPrior);
					cDist[0] = 1.0/(1.0 + exp(cDist[1] - cDist[0]));
					//cout << cDist[0] << endl;
					if(uni < cDist[0]){
						c_vec[k] = 0;
						}else{
						c_vec[k] = 1;
					}

				}
			}else{
				for (k = 0; k < K - 2; k++) componentsIndex[k] = k + 1; // 1 2 ... K - 2
				random_shuffle ( componentsIndex.begin(), componentsIndex.end() ); //random permutation
				// update first 5 elements of random permutation
				for (j = 0; j < 5; j++){
					k = componentsIndex[j];
					for (i = 0; i < K; i++){
						c_new[i] = c_vec[i];
					}
					uni = unifRand();
					c_new[k] = 0;
					cDist[0] = log_f_c_new() + log(1.0-gPrior);
					c_new[k] = 1;
					cDist[1] = log_f_c_new() + log(gPrior);
					cDist[0] = 1.0/(1.0 + exp(cDist[1] - cDist[0]));
					//cout << cDist[0] << endl;
					if(uni < cDist[0]){
						c_vec[k] = 0;
						}else{
						c_vec[k] = 1;
					}

				}				
			}
			c_sum = 0;
			for (k = 0; k < K; k++){
				c_sum += c_vec[k];
				gamma_prior[k] = 1.0;
			}
			permutation_current_collapsed(); //update permutation
//			Update gPrior:
			gPrior = beta_rand(c_sum + 0.5, K - c_sum + 0.5);
			//gPrior = 0.99;
			
			if (iter % saveStep == 0){
				if(iter > myBurn){
					nIterations++;
				}
//###########################################################################################
//###########################################################################################
				cum_sum1 = 0.0;
				cum_sum2 = 0.0;
				if(c_sum<K){
					for (k = 0; k < K - c_sum; k++){
						alpha_gd[k] = alpha_prior[permutation[k]] + sum_1[permutation[k]] + sum_2[permutation[k]];
						beta_gd[k] = tot_sum1 - (alpha_gd[k] + cum_sum1);
						cum_sum1 += alpha_gd[k];
						cum_sum2 += alpha_gd[k] - sum_2[permutation[k]];
					}
				}
				if (c_sum>0){
					for (k = K - c_sum; k < K ; k++){
						alpha_gd[k] = alpha_prior[permutation[k]] + sum_1[permutation[k]];
						beta_gd[k] = tot_sum2 - (alpha_gd[k] + cum_sum2);
						cum_sum2 += alpha_gd[k];
						w_temp[k - K + c_sum] = gamma_prior[permutation[k -K + c_sum]] + sum_2[permutation[k]];			
					}	
				}

				gen_dir(alpha_gd,beta_gd,K);
				permutation_current();	
				if (c_sum>0){
					dir(w_temp,c_sum);
					for (k = 0; k < K - c_sum; k++){
						w_temp[k] = u[k];
					}
					for (k = K-c_sum; k < K; k++){
						w_temp[k] = v[k - K + c_sum]*u_sum_c1;
					}
					for (k = 0; k < K; k++){
						w[k] = w_temp[inv_permutation[k]];
					}
				}else{
					for (k = 0; k< K; k++){
						w[k] = theta[k];
					}
				}
//###########################################################################################
//###########################################################################################
				for (i = 0; i<K; ++i){
					de_file << c_vec[i] << " ";
					if (switchCondition == 0){
						theta_out << theta[i] << " ";
						w_out << w[i] << " ";
					}else{
						theta_out << w[i] << " ";
						w_out << theta[i] << " ";
					}
					if(iter > myBurn){
						if (switchCondition == 0){
							thetaMEAN[i] += theta[i];
							wMEAN[i] += w[i];
						}else{
							thetaMEAN[i] += w[i];
							wMEAN[i] += theta[i];
						}
						stateMEAN[i] += c_vec[i];
					}
				}
				de_file << endl;
				theta_out << endl;
				w_out << endl;
				if (iter % 1000 == 0){
					t_end = omp_get_wtime( );
					cout << "Iteration " << iter <<": RJ rate = "<<setprecision(3)<<100.0*rj_rate/iter <<"%, time: "<< setprecision(5)<< t_end - t_start<<" sec."<< "\n\n";
					cout << sum_1[0] << ", " << sum_2[0]<<"\n";
				}
			}
	//##################################################################################################################
		}
		ar_out << rj_rate/iter << endl;

	}

	ofstream de_file ( output_file.data() );
	ofstream theta_out ( theta_file.data() );
	ofstream w_out ( w_file.data() );


	for (i = 0; i<K; ++i){
		thetaMEAN[i] /= (nIterations + 0.0);
		wMEAN[i] /= (nIterations + 0.0);
		stateMEAN[i] /= (nIterations + 0.0);
		de_file << stateMEAN[i] << endl;
		theta_out << thetaMEAN[i] << endl;
		w_out << wMEAN[i] << endl;
	}

	cout << "Sampling done."<<endl;

return(0);
}





void gen_dir( vector<double>& a_pars, vector<double>& b_pars, int n_comps){
	int i;
	//double *theta_vec = new double[n_comps]; 
	double sum_theta, g1, g2; // g1 and g2 are the two gamma rv's in order to simulate a beta rv
	int seed = rand();
	boost::mt19937 rng(seed);
	i = 0;
	boost::gamma_distribution<double> gd1(a_pars[i]), gd2(b_pars[i]);
	boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma1( rng, gd1 ), var_gamma2( rng, gd2 );
	g1 = var_gamma1(); // g1 ~ Gamma(a_pars[], 1)
	g2 = var_gamma2(); // g2 ~ Gamma(b_pars[], 1)	
	u[i] = g1/(g1 + g2); // theta[] ~ Beta(a_pars[], b_pars[])
	sum_theta = u[i];
	for (i = 1; i< n_comps-1; ++i){
		boost::gamma_distribution<double> gd1(a_pars[i]), gd2(b_pars[i]);
		boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma1( rng, gd1 ), var_gamma2( rng, gd2 );
		g1 = var_gamma1(); // g1 ~ Gamma(a_pars[], 1)
		g2 = var_gamma2(); // g2 ~ Gamma(b_pars[], 1)	
		u[i] = g1*(1.0-sum_theta)/(g1 + g2); // theta[] ~ Beta(a_pars[], b_pars[])
		sum_theta += u[i];
	}
	u[n_comps-1] = 1.0-sum_theta;

}

/*
*****************************************************************************************
 	Function to simulate Beta(a1,a2) r.v						*
*****************************************************************************************								
*/

double  beta_rand( double a1, double a2){
	double theta_vec, g1, g2; // g1 and g2 are the two gamma rv's in order to simulate a beta rv
	int seed = rand();
	boost::mt19937 rng(seed);
	boost::gamma_distribution<double> gd1(a1), gd2(a2);
	boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma1( rng, gd1 ), var_gamma2( rng, gd2 );
	g1 = var_gamma1(); // g1 ~ Gamma(a_pars[], 1)
	g2 = var_gamma2(); // g2 ~ Gamma(b_pars[], 1)	
	theta_vec = g1/(g1 + g2); // theta[] ~ Beta(a_pars[], b_pars[])
	return theta_vec;
}

/*
*****************************************************************************************
 	Function to compute the log pdf Beta(a1,a2) 					*
*****************************************************************************************								
*/


double log_beta(double b, double a1, double a2){
	double pdf;
	pdf = boost::math::lgamma(a1) + boost::math::lgamma(a2) - boost::math::lgamma(a1+a2); //beta function
	pdf = (a1-1.0)*log(b) + (a2-1.0)*log(1.0-b) - pdf;
	return pdf;
}


/*
*****************************************************************************************
 	Function to compute the log f(xi_alloc,z_alloc|c), given xi_alloc, z_alloc, c	*
*****************************************************************************************								
*/


double log_f_c_old(){
	double pdf, u1, talHeads, talHeads2, sum_lg_c0, lg_sum_c1, lg_sum_c1_2, sum_lg_c1, lg_sum_c1_3, sum_lg_c1_3;
	int k;
	sum_lg_c0 = 0.0;
	lg_sum_c1 = 0.0;
	lg_sum_c1_2 = 0.0;
	sum_lg_c1 = 0.0;
	lg_sum_c1_3 = 0.0; 
	sum_lg_c1_3 = 0.0;
	c_sum = 0;
	for(k = 0; k < K; k++){
		talHeads2 = alpha_prior[k] + sum_1[k]; 
		talHeads = talHeads2 + sum_2[k];
		if(c_vec[k] == 0){
			sum_lg_c0 +=  boost::math::lgamma(talHeads);
		}else{
			lg_sum_c1 += talHeads;
			lg_sum_c1_2 += talHeads2;
			sum_lg_c1 +=  boost::math::lgamma(talHeads2);
			u1 = gamma_prior[0] + sum_2[k];
			sum_lg_c1_3 += boost::math::lgamma(u1);
			lg_sum_c1_3 += u1;
			c_sum++;
		}	
	}
	lg_sum_c1 = boost::math::lgamma(lg_sum_c1);
	lg_sum_c1_2 = boost::math::lgamma(lg_sum_c1_2);
	lg_sum_c1_3 = boost::math::lgamma(lg_sum_c1_3);
	pdf = lg_sum_c1 + sum_lg_c0 + sum_lg_c1 + sum_lg_c1_3 - lg_sum_c1_2 - lg_sum_c1_3; 
	return pdf;
}

/*
*****************************************************************************************
 	Function to compute the log f(xi,z|c_new) 					*
*****************************************************************************************								
*/

double log_f_c_new(){
	double pdf, u1, talHeads, talHeads2, sum_lg_c0, lg_sum_c1, lg_sum_c1_2, sum_lg_c1, lg_sum_c1_3, sum_lg_c1_3;
	int k;

	sum_lg_c0 = 0.0;
	lg_sum_c1 = 0.0;
	lg_sum_c1_2 = 0.0;
	sum_lg_c1 = 0.0;
	lg_sum_c1_3 = 0.0; 
	sum_lg_c1_3 = 0.0;
	for(k = 0; k < K; k++){
		talHeads2 = alpha_prior[k] + sum_1[k]; 
		talHeads = talHeads2 + sum_2[k];
		if(c_new[k] == 0){
			sum_lg_c0 +=  boost::math::lgamma(talHeads);
		}else{
			lg_sum_c1 += talHeads;
			lg_sum_c1_2 += talHeads2;
			sum_lg_c1 +=  boost::math::lgamma(talHeads2);
			u1 = gamma_prior[0] + sum_2[k];
			sum_lg_c1_3 += boost::math::lgamma(u1);
			lg_sum_c1_3 += u1;
		}	
	}
	lg_sum_c1 = boost::math::lgamma(lg_sum_c1);
	lg_sum_c1_2 = boost::math::lgamma(lg_sum_c1_2);
	lg_sum_c1_3 = boost::math::lgamma(lg_sum_c1_3);
	pdf = lg_sum_c1 + sum_lg_c0 + sum_lg_c1 + sum_lg_c1_3 - lg_sum_c1_2 - lg_sum_c1_3; 
	return pdf;
}



void dir( vector<double>& a_pars, int n_comps){
	int i;

	double sum_theta; 
	int seed = rand();
	boost::mt19937 rng(seed);

	i = 0;
	//boost::gamma_distribution<double> gd1(a_pars[i]);
	//boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma1( rng, gd1 );
	sum_theta = 0.0;
	for (i = 0; i< n_comps; ++i){
		boost::gamma_distribution<double> gd1(a_pars[i]);
		boost::variate_generator<boost::mt19937&,boost::gamma_distribution<> > var_gamma1( rng, gd1 );
		v[i] = var_gamma1(); // g1 ~ Gamma(a_pars[], 1)
		sum_theta += v[i];
	}
	for (i = 0; i< n_comps; ++i){
		v[i] = v[i]/sum_theta; // theta[] ~ Beta(a_pars[], b_pars[])
	}

}



/*
*****************************************************************************************
 	Function to simulate U(0,1) r.v's						*
*****************************************************************************************								
*/

//double unifRand()
//{
//	int seed = rand();
//	boost::mt19937 rng(seed);
//	boost::uniform_01<double> gd1;
//	boost::variate_generator<boost::mt19937&,boost::uniform_01<> > uni( rng, gd1 );
  //  return uni();
//}


double unifRand()
{
  return rand() / double(RAND_MAX);
}



double allocateCollapsed(int sample1, int sample2){
int z,i,k,sample_start, sample_end;
double uni, z_prob_sum, z_dist, u1, u2, talHeads;

sample_start = cum_sum_1[sample1-1] - r[sample1-1];
sample_end = cum_sum_1[sample1-1];



// sum_1 = {s_k(xi); k = 1, ..., K}
// sum_2 = {s_k(z); k = 1, ..., K}
	for (i = sample_start; i < sample_end; i++){
		z = xi_alloc[i];
		//cout<<", "<< i << ", previous alloc_i = " << z << "  ---------"<< endl;
		u1 = 0.0;
		u2 = 0.0;
		z_prob_sum = 0.0;
		for (k = 0; k < K; k++){
			talHeads = alpha_prior[k] + sum_1[k];
			if(c_vec[k] == 0){
				theta[k] = talHeads + sum_2[k];
			}else{
				theta[k] = talHeads;
				u1 += talHeads + sum_2[k];
				u2 += talHeads;
			}
		}
		theta[z] = theta[z] - 1.0;
		if(c_vec[z] == 1){
			u1 = u1 - 1.0;
			u2 = u2 - 1.0;
		}
		for (k = 0; k < y_1[i]; k++){
			if(c_vec[c_1[i][k]] == 0){
				theta[c_1[i][k]] *= x_1[i][k];
			}else{
				theta[c_1[i][k]] *= u1*x_1[i][k]/u2;
			}
			z_prob_sum += theta[c_1[i][k]];			
		}

		k = 0;
		z_prob[k] = theta[c_1[i][k]]/z_prob_sum;
		z_dist = z_prob[k];
		uni = unifRand();
		xi_alloc[i] = c_1[i][k];
		if (z_prob[k] != z_prob[k]){
			cout <<"oops"<<endl;z_prob[k] = threshold;
			//cout << " i = "<< i << ", f_i_k = " << x_1[i][k] << ", sum = "<< z_prob_sum << ", theta = " << theta[c_1[i][k]] <<endl;	
			//cin.get();
		}
		while(uni > z_dist && k < y_1[i] - 1 ){
			k++;
			z_prob[k] = theta[c_1[i][k]]/z_prob_sum;
			if (z_prob[k] != z_prob[k]){cout <<"oops"<<endl;z_prob[k] = threshold;
					cout <<"oops"<<endl;z_prob[k] = threshold;
					//cout << " i = "<< i << ", f_i_k = " << x_1[i][k]<< ", sum = "<< z_prob_sum <<endl;	
					//cin.get();
			}
			z_dist = z_dist + z_prob[k];
//			cout <<"k = "<<k<<", u = "<< uni<<", z_dist = " <<z_dist<<", comp = "<<c_1[i][k]<<endl;
//			cout << " i = "<< i << ", f_i_k = " << x_1[i][k]<< ", sum = "<< z_dist <<endl;	
		}
//		cout << "andronikos: " << k <<endl;
		if ( k >= y_1[i]){cout<<"shit"<<endl;k = y_1[i] - 1;}
//		cout << "talking heads: " << c_1[i][k] <<endl;
		//if(c_1[i][k] > K ){xi_alloc[i] = 0;}else{
		xi_alloc[i] = c_1[i][k];
//		cout << "xi_alloc: " << xi_alloc[i] << endl;
		//}
		sum_1[xi_alloc[i]]++; 
//		cout << "edw: " << sum_1[xi_alloc[i]] << endl;
		sum_1[z] += -1; 
//		cout << "ekei: " << sum_1[z] << endl;
	}

//cout << "edw"<<endl;
// condition B
sample_start = cum_sum_2[sample2-1] - s[sample2-1];
sample_end = cum_sum_2[sample2-1] ;

	for (i = sample_start; i < sample_end; i++){
		//cout<< "B: "<< i << endl;
		z = z_alloc[i];
		u1 = 0.0;
		u2 = 0.0;
		z_prob_sum = 0.0;
		for (k = 0; k < K; k++){
			talHeads = alpha_prior[k] + sum_1[k];
			if(c_vec[k] == 0){
				w[k] = talHeads + sum_2[k];
			}else{
				w[k] = gamma_prior[0] + sum_2[k];
				u1 += talHeads + sum_2[k];
				u2 += w[k];
			}
		}
		w[z] = w[z] - 1.0;
		if(c_vec[z] == 1){
			u1 = u1 - 1.0;
			u2 = u2 - 1.0;
		}
		for (k = 0; k < y_2[i]; k++){
			if(c_vec[c_2[i][k]] == 0){
				w[c_2[i][k]] *= x_2[i][k];
			}else{
				w[c_2[i][k]] *= u1*x_2[i][k]/u2;
			}
			z_prob_sum += w[c_2[i][k]];			
		}

		k = 0;
		z_prob[k] = w[c_2[i][k]]/z_prob_sum;
		z_dist = z_prob[k];
		uni = unifRand();
		z_alloc[i] = c_2[i][k];
		//cout <<"k = "<<k<<", u = "<< uni<<", z_dist = " <<z_dist <<", comp = "<<c_1[i][k]<<endl;
		if (z_prob[k] != z_prob[k]){
			cout <<"oops"<<endl;z_prob[k] = threshold;
			//cout << " i = "<< i << ", f_i_k = " << x_1[i][k] << ", sum = "<< z_prob_sum << ", theta = " << theta[c_1[i][k]] <<endl;	
			//cin.get();
		}
		while(uni > z_dist && k < y_2[i]-1 ){
			k++;
			z_prob[k] = w[c_2[i][k]]/z_prob_sum;
			if (z_prob[k] != z_prob[k]){cout <<"oops"<<endl;z_prob[k] = threshold;
					cout <<"oops"<<endl;z_prob[k] = threshold;
					//cout << " i = "<< i << ", f_i_k = " << x_1[i][k]<< ", sum = "<< z_prob_sum <<endl;	
					//cin.get();
			}
			z_dist = z_dist + z_prob[k];
			//cout <<"k = "<<k<<", u = "<< uni<<", z_dist = " <<z_dist<<", comp = "<<c_1[i][k]<<endl;
		}

		if ( k >= y_2[i]){cout<<"shit"<<endl;k = y_2[i] - 1;}
		z_alloc[i] = c_2[i][k];
		sum_2[z_alloc[i]]++; 
		sum_2[z] += -1; 
	}

//cout << "cluster: "<< myClustID<< ": " << sum_1[K-1] <<", "<< sum_2[K-1]<<endl;
//exit(1);



return(1.0);
}


double allocate(int sample1, int sample2){

int z,i,k,sample_start, sample_end;
double uni, z_prob_sum, z_dist;

//cout<< sample1 <<" "<< sample2<<endl;

//compute start end end of 1st sample
sample_start = cum_sum_1[sample1-1] - r[sample1-1];
sample_end = cum_sum_1[sample1-1];
//cout<<endl;
//cout<<sample_start<<"  "<<sample_end<<endl;
//cin.get();
// set all partial sums to zero
std::fill (sum_1.begin(),sum_1.end(),0.0);
std::fill (sum_2.begin(),sum_2.end(),0.0);
// allocate sample1

p_alloc = 0.0;
log_likelihood = 0.0;

for (i = sample_start; i < sample_end; i++){
	z_prob_sum = 0.0;
	for (k = 0; k < y_1[i]; k++){
		z_prob[k] = theta[c_1[i][k]]*x_1[i][k];
		z_prob_sum += z_prob[k];
	}
	k = 0;
	z_prob[k] /= z_prob_sum;
	z_dist = z_prob[k];
	uni = unifRand();
	z = c_1[i][k];
	//cout <<"k = "<<k<<", u = "<< uni<<", z_dist = " <<z_dist <<", comp = "<<c_1[i][k]<<endl;
	if (z_prob[k] != z_prob[k]){
		cout <<"oops"<<endl;z_prob[k] = threshold;
		//cout << " i = "<< i << ", f_i_k = " << x_1[i][k] << ", sum = "<< z_prob_sum << ", theta = " << theta[c_1[i][k]] <<endl;	
		//cin.get();
	}
	while(uni > z_dist && k < y_1[i] ){
		k++;
		z_prob[k] /= z_prob_sum;
		if (z_prob[k] != z_prob[k]){cout <<"oops"<<endl;z_prob[k] = threshold;
				cout <<"oops"<<endl;z_prob[k] = threshold;
				//cout << " i = "<< i << ", f_i_k = " << x_1[i][k]<< ", sum = "<< z_prob_sum <<endl;	
				//cin.get();
		}
		z_dist = z_dist + z_prob[k];
		//cout <<"k = "<<k<<", u = "<< uni<<", z_dist = " <<z_dist<<", comp = "<<c_1[i][k]<<endl;
	}

	if ( k >= y_1[i]){cout<<"shit"<<endl;k = y_1[i] - 1;}
	//if ( k >= y_1[i]){cout<<"shit"<<endl;k = 1;}
	z = c_1[i][k];

	xi_alloc[i] = z;
	//cout <<"z = "<<z<<endl;
	//cout <<"-------------------------------- "<<endl;
	//cin.get();
	sum_1[z]++; 
	p_alloc += log(z_prob[k]);
	log_likelihood += log(x_1[i][k]);
}


// allocate sample2

sample_start = cum_sum_2[sample2-1] - s[sample2-1];
sample_end = cum_sum_2[sample2-1] ;
//cout<<endl;
//cout<<sample_start<<"  "<<sample_end<<endl;
//cin.get();
for (i = sample_start; i < sample_end; i++){
	z_prob_sum = 0.0;
	for (k = 0; k < y_2[i]; k++){
		z_prob[k] = w[c_2[i][k]]*x_2[i][k];
		z_prob_sum += z_prob[k];
	}
	k = 0;
	z_prob[k] /= z_prob_sum;
	z_dist = z_prob[k];
	uni = unifRand();
	z = c_2[i][k];
	//cout <<"z_dist = "<< z_dist<<endl;
	if (z_prob[k] != z_prob[k]){cout <<"oops"<<endl; z_prob[k] = threshold;}
	while(uni > z_dist && k < y_2[i]){
		k++;
		z_prob[k] /= z_prob_sum;
		if (z_prob[k] != z_prob[k]){cout <<"oops"<<endl;z_prob[k] = threshold;}
		z_dist = z_dist + z_prob[k];
		//cout <<"z_dist = "<< z_dist<<endl;
		//cout <<"k = "<<k<<", u = "<< uni<<", z_dist = " <<z_dist<<endl;
	}
	if ( k >= y_2[i]){cout<<"shit"<<endl;k = y_2[i] - 1;}
	z = c_2[i][k];
	z_alloc[i] = z;
	sum_2[z]++; 
	p_alloc += log(z_prob[k]);
	//log_likelihood += log(w[z]) + log(x_2[i][k]);
	log_likelihood += log(x_2[i][k]);
	//cout<<z<<","<<log(theta[z])<<endl;
}

//cout<< "P_alloc = "<<	p_alloc<<endl;
//cout<< "logLikelihood 1= "<<	log_likelihood<<endl;

//sum_1[K-1] = r[sample1 - 1]/papasA[sample1 - 1] - r[sample1 - 1];
//sum_2[K-1] = s[sample2 - 1]/papasB[sample2 - 1] - s[sample2 - 1];

sum_1[K-1] = ty1;
sum_2[K-1] = ty2;


//cout << sum_1[K-1] <<", "<< sum_2[K-1]<<endl;
//exit(1);

return(1.0);

}


//this function computes the loglikelihood given that function allocate has been previously called.

void log_like(){

int k;

for (k = 0; k < K; k++){
	log_likelihood += sum_1[k]*log(theta[k]) + sum_2[k]*log(w[k]);
}
//cout<< "logLikelihood = "<<	log_likelihood<<endl;

}


void permutation_current(){

int k, iter1,iter2;

iter1 = 0;
iter2 = 0;
u_sum_c1 = 0.0;
for (k = 0; k<K; k++){
	if (c_vec[k] == 0){
		iter1++;
		permutation[iter1-1] = k;
		inv_permutation[k] = iter1 - 1;
		theta[k] = u[inv_permutation[k]];
	}else{
		iter2++;
		permutation[K - c_sum - 1 + iter2] = k;
		inv_permutation[k] = K - c_sum - 1 + iter2;	
		theta[k] = u[inv_permutation[k]];
		u_sum_c1 += theta[k];
	}
}
}



void permutation_current_collapsed(){
int k, iter1,iter2;
iter1 = 0;
iter2 = 0;
for (k = 0; k<K; k++){
	if (c_vec[k] == 0){
		iter1++;
		permutation[iter1-1] = k;
		inv_permutation[k] = iter1 - 1;
	}else{
		iter2++;
		permutation[K - c_sum - 1 + iter2] = k;
		inv_permutation[k] = K - c_sum - 1 + iter2;	
	}
}
}





/*
####################################################################################
*/

void permutation_prop(){
int k, iter1,iter2;
iter1 = 0;
iter2 = 0;
u_sum_c1_new = 0.0;

for (k = 0; k<K; k++){
	if (c_new[k] == 0){
		iter1++;
		permutation_new[iter1-1] = k;
		inv_permutation_new[k] = iter1 - 1;
	}else{
		iter2++;
		permutation_new[K - c_sum_new - 1 + iter2] = k;
		inv_permutation_new[k] = K - c_sum_new - 1 + iter2;	
		u_sum_c1_new += theta[k];
	}
	u_temp1[k] = u[inv_permutation[k]];
}
for (k = 0; k<K; k++){
	u_temp[k] = u_temp1[permutation_new[k]];
}
}


void permutation_prop_collapsed(){
int k, iter1,iter2;
iter1 = 0;
iter2 = 0;
for (k = 0; k<K; k++){
	if (c_new[k] == 0){
		iter1++;
		permutation_new[iter1-1] = k;
		inv_permutation_new[k] = iter1 - 1;
	}else{
		iter2++;
		permutation_new[K - c_sum_new - 1 + iter2] = k;
		inv_permutation_new[k] = K - c_sum_new - 1 + iter2;	
	}
}
}





/*
####################################################################################
*/


// k1 denotes c_sum
double bp(int k1){
	if (k1 == 2){
		return(0.0);
	}else{
		//return(log(K-k1+0.0) - log(K+0.0));		
		return(log(K-k1+0.0) - log(K-2.0));		
	}
}

double dp(int k1){
	//return(log(k1+0.0) - log(K+0.0));		
	return(log(k1-2.0) - log(K-2.0));		
}


/*
################################################################################
*/


double allocate_new(int sample1, int sample2){

int z,i,k,sample_start, sample_end;
double uni, z_prob_sum, z_dist;

//compute start end end of 1st sample
sample_start = cum_sum_1[sample1-1] - r[sample1-1];
sample_end = cum_sum_1[sample1-1];
//cout<<endl;
//cout<<sample_start<<"  "<<sample_end<<endl;

// set all partial sums to zero
std::fill (sum_1_new.begin(),sum_1_new.end(),0.0);
std::fill (sum_2_new.begin(),sum_2_new.end(),0.0);
// allocate sample1

p_alloc_new = 0.0;
log_likelihood_new = 0.0;

for (i = sample_start; i < sample_end; i++){
        z_prob_sum = 0.0;
        for (k = 0; k < y_1[i]; k++){
                z_prob[k] = theta[c_1[i][k]]*x_1[i][k];
                z_prob_sum += z_prob[k];
        }
        k = 0;
        z_prob[k] /= z_prob_sum;
        z_dist = z_prob[k];
        uni = unifRand();
        z = c_1[i][k];
        //cout <<"z_dist = "<< z_dist<<endl;
        if (z_prob[k] != z_prob[k]){cout <<"oops"<<endl;z_prob[k] = threshold;}
        while(uni > z_dist && k < y_1[i]){
                k++;
                z_prob[k] /= z_prob_sum;
                if (z_prob[k] != z_prob[k]){cout <<"oops"<<endl;z_prob[k] = threshold;}
                z_dist = z_dist + z_prob[k];
                //cout <<"z_dist = "<< z_dist<<endl;
                //cout <<"k = "<<k<<", u = "<< uni<<", z_dist = " <<z_dist<<endl;
        }
        if ( k >= y_1[i]){cout<<"shit"<<endl;k = y_1[i] - 1;}
        z = c_1[i][k];
        sum_1_new[z]++; 
        p_alloc_new += log(z_prob[k]);
        //log_likelihood += log(theta[z]) + log(x_1[i][k]);
        log_likelihood_new += log(x_1[i][k]);
        //cout<<z<<","<<log(theta[z])<<endl;
}
// allocate sample2

sample_start = cum_sum_2[sample2-1] - s[sample2-1];
sample_end = cum_sum_2[sample2-1];
//cout<<endl;
//cout<<sample_start<<"  "<<sample_end<<endl;

for (i = sample_start; i < sample_end; i++){
        z_prob_sum = 0.0;
        for (k = 0; k < y_2[i]; k++){
                z_prob[k] = w_temp[c_2[i][k]]*x_2[i][k];
                z_prob_sum += z_prob[k];
        }
        k = 0;
        z_prob[k] /= z_prob_sum;
        z_dist = z_prob[k];
        uni = unifRand();
        z = c_2[i][k];
        //cout <<"z_dist = "<< z_dist<<endl;
        if (z_prob[k] != z_prob[k]){cout <<"oops"<<endl;z_prob[k] = threshold;}
        while(uni > z_dist && k < y_2[i]){
                k++;
                z_prob[k] /= z_prob_sum;
                if (z_prob[k] != z_prob[k]){cout <<"oops"<<endl;z_prob[k] = threshold;}
                z_dist = z_dist + z_prob[k];
                //cout <<"z_dist = "<< z_dist<<endl;
                //cout <<"k = "<<k<<", u = "<< uni<<", z_dist = " <<z_dist<<endl;
        }
        if ( k >= y_2[i]){cout<<"shit"<<endl;k = y_2[i] - 1;}
        z = c_2[i][k];
        sum_2_new[z]++; 
        p_alloc_new += log(z_prob[k]);
        //log_likelihood += log(w[z]) + log(x_2[i][k]);
        log_likelihood_new += log(x_2[i][k]);
        //cout<<z<<","<<log(theta[z])<<endl;
}

//cout<< "P_alloc = "<< p_alloc<<endl;
//cout<< "logLikelihood 1= "<<  log_likelihood<<endl;

//sum_1_new[K-1] = r[sample1 - 1]/papasA[sample1 - 1] - r[sample1 - 1];
//sum_2_new[K-1] = s[sample2 - 1]/papasB[sample2 - 1] - s[sample2 - 1];


sum_1_new[K-1] = ty1;
sum_2_new[K-1] = ty2;


return(1.0);

}


//this function computes the loglikelihood given that function allocate has been previously called.

void log_like_new(){
int k;

for (k = 0; k < K; k++){
        log_likelihood_new += sum_1_new[k]*log(theta[k]) + sum_2_new[k]*log(w_temp[k]);
}
//cout<< "logLikelihood = "<<   log_likelihood<<endl;

}




/*		Birth acceptance probability
#################################################################################################
*/

double birth_accept_prob(double ll_old, double ll_new,double p_all_old, double p_all_new,int c_s_old,int c_s_new,double rv){

	double lap, jacobian, prop_ratio, prior_ratio, prior_v_ratio, ll_ratio, birth_prop, death_prop, mixBeta;
	int i;

//	prior_ratio = 0.0;					//this is the 1/(2^{K} - K) prior for c
	

	if (c_s_old == 0){
		jacobian = 0.0;
		birth_prop = bp(c_s_old) + log(2.0) - log(K+0.0) - log(K - 1.0);
		death_prop = dp(c_s_new); 
		prior_v_ratio = log(2.0);
		prior_ratio = 2.0*log(gPrior) - 2.0*log(1.0 - gPrior);
		if (prior_c > 0){
			prior_ratio = - log(K - 1.0) - log(K + 0.0); 		//this is used if prior_c = 1
		}
	}else{
		jacobian = (c_s_old - 1.0)*log(1.0 - rv);
		//birth_prop = bp(c_s_old) - log(K - c_s_old + 0.0);
		birth_prop = bp(c_s_old) - log(K - c_s_old);
		//death_prop = dp(c_s_new) - log(c_s_old + 1.0);
		death_prop = dp(c_s_new) - log(c_s_old - 1.0);
		prior_v_ratio = log(c_s_old);
		//prior_v_ratio = (kRest - 1.0)*log(1.0-rv);
		prior_ratio = log(gPrior) - log(1.0 - gPrior);
		if (prior_c > 0){
			prior_ratio = log(c_s_old + 1) - log(K - c_s_old + 0.0); //this is used if prior_c = 1
			//prior_ratio = log(c_s_old - 2.0) - log(K - c_s_old + 0.0); //this is used if prior_c = 1
		}
	}
	mixBeta = 0.0;
	for(i = 0; i<nMIX; i++){
		mixBeta += exp(log_beta(rv, beta_param1, beta_paramNEW[i]))/(nMIX + 0.0);
	}
	mixBeta = log(mixBeta);	
	prop_ratio = p_all_old - p_all_new + jacobian + death_prop - (birth_prop + mixBeta);  
	ll_ratio = ll_new - ll_old;

	lap = ll_ratio + prop_ratio + prior_ratio + prior_v_ratio;
//	cout << lap <<endl;
	return(lap);
}





/*		Birth acceptance probability for the Metropolis-Hastings step
#################################################################################################
*/

double birth_mh(double ll_old, double ll_new, int c_s_old,int c_s_new){

	double lap, prop_ratio, prior_ratio, ll_ratio, birth_prop, death_prop;

//	prior_ratio = 0.0;					//this is the 1/(2^{K} - K) prior for c
	
	if (c_s_old == 0){
		birth_prop = bp(c_s_old) + log(2.0) - log(K+0.0) - log(K - 1.0);
		death_prop = dp(c_s_new); 
		prior_ratio = 2.0*log(gPrior) - 2.0*log(1.0 - gPrior);
		if (prior_c > 0){
			prior_ratio = - log(K - 1.0) - log(K + 0.0); 		//this is used if prior_c = 1
		}
	}else{
		birth_prop = bp(c_s_old) - log(K - c_s_old);
		death_prop = dp(c_s_new) - log(c_s_old - 1.0);
		prior_ratio = log(gPrior) - log(1.0 - gPrior);
		if (prior_c > 0){
			prior_ratio = log(c_s_old + 1) - log(K - c_s_old + 0.0); //this is used if prior_c = 1
		}
	}
	prop_ratio = death_prop - birth_prop;  
	ll_ratio = ll_new - ll_old;
	lap = ll_ratio + prop_ratio + prior_ratio;
	return(lap);
}

















