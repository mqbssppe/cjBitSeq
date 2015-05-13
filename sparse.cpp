#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib> 
#include <string> 
#include <sstream>
#include <iomanip>
#include <ctime>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace boost::numeric::ublas;

extern "C" {
	void sparse (int *Kmax) {


		int i, j, place = 0, num, k1, k2;
		double row;
		int K = Kmax[0];
		//int K = 91;
		//K = 91;
		mapped_matrix<double> m (K, K, K * K);
		ifstream inFile;
		string line, tid;
		int zero;

		inFile.open("conditionA_1.prob", ifstream::in);
		cout<< "              Reading conditionA_1.prob file: " << endl;
		if (!inFile) {
			cerr << "Unable to open file .alignment file\n";
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
		inFile >> tid >> num;
		num = num - 1;
		std::vector<int> c_row(num);
		for ( j = 0; j < num; ++j){
			inFile >> c_row[j] >> row;
		}

		for (k1 = 0; k1 < num; k1++){
			for (k2 = 0; k2 < num; k2++){
				m(c_row[k1],c_row[k2]) = 1;
			}
		}
		inFile >> zero >> row;
		inFile.ignore(10000000,'\n');
		i = i+1;		
		//inFile.good()
		while( !inFile.eof() ){
			inFile >> tid >> num;
			num = num - 1;
			std::vector<int> c_row(num);
			//cout << "row " << i << ", Nmaps = " << num <<". ";
			for ( j = 0; j < num; ++j){
				inFile >> c_row[j] >> row;
			//	cout << c_row[j] << " ";
			}
			for (k1 = 0; k1 < num; k1++){
				for (k2 = 0; k2 < num; k2++){
					m(c_row[k1],c_row[k2]) = 1;
				}
			}
			inFile >> zero >> row;
			inFile.ignore(10000000,'\n');
			i = i+1;
			if(i%1000000 == 0){cout << "read " << i << "lines"<< endl;}
		}
		inFile.close();
		cout << "done: read " << i << "lines"<< endl;
		inFile.open("conditionB_1.prob", ifstream::in);
		cout<< "              Reading conditionB_1.prob file: " << endl;
		if (!inFile) {
			cerr << "Unable to open file .alignment file\n";
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
		inFile >> tid >> num;
		num = num - 1;
		std::vector<int> c_row2(num);
		for ( j = 0; j < num; ++j){
			inFile >> c_row2[j] >> row;
		}

		for (k1 = 0; k1 < num; k1++){
			for (k2 = 0; k2 < num; k2++){
				m(c_row2[k1],c_row2[k2]) = 1;
			}
		}
		inFile >> zero >> row;
		inFile.ignore(10000000,'\n');
		i = i+1;		
		//inFile.good()
		while( !inFile.eof() ){
			inFile >> tid >> num;
			num = num - 1;
			std::vector<int> c_row2(num);
			//cout << "row " << i << ", Nmaps = " << num <<". ";
			for ( j = 0; j < num; ++j){
				inFile >> c_row2[j] >> row;
			//	cout << c_row[j] << " ";
			}
			for (k1 = 0; k1 < num; k1++){
				for (k2 = 0; k2 < num; k2++){
					m(c_row2[k1],c_row2[k2]) = 1;
				}
			}
			inFile >> zero >> row;
			inFile.ignore(10000000,'\n');
			i = i+1;
			if(i%1000000 == 0){cout << "read " << i << "lines"<< endl;}
		}
		inFile.close();
		cout << "done: read " << i << "lines"<< endl;





		ofstream outFile ( "triplet_sparse_matrix.txt" );
		cout << "writing sparse matrix rows and columns: " << endl;
		for (i = 1; i < K; i++){
			for (j = 1; j < K; j++){
				if(m(i,j)>0){
					outFile << i << " " << j << endl;
				}
			}
			if(i%1000 == 0){cout << "Done " << i << "transcripts"<< endl;}

		}
		cout << "done." << endl;


	}

}






