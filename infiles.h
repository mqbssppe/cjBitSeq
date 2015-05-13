#ifndef infiles
#define infiles
#include <string>
const char vb_file[] = "data.tr";
//const int initial_state = 1;	// 0: initialization from zero DE. 99: initialization from previousState.txt. In any other case initialization from full DE.
const int prior_c = 0;		// if it is set to zero we have the non-informative prior.
const int burn = 500;
const double beta_param1 = 1.0;
const int iterations = 5001;
const int saveStep = 5;
const int nChains = 6;

const std::string output_file = "state_vector.txt";
const std::string theta_file = "theta.txt";
const std::string w_file = "w.txt";
const std::string ar_file = "acceptance.txt";
const std::string prob_file_1 = "../../conditionA";
const std::string prob_file_2 = "../../conditionB";
const std::string tr_info_file = "data.tr";
 
const int n_1 = 1;
const int n_2 = 1;


#endif

