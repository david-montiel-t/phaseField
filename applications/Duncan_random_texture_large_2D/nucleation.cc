// #include <cmath>
// #include"../../include/variableContainer.h"
using namespace std;
// =================================================================================
// NUCLEATION FUNCTIONS
// =================================================================================
// =================================================================================
// Nucleation probability
// =================================================================================
template <int dim, int degree>
double customPDE<dim,degree>::getNucleationProbability(variableValueContainer variable_value, double dV, dealii::Point<dim> p, unsigned int variable_index) const
{
  //Number of order parameters
  unsigned int number_ops_ini = 91;
  
  //Number of empty ops (for nucleation)
  unsigned int number_ops_nuc = 50;
  
  //Number of ops (total)
  unsigned int number_ops_total = number_ops_ini + number_ops_nuc;
	// //Supersaturation factor
    // double ssf;
    // if (dim ==2) ssf=variable_value(0)-calmin;
    // if (dim ==3) ssf=(variable_value(0)-calmin)*(variable_value(0)-calmin);
	
    // Calculate the nucleation rate
	// double J=k1*exp(-k2/(std::max(ssf,1.0e-6)))*exp(-tau/(this->currentTime));
	// double retProb=1.0-exp(-J*userInputs.dtValue*((double)userInputs.steps_between_nucleation_attempts)*dV);

    double nuclProb = 0.0;

    double mult2Op = variable_value(number_ops_total+1);
    double mult3Op = variable_value(number_ops_total+1);
  
    nuclProb = 1.0-exp(-A*(mult2Op+B*mult3Op)*userInputs.dtValue*((double)userInputs.steps_between_nucleation_attempts)*dV);

    return nuclProb;
}
