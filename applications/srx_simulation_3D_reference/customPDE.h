#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:
    // Constructor
    customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) , userInputs(_userInputs) {};

    // Function to set the initial conditions (in ICs_and_BCs.h)
    void setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC);

    // Function to set the non-uniform Dirichlet boundary conditions (in ICs_and_BCs.h)
    void setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC);

private:
	#include "../../include/typeDefs.h"

	const userInputParameters<dim> userInputs;

	// Function to set the RHS of the governing equations for explicit time dependent equations (in equations.h)
    void explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

    // Function to set the RHS of the governing equations for all other equations (in equations.h)
    void nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Function to set the LHS of the governing equations (in equations.h)
	void equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Function to set postprocessing expressions (in postprocess.h)
	#ifdef POSTPROCESS_FILE_EXISTS
	void postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
					const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;
	#endif

	// Function to set the nucleation probability (in nucleation.h)
	#ifdef NUCLEATION_FILE_EXISTS
  double getNucleationProbability(variableValueContainer variable_value, double dV, dealii::Point<dim> p, unsigned int variable_index) const;
	#endif

	// ================================================================
	// Methods specific to this subclass
	// ================================================================
  // Method to place the nucleus and calculate the mobility modifier in residualRHS
  //void seedNucleus(const dealii::Point<dim, dealii::VectorizedArray<double> > & q_point_loc,
  //                 dealii::VectorizedArray<double> & source_term,
  //                 dealii::VectorizedArray<double> & gamma) const;
  void seedNucleus(const dealii::Point<dim, dealii::VectorizedArray<double> > & q_point_loc,
                   std::vector<dealii::VectorizedArray<double> > & source_term,
                   std::vector<dealii::VectorizedArray<double> > & gamma) const;

	// ================================================================
	// Model constants specific to this subclass
	// ================================================================

  double m0 =userInputs.get_model_constant_double("m0");
  double MnV = userInputs.get_model_constant_double("MnV");
  double KnV = userInputs.get_model_constant_double("KnV");
  double alpha = userInputs.get_model_constant_double("alpha");
  double fsbs = userInputs.get_model_constant_double("fsbs");
  double A = userInputs.get_model_constant_double("A");
  double B = userInputs.get_model_constant_double("B");
  double globalSeedingTime = userInputs.get_model_constant_double("globalSeedingTime");
  double interface_coeff = std::sqrt(2.0*KnV/m0);
 
	// ================================================================

  //Dislocation densities
  double rho_i[160] = {1.01299, 1.20493, 0.92342, 0.64405,
    1.13669, 1.09403, 0.68884, 0.49050,
    0.78054, 1.23479, 1.06418, 0.92982,
    0.27084, 1.01513, 0.84878, 1.24972,
    1.07911, 1.03219, 1.23905, 1.22412,
    1.24972, 0.60566, 0.35615, 0.51183,
    1.23905, 1.21986, 1.20067, 0.42226,
    1.10470, 1.22626, 0.88077, 0.78054,
    1.12176, 1.10043, 1.25398, 1.23265,
    1.24972, 1.24972, 1.24972, 1.16654,
    1.23479, 1.05991, 1.24332, 0.90210,
    0.41799, 1.13669, 1.25398, 1.25398,
    1.23479, 1.18360, 0.48197, 1.13029,
    1.18574, 0.99593, 1.23052, 1.20280,
    0.72722, 1.05991, 0.59074, 1.10043,
    1.25185, 0.99807, 1.24972, 1.24758,
    1.24758, 1.16015, 1.12816, 1.07697,
    1.20920, 1.17081, 1.10256, 1.20280,
    1.19853, 1.25398, 1.17507, 0.72936,
    0.42013, 1.14948, 0.99380, 0.95541,
    1.22626, 1.22839, 0.70590, 0.80187,
    0.98954, 0.55235, 1.23265, 1.03432,
    0.49264, 0.38387, 1.14095, 1.14735,
    0.63126, 1.05565, 1.17721, 1.22839,
    1.25398, 0.70376, 1.23692, 1.24118,
    1.23905, 1.16228, 0.62059, 1.24972,
    0.21966, 1.12602, 0.36468, 0.94049,
    1.14308, 0.82959, 1.10256, 1.16015,
    0.61206, 1.17081, 1.25185, 0.36041,
    0.84665, 1.00233, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000};
};
