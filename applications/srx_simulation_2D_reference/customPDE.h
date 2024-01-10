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
  double rho_i[142] = {0.99511, 1.18365, 0.90712, 0.63268,
    1.11661, 1.07471, 0.67667, 0.48184,
    0.76675, 1.21298, 1.04538, 0.91340,
    0.26606, 0.99720, 0.83379, 1.22765,
    1.06005, 1.01396, 1.21717, 1.20251,
    1.22765, 0.59497, 0.34986, 0.50279,
    1.21717, 1.19832, 1.17946, 0.41480,
    1.08519, 1.20460, 0.86522, 0.76675,
    1.10195, 1.08100, 1.23184, 1.21089,
    1.22765, 1.22765, 1.22765, 1.14594,
    1.21298, 1.04119, 1.22136, 0.88617,
    0.41061, 1.11661, 1.23184, 1.23184,
    1.21298, 1.16270, 0.47346, 1.11033,
    1.16480, 0.97835, 1.20879, 1.18156,
    0.71438, 1.04119, 0.58030, 1.08100,
    1.22974, 0.98044, 1.22765, 1.22555,
    1.22555, 1.13966, 1.10823, 1.05795,
    1.18784, 1.15013, 1.08309, 1.18156,
    1.17737, 1.23184, 1.15432, 0.71648,
    0.41271, 1.12918, 0.97625, 0.93854,
    1.20460, 1.20670, 0.69343, 0.78770,
    0.97206, 0.54259, 1.21089, 1.01605,
    0.48394, 0.37709, 1.12080, 1.12709,
    0.62011, 1.03700, 1.15642, 1.20670,
    1.23184, 0.69134, 1.21508, 1.21927,
    1.21717, 1.14175, 0.60963, 1.22765,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000};
};
