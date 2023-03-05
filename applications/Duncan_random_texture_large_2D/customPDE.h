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
  double rho_i[141] = {0.74379, 0.97468, 0.95295, 0.94644,
    0.86263, 0.75147, 0.74440, 0.74710,
    1.11225, 0.89674, 1.00449, 1.09348,
    1.06179, 1.11087, 0.98390, 1.08104,
    0.99111, 1.10806, 0.97158, 1.71975,
    0.75387, 0.86920, 0.95554, 0.96139,
    0.97955, 0.75236, 1.06638, 0.74454,
    1.65525, 1.11892, 1.54539, 0.95616,
    0.89082, 0.74882, 0.75473, 0.74462,
    0.74412, 0.75278, 1.75593, 0.74409,
    1.77886, 0.97149, 0.75363, 1.61215,
    0.74882, 1.10649, 0.89006, 0.75405,
    0.86556, 1.83913, 0.79553, 0.90953,
    0.75431, 0.94140, 0.89288, 0.75396,
    1.07205, 1.73474, 0.74398, 1.07921,
    1.05769, 0.96830, 1.10695, 0.89845,
    1.08283, 0.93568, 0.76077, 0.75492,
    1.12178, 1.09644, 1.01612, 0.94505,
    0.89804, 1.09060, 1.10449, 0.74400,
    0.81443, 0.95761, 1.11311, 1.11660,
    1.84322, 1.53383, 1.10812, 1.10445,
    0.76093, 0.74432, 0.74269, 0.97696,
    0.74360, 0.75321, 0.75402, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000, 0.00000, 0.00000, 0.00000,
    0.00000};

};
